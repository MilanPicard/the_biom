import dash
import os
import sys
import dash_bootstrap_components as dbc
import data_manager
import logging
from dash import Output, Input, State
from dash import dcc
from flask import send_from_directory
import controller  # Ensure all callbacks are registered
controller.Controller.declare_callback()  # Register all Dash callbacks at top level
logging.getLogger('werkzeug').setLevel(logging.ERROR)

if "THE_BIOM_MODE" in os.environ and os.environ["THE_BIOM_MODE"]=="PROD":
    from cel_app import celery_app
    background_callback_manager = dash.CeleryManager(celery_app)
elif "THE_BIOM_MODE" in os.environ and os.environ["THE_BIOM_MODE"]=="DEV_DISK":
    import diskcache
    cache = diskcache.Cache("./cache")
    background_callback_manager = dash.DiskcacheManager(cache)
else:
    background_callback_manager = None
app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.BOOTSTRAP],
                assets_ignore="lib/.*",
                use_pages=True,
                suppress_callback_exceptions=True,
                background_callback_manager=background_callback_manager)

about_page = dash.page_registry["pages.about"]
main_page = dash.page_registry["pages.main_app"]
app.title = "THe Biom"
mouse_up_event = {"event":"mouseup","props":["target","buttons","offsetX","offsetY","type"]}
mouse_down_event = {"event":"mousedown","props":["target","buttons","offsetX","offsetY","type"]}
mouse_out_event = {"event":"mouseout","props":["buttons","offsetX","offsetY","type"]}
mouse_move_event = {"event":"mousemove","props":["buttons","offsetX","offsetY","target","type"]}

# Client-side callback for clipboard functionality
app.clientside_callback(
    """
    function(n_clicks, text_to_copy) {
        console.log('Clientside callback triggered:', {n_clicks: n_clicks, text_to_copy: text_to_copy});
        
        // Only proceed if we have a click and text_to_copy
        if (!n_clicks || n_clicks === 0) {
            console.log('No click detected, returning');
            return window.dash_clientside.no_update;
        }
        
        // Ensure we have a valid text_to_copy (even if it's just a space)
        var textToCopy = text_to_copy || " ";
        console.log('Text to copy:', textToCopy, 'Length:', textToCopy.length);
        
        // Use the modern Clipboard API
        if (navigator.clipboard && window.isSecureContext) {
            console.log('Using modern Clipboard API');
            navigator.clipboard.writeText(textToCopy).then(function() {
                console.log('Text copied to clipboard successfully:', textToCopy);
            }).catch(function(err) {
                console.error('Failed to copy text: ', err);
                // Fallback to older method
                fallbackCopyTextToClipboard(textToCopy);
            });
        } else {
            // Fallback for older browsers or non-secure contexts
            console.log('Using fallback clipboard method');
            fallbackCopyTextToClipboard(textToCopy);
        }
        
        return window.dash_clientside.no_update;
    }
    
    function fallbackCopyTextToClipboard(text) {
        console.log('Using fallback copy method for text:', text);
        var textArea = document.createElement("textarea");
        textArea.value = text;
        textArea.style.top = "0";
        textArea.style.left = "0";
        textArea.style.position = "fixed";
        document.body.appendChild(textArea);
        textArea.focus();
        textArea.select();
        try {
            var successful = document.execCommand('copy');
            if (successful) {
                console.log('Text copied to clipboard using fallback method:', text);
            } else {
                console.error('Fallback copy command failed');
            }
        } catch (err) {
            console.error('Fallback copy failed: ', err);
        }
        document.body.removeChild(textArea);
    }
    """,
    Output('dummy_div', 'children'),
    Input({"type": "legend-action", "action": "copy"}, 'n_clicks'),
    Input('clipboard_text_store', 'data'),
    prevent_initial_call=True
)

# Client-side callback for gProfiler functionality
app.clientside_callback(
    """
    function(n_clicks, gprofiler_url) {
        // Store the last processed n_clicks to prevent duplicate calls
        if (!window.lastGProfilerClicks) {
            window.lastGProfilerClicks = 0;
        }
        
        if (n_clicks && n_clicks > window.lastGProfilerClicks && gprofiler_url) {
            window.lastGProfilerClicks = n_clicks;
            
            // Open gProfiler URL in new tab
            window.open(gprofiler_url, '_blank');
            console.log('Opening gProfiler URL:', gprofiler_url);
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('dummy_div', 'children', allow_duplicate=True),
    Input({"type": "legend-action", "action": "gprofiler"}, 'n_clicks'),
    Input('gprofiler_url_store', 'data'),
    prevent_initial_call=True
)

# Custom route to serve files from the data directory
@app.server.route('/download/data/<path:filename>')
def download_data(filename):
    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
    return send_from_directory(data_dir, filename, as_attachment=True)

app.layout = dbc.Container([
    dcc.Store(id="pathway_to_highlight", data=None),
    dbc.Row([
        dbc.Col([
            dash.html.A(
                dash.html.H1("THe Biom", id="title", className="logo-title"),
                href=main_page["relative_path"],
                className="logo-link",
                id="logo_link"
            ),
            dash.html.Span("|", className="header-delimiter"),
            dash.dcc.Link(
                dash.html.Span("About", id="about_link_click", n_clicks=0),
                href=about_page["relative_path"],
                className="about-link",
                id="about_link_nav",
                style={"marginLeft": "12px", "fontSize": "1.08rem", "fontWeight": 500, "color": "#1a237e", "padding": 0, "background": "none", "border": "none", "boxShadow": "none"}
            ),
            dash.html.Span("|", className="header-delimiter"),
            dash.dcc.Link(
                dash.html.Span("Download", id="download_link_click"),
                href="/download",
                className="about-link",
                id="download_link_nav",
                style={"marginLeft": "12px", "fontSize": "1.08rem", "fontWeight": 500, "color": "#1a237e", "padding": 0, "background": "none", "border": "none", "boxShadow": "none"}
            )
        ], width="auto", style={"display": "flex", "alignItems": "center", "justifyContent": "flex-start", "paddingLeft": "24px", "paddingRight": "40px"})
    ], style={"flexGrow":0, "flexShrink":1, "width": "100%", "marginBottom": "8px"}),
    dash.page_container
], fluid=True, style={"height":"100vh","display":"flex","flexDirection":"column"}, id="main_container")

# Add a clientside callback to highlight the pathway in the detail graph after selection
app.clientside_callback(
    """
    function(pathway_id, elements) {
        if(pathway_id && window.tapMultiSignPathway) {
            // Wait a bit to ensure the graph is rendered
            setTimeout(function() {
                window.tapMultiSignPathway(pathway_id);
            }, 300);
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("dummy_div", "children", allow_duplicate=True),
    Input("pathway_to_highlight", "data"),
    Input("detail_graph", "elements"),
    prevent_initial_call=True
)

# Update the clientside callback to use logo_click n_clicks
app.clientside_callback(
    """
    function(n_clicks, pathname) {
        if (n_clicks > 0 && pathname === '/') {
            window.location.reload();
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('dummy_div', 'children', allow_duplicate=True),
    Input('logo_click', 'n_clicks'),
    Input({'type': 'location', 'index': dash.dash._ID_LOCATION}, 'pathname'),
    prevent_initial_call=True
)

# Update the clientside callback to use about_link_click n_clicks
app.clientside_callback(
    """
    function(n_clicks, pathname) {
        if (n_clicks > 0 && pathname === '/about') {
            window.location.reload();
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('dummy_div', 'children', allow_duplicate=True),
    Input('about_link_click', 'n_clicks'),
    Input({'type': 'location', 'index': dash.dash._ID_LOCATION}, 'pathname'),
    prevent_initial_call=True
)

application = app.server
if __name__ == '__main__':
    import controller

    if len(sys.argv)>=4 and all([os.path.exists(i) for i in sys.argv[1:4]]):
        signatures = sys.argv[1]
        expressions = sys.argv[2]
        pathways = sys.argv[3]
    elif ("THE_BIOM_SIGNATURES" in os.environ and "THE_BIOM_EXPRESSIONS" in os.environ and "THE_BIOM_PATHWAYS" in os.environ):
        signatures = os.environ["THE_BIOM_SIGNATURES"].strip()
        expressions = os.environ["THE_BIOM_EXPRESSIONS"].strip()
        pathways = os.environ["THE_BIOM_PATHWAYS"].strip()
    else:
        signatures =os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","signatures","THe_Biom_DEV_dataset.csv")
        expressions = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","activations","fake_data.csv")
        pathways = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","pathways","fake_pathways")
    dm = data_manager.DataManager(signatures,expressions,pathways)
    ctrl = controller.Controller(dm)

    # app.run(debug=False,host="172.17.11.246")
    app.run(debug=False)