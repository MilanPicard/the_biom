from dash import Dash, dcc, html, Input, Output,callback
import dash_bootstrap_components as dbc


def menu(diseases,stages,diseases_cmap,stages_cmap):
    colors_str = dict([(d,f"rgba({','.join([str(i*255) for i in diseases_cmap[d][:3]])},{diseases_cmap[d][3]})") for d in diseases])
    disease_checklist = dcc.Checklist(
        [{"label":html.Div(d,style={
            'color':colors_str[d]
            }),"value":d,
        #   
            } for d in diseases],
        diseases,id="disease_filter",labelStyle={"display": "flex", "align-items": "center"})
    stage_checklist = dcc.Checklist(stages,stages,id="stage_filter")
    
    return [html.H1("menu"),
            dbc.Accordion([
                dbc.AccordionItem([
                    html.P("Diseases"),
                    disease_checklist
                ],title="disease_filter"),
                dbc.AccordionItem([
                    html.P("Stages"),
                    stage_checklist
                ],title="stage_filter"),
            ]
            )]