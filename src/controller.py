import re
import time
# from dash_extensions.enrich import Dash, html, Input, Output, State, callback,ctx,ALL
from dash import Dash, html, dcc, Input, Output, State, callback,ctx,ALL,MATCH
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import clientside_callback,ClientsideFunction, Input, Output
# from dash_extensions.enrich import clientside_callback,ClientsideFunction, Input, Output
from dash import Patch
import dash
import numpy as np
import plotly.graph_objs as go
import pandas as pd
import plotly.express as px
import plotly
import plotly.subplots
import detail_box_plot
import overview as ov
import detail_graph
import utils
import stats
from  data_manager import DataManager
import cProfile,pstats,io
from dash.dependencies import ALL
from detail_graph import get_multi_signature_genes, get_last_green_genes, update_node_highlight_styles, legend_highlight_stylesheet

# Feature flag to enable/disable the toast for large pathway selection
ENABLE_PATHWAY_TOAST = False

# Configuration for pathway-signature relationship
PATHWAY_SIGNATURE_COVERAGE_THRESHOLD = 0.5  # 50% of required signatures must be available

class Controller(object):
    _instance = None
    def __new__(cls,dm):
        if cls._instance is None:
            cls._instance = super(Controller, cls).__new__(cls)
            cls._instance.dm = dm
        return cls._instance
    @staticmethod
    def declare_callback(background=False):
        @callback(
                Output("genes_menu_select","options"),
                Output("pathways_menu_select","options"),
                Input("filters_dropdown","value")
        )
        def update_selectable_genes(selected_filter):
            return [{"label":"None","value":"None"}]+[{"label":f"{' '.join(j['GeneSymbolID'])}","title" : f"{j['counts']} signatures","value":i} for i,j in DataManager.get_instance().get_genes(selected_filter).items()],[{"label":"None","value":"None"}]+[{"label":j['PathwayDisplayName'],"title" : f"{j['counts']} genes","value":i} for i,j in DataManager.get_instance().get_pathways([],selected_filter).items()]
        @callback(
                Output('selected_genes_store','data'),
                Output('pathway_to_highlight', 'data'),
                Output('pathway_signature_toast_store', 'data'),
                Output('pending_pathway_store', 'data'),
                Output('previous_selection_store', 'data'),
                Output('pathways_menu_select', 'value'),
                Output('revert_in_progress_store', 'data'),
                Output('block_ui_store', 'data'),
                Input("genes_menu_select","value"),
                Input({"type":"selected_gene_button","gene":ALL},"n_clicks"),
                Input("detail_graph","tapNodeData"),
                Input("detail_graph","selectedNodeData"),
                Input("filters_dropdown","value"),
                Input("pathways_menu_select","value"),
                Input({"type":"selected_pathway_button","pathway":ALL},"n_clicks"),
                Input("mono_graph","selectedNodeData"),
                Input("mono_graph","tapNodeData"),
                Input({"type": "close_boxplot", "gene": ALL}, "n_clicks"),
                Input('pending_pathway_store', 'data'),
                Input('previous_selection_store', 'data'),
                Input('revert_in_progress_store', 'data'),
                Input('block_ui_store', 'data'),
                Input('proceed_yes_store', 'data'),
                State('selected_genes_store','data'),
                State('pending_pathway_store', 'data'),
                State('previous_selection_store', 'data'),
                State('pathways_menu_select', 'value'),
                State('revert_in_progress_store', 'data'),
                State('block_ui_store', 'data'),
                State('proceed_yes_store', 'data'),
                prevent_initial_call=True
        )
        def add_remove_gene(menu_select,button,multip_tap,fromGraph,selected_filter,menu_pathway,pathway_button,fromMonoGraph,monotap,close_clicks,pending_pathway_input,previous_selection_input,revert_in_progress,block_ui,proceed_yes_pathway,current,pending_pathway_state,previous_selection_state,current_dropdown_value,revert_in_progress_state,block_ui_state,proceed_yes_state):
            import copy
            # Prevent interference from signature selection
            if (ctx.triggered_id == "detail_graph" and 
                "detail_graph.selectedNodeData" in ctx.triggered_prop_ids and 
                fromGraph is not None and 
                len(fromGraph) > 0):
                # Check if this is triggered by signature selection (nodes with "_" in ID are signatures)
                signature_nodes = [node for node in fromGraph if "_" in node.get("id", "")]
                if len(signature_nodes) > 0:
                    raise dash.exceptions.PreventUpdate()
            
            if revert_in_progress or revert_in_progress_state:
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, False
            if current is None:
                current = {"selected": [], "from_pathways": {"ids": [], "genes": [], "signatures": []}, "covered_signatures": [], "pathway_signature_map": {}}
            highlight_pathway_id = None
            toast_data = {}
            # Handle Yes button for large pathway ONLY if triggered by proceed_yes_store
            if proceed_yes_pathway and ctx.triggered_id == "proceed_yes_store":
                if proceed_yes_pathway not in current["from_pathways"]["ids"]:
                    to_add = DataManager.get_instance().get_genes(selected_filter, pathway=proceed_yes_pathway)
                    pathway_gene_ids = list(to_add.keys()) if isinstance(to_add, dict) else list(to_add)
                    if pathway_gene_ids:
                        # Get signatures associated with this pathway
                        pathway_signatures = DataManager.get_instance().get_pathway_signatures(
                            proceed_yes_pathway, selected_filter
                        )
                        
                        gene_inter = DataManager.get_instance().get_genes_intersections([], [], [], pathway_gene_ids, selected_filter)[1]
                        current["from_pathways"]["ids"].append(proceed_yes_pathway)
                        current["from_pathways"]["genes"].append(to_add)
                        if "signatures" not in current["from_pathways"]:
                            current["from_pathways"]["signatures"] = []
                        current["from_pathways"]["signatures"].append(pathway_signatures)
                        
                        # Update pathway signature map
                        if "pathway_signature_map" not in current:
                            current["pathway_signature_map"] = {}
                        current["pathway_signature_map"][proceed_yes_pathway] = pathway_signatures
                        
                        if gene_inter is not None:
                            current["covered_signatures"] = gene_inter["id"].explode().unique().tolist()
                        else:
                            current["covered_signatures"] = []
                        highlight_pathway_id = proceed_yes_pathway
                # Always reset proceed_yes_store to None after handling
                return current, highlight_pathway_id, {}, None, dash.no_update, current_dropdown_value, False, False
            if pending_pathway_input == 'REVERT' and previous_selection_input:
                return (
                    copy.deepcopy(previous_selection_input.get('selected_genes_store', {"selected": [], "from_pathways": {"ids": [], "genes": [], "signatures": []}, "covered_signatures": [], "pathway_signature_map": {}})),
                    None,
                    {},
                    None,
                    None,
                    previous_selection_input.get('pathways_menu_select', 'None'),
                    False,
                    False
                )
            previous_state = {
                'selected_genes_store': copy.deepcopy(current),
                'pathways_menu_select': current_dropdown_value
            }
            # Handle gene removal first
            if ctx.triggered and isinstance(ctx.triggered_id, dict):
                if ctx.triggered_id["type"] == "close_boxplot":
                    # Check if a close boxplot button was actually clicked
                    if close_clicks and any(click is not None for click in close_clicks):
                        gene_to_remove = ctx.triggered_id["gene"]
                        if gene_to_remove in current["selected"]:
                            current["selected"].remove(gene_to_remove)
                            if len(current["selected"]) == 0:
                                current["covered_signatures"] = []
                        return current, None, {}, None, previous_state, current_dropdown_value, False, False
                    else:
                        raise dash.exceptions.PreventUpdate()
                elif ctx.triggered_id["type"] == "selected_gene_button":
                    if ctx.triggered_id["gene"] in current["selected"]:
                        current["selected"].remove(ctx.triggered_id["gene"])
                        if len(current["selected"]) == 0:
                            current["covered_signatures"] = []
                    return current, None, {}, None, previous_state, current_dropdown_value, False, False
                elif ctx.triggered_id["type"] == "selected_pathway_button":
                    try:
                        index = current["from_pathways"]["ids"].index(ctx.triggered_id["pathway"])
                        current["from_pathways"]["ids"].remove(ctx.triggered_id["pathway"])
                        del current["from_pathways"]["genes"][index]
                        if "signatures" in current["from_pathways"] and len(current["from_pathways"]["signatures"]) > index:
                            del current["from_pathways"]["signatures"][index]
                        if "pathway_signature_map" in current and ctx.triggered_id["pathway"] in current["pathway_signature_map"]:
                            del current["pathway_signature_map"][ctx.triggered_id["pathway"]]
                        if len(current["selected"]) == 0 and len(current["from_pathways"]["ids"]) == 0:
                            current["covered_signatures"] = []
                    except ValueError:
                        pass
                    return current, None, {}, None, previous_state, current_dropdown_value, False, False
            if ctx.triggered_id == "filters_dropdown":
                current["selected"] = []
                current["from_pathways"]["ids"] = []
                current["from_pathways"]["genes"] = []
                current["from_pathways"]["signatures"] = []
                current["covered_signatures"] = []
                current["pathway_signature_map"] = {}
                return current, None, {}, None, previous_state, current_dropdown_value, False, False
            added = []
            match ctx.triggered_id:
                case "genes_menu_select":
                    if menu_select == "None":
                        current["selected"] = []
                        current["covered_signatures"] = []
                    elif menu_select is not None and menu_select.startswith("ENSG") and menu_select not in current["selected"]:
                        current["selected"].append(menu_select)
                        added = [menu_select]
                    else:
                        raise dash.exceptions.PreventUpdate()
                case "detail_graph":
                    if "detail_graph.tapNodeData" in ctx.triggered_prop_ids:
                        gene = multip_tap
                        if (not "is_pathways" in gene or not gene["is_pathways"]) and gene["id"].startswith("ENSG") and gene["id"] not in current["selected"]:
                            current["selected"].append(gene["id"])
                            added = [gene["id"]]
                        else:
                            raise dash.exceptions.PreventUpdate()
                    elif "detail_graph.selectedNodeData" in ctx.triggered_prop_ids and fromGraph is not None and len(fromGraph) > 0:
                        # Handle gene selection from detail graph (signature selection is already prevented above)
                        added_any = False
                        for gene in fromGraph:
                            if (not "is_pathways" in gene or not gene["is_pathways"]) and gene["id"].startswith("ENSG") and gene["id"] not in current["selected"]:
                                current["selected"].append(gene["id"])
                                added_any = True
                                added.append(gene["id"])
                        if not added_any:
                            raise dash.exceptions.PreventUpdate()
                    else:
                        raise dash.exceptions.PreventUpdate()
                case "mono_graph":
                    if fromMonoGraph is not None and len(fromMonoGraph) > 0:
                        for gene in fromMonoGraph:
                            if (not "is_pathways" in gene or not gene["is_pathways"]) and gene["id"].startswith("ENSG") and gene["id"] not in current["selected"]:
                                current["selected"].append(gene["id"])
                                added.append(gene["id"])
                    else:
                        raise dash.exceptions.PreventUpdate()
                case "pathways_menu_select":
                    if menu_pathway is not None and menu_pathway != "None":
                        if menu_pathway in current["from_pathways"]["ids"]:
                            index = current["from_pathways"]["ids"].index(menu_pathway)
                            current["from_pathways"]["ids"].remove(menu_pathway)
                            del current["from_pathways"]["genes"][index]
                            if len(current["selected"]) == 0 and len(current["from_pathways"]["ids"]) == 0:
                                current["covered_signatures"] = []
                            return current, None, {}, None, previous_state, current_dropdown_value, False, False
                        else:
                            to_add = DataManager.get_instance().get_genes(selected_filter, pathway=menu_pathway)
                            pathway_gene_ids = list(to_add.keys()) if isinstance(to_add, dict) else list(to_add)
                            if pathway_gene_ids:
                                gene_inter = DataManager.get_instance().get_genes_intersections([], [], [], pathway_gene_ids, selected_filter)[1]
                                n_signatures = len(gene_inter["id"].explode().unique().tolist()) if gene_inter is not None else 0
                                if ENABLE_PATHWAY_TOAST and n_signatures >= 5:
                                    toast_data = {
                                        "is_open": True,
                                        "message": [
                                            "This action will take about 1 minute, do you want to proceed?",
                                            html.Div([
                                                dbc.Button("Yes", id="proceed_yes_btn", color="primary", style={"marginRight": "1rem", "minWidth": "80px"}),
                                                dbc.Button("No", id="proceed_no_btn", color="secondary", style={"minWidth": "80px"})
                                            ], style={"marginTop": "1.2rem", "display": "flex", "justifyContent": "center", "gap": "0.5rem"})
                                        ]
                                    }
                                    return (
                                        dash.no_update, dash.no_update, toast_data, menu_pathway, previous_state, current_dropdown_value, False, False
                                    )
                                else:
                                    # Proceed as if user clicked Yes (add pathway immediately)
                                    if menu_pathway not in current["from_pathways"]["ids"]:
                                        # Get signatures associated with this pathway
                                        pathway_signatures = DataManager.get_instance().get_pathway_signatures(
                                            menu_pathway, selected_filter
                                        )
                                        
                                        current["from_pathways"]["ids"].append(menu_pathway)
                                        current["from_pathways"]["genes"].append(to_add)
                                        if "signatures" not in current["from_pathways"]:
                                            current["from_pathways"]["signatures"] = []
                                        current["from_pathways"]["signatures"].append(pathway_signatures)
                                        
                                        # Update pathway signature map
                                        if "pathway_signature_map" not in current:
                                            current["pathway_signature_map"] = {}
                                        current["pathway_signature_map"][menu_pathway] = pathway_signatures
                                        
                                        if gene_inter is not None:
                                            current["covered_signatures"] = gene_inter["id"].explode().unique().tolist()
                                        else:
                                            current["covered_signatures"] = []
                                        highlight_pathway_id = menu_pathway
                                    return current, highlight_pathway_id, {}, None, previous_state, current_dropdown_value, False, False
                    else:
                        raise dash.exceptions.PreventUpdate()
            if ctx.triggered_id == "genes_menu_select" and len(added) > 0:
                gene_inter = DataManager.get_instance().get_genes_intersections([], [], [], added, selected_filter)[1]
                if gene_inter is not None:
                    current["covered_signatures"] = gene_inter["id"].explode().unique().tolist()
                else:
                    current["covered_signatures"] = []
            return current, highlight_pathway_id, toast_data, None, previous_state, current_dropdown_value, False, False
        @callback(
            Output("selected_genes_div","children"),
            Input("selected_genes_store","data"),
            State("selected_genes_div","children"),
            prevent_initial_call=False
        )
        def update_genes_buttons(store_data,cur_children):
            return update_button(store_data["selected"], cur_children,attr="gene",get_symbol=DataManager.get_instance().get_symbol) 
        @callback(
            Output("selected_pathways_div","children"),
            Input("selected_genes_store","data"),
            State("selected_pathways_div","children"),
            prevent_initial_call=False
        )
        def update_pathways_buttons(store_data,cur_children):
            return update_button(store_data["from_pathways"]["ids"], cur_children,attr="pathway",get_symbol=DataManager.get_instance().get_pathway_label) 
        def update_button(store_data, cur_children,attr,get_symbol):
            patched = Patch()
            if cur_children is None:
                cur_children = []
            n = len(cur_children)
            if len(store_data)==0 and n>0:
                patched.clear()
            else:
                if n!=len(store_data):
                    color_scheme = detail_graph.get_color_scheme(store_data)
                    if(len(store_data)<len(cur_children)):
                        found=False
                        for i in range(len(store_data)):
                            if cur_children[i]["props"]["id"][attr]!=store_data[i]:
                                del patched[i]
                                del cur_children[i]
                                found = True
                                break
                        if not found:
                            del patched[len(cur_children)-1]
                            del cur_children[-1]
                        n-=1
                    else:
                        for i in range(len(cur_children),len(store_data)):
                            if attr == "gene":
                                # For genes, get the symbol
                                try:
                                    symbol = get_symbol(store_data[i])
                                    btn = html.Button(symbol,id={"type":f"selected_{attr}_button",attr:store_data[i]},className="btn",style={})
                                except KeyError:
                                    # If it's not a gene ID, skip it
                                    continue
                            else:
                                # For pathways, get the label
                                btn = html.Button(get_symbol(store_data[i]),id={"type":f"selected_{attr}_button",attr:store_data[i]},className="btn",style={})
                            patched.append(btn)
                            n+=1
                    for i,g in enumerate(store_data):
                        if attr == "gene" and not g.startswith("ENSG"):
                            continue
                        tc = utils.get_text_color(color_scheme[i%len(color_scheme)])
                        patched[i]["props"]["style"]={
                            "--bs-btn-bg":color_scheme[i%len(color_scheme)],
                            "--bs-btn-hover-bg":color_scheme[i%len(color_scheme)],
                            "--bs-btn-hover-color":tc,
                            "--bs-btn-color":tc}
                else:
                    return patched
            return patched         


        @callback(
                Output('overview_graph','elements'),
                Output('overview_graph_layout','data'),
                Input("filters_dropdown","value"),
                State("overview_graph","elements")
                )
        def update_overview(
                            selected_filter,
                            cur_elems):
            overview_graph_layout = dash.no_update
            if any(["Filter" in c["data"] and c["data"]["Filter"]!=selected_filter for c in cur_elems]):
                cur_elems = ov.get_elements(DataManager.get_instance(),selected_filter=selected_filter)
                overview_graph_layout = {"rerun":True}
            return cur_elems,overview_graph_layout

        @callback(Output('data_overview_selected','value'),
                Input('overview_graph','selectedNodeData'), prevent_initial_call=True
                )
        def update_overview_selected_data(data):
            if data is None:
                return ""
            return ";".join([i["id"] for i in data])


        @callback(Output('detail_graph','elements'),
                Output('detail_graph','stylesheet'),
                Output('detail_graph','layout'),
                Output('detail_graph_pos','data'),
                Output('multi_legend',"data"),
                Output('multi_export',"data"),
                Input('overview_graph','selectedNodeData'),
                Input('selected_genes_store','data'),
                Input("fake_graph_size","data"),
                Input("filters_dropdown","value"),
                Input("show_pathways_store", "data"),
                State('detail_graph','elements'),
                State("detail_graph_pos","data"),            
                State('detail_graph','stylesheet'),
                State('legend_active_store', 'data'),
                    prevent_initial_call=True,
                        cancel=[
                Input('overview_graph','selectedNodeData'),
                Input("fake_graph_size","data"),
                Input("filters_dropdown","value"),
                                ],    background=background,
                )
        def display_detail_graph(
            overview_nodes,
            menu_genes,fake_graph_size,selected_filter,show_pathways,existing_elements,detail_pos_store,current_stylesheets,legend_active_data):
            selected = [n["id"] for n in overview_nodes] if overview_nodes is not None else []
            if (ctx.triggered_id == "selected_genes_store"
                and menu_genes is not None
                and "covered_signatures" in menu_genes
                and len(menu_genes["covered_signatures"]) > 0):
                if any([i not in selected for i in menu_genes["covered_signatures"]]):
                    raise dash.exceptions.PreventUpdate
            signatures = ";".join(selected) if overview_nodes is not None else None
            if ctx.triggered_id=="fake_graph_size" :
                if (fake_graph_size is None or "just_redraw" not in fake_graph_size or not fake_graph_size["just_redraw"]):
                    raise dash.exceptions.PreventUpdate()
                else:
                    return detail_graph.redraw(existing_elements,detail_pos_store,1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],current_stylesheets)
            if(all([
                  signatures is None or signatures ==""])):
                return [],[],{"name":"preset"},{},dash.no_update,dash.no_update            
            diseases, comparisons, signatures, genes_set = get_detail_subset(None, [], signatures, menu_genes,selected_filter)
            # Determine highlight_gene_ids from legend state
            highlight_gene_ids = []
            if legend_active_data and 'active_ids' in legend_active_data:
                for active_id in legend_active_data['active_ids']:
                    import ast
                    item_dict = ast.literal_eval(active_id)
                    if item_dict.get('category') == 'genes' and item_dict.get('name') == 'genes in multiple signatures':
                        highlight_gene_ids = get_last_green_genes()
            if len(diseases)!=0 or len(signatures)!=0:
                r = detail_graph.display_detail_graph([],signatures,genes_set,existing_elements,detail_pos_store if detail_pos_store is not None else dict(),1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],selected_filter,comparisons,show_pathways=show_pathways,highlight_gene_ids=highlight_gene_ids)
                return r
            else:
                return existing_elements,[],{"name":"preset"},{},dash.no_update,dash.no_update

        def get_detail_subset(diseases, comparisons, signatures, menu_genes,selected_filter):
            if diseases is None:
                diseases = []
            if signatures is None:
                signatures = []
            if len(signatures) > 0:
                if isinstance(signatures, str):
                    signatures = signatures.split(";")
                # Only process signature IDs (those containing "_")
                for s in signatures:
                    if "_" in s:  # Only process signature IDs
                        cancer, comp, fil = s.split("_")
                        if cancer not in diseases:
                            diseases.append(cancer)
                        if comp not in comparisons:
                            comparisons.append(comp)
                        
            genes_set = set()
            for p in menu_genes["from_pathways"]["genes"]:
                genes_set.update(p)
            genes_set = menu_genes["selected"] + sorted(list(genes_set.difference(menu_genes["selected"])))
            
            diseases = list(filter(lambda a: len(a) > 0, diseases))
            if isinstance(signatures, str):
                signatures = list(filter(lambda a: len(a) > 0, signatures.split(";")))
            return diseases, comparisons, signatures, genes_set
        @callback(
                Output('mono_graph','elements'),
                Output('mono_graph','stylesheet'),
                Output('mono_graph','layout'),
                Output('mono_graph_pos','data'),
                Output('mono_legend',"data"),
                Output('mono_export',"data"),
                Output('mono_tab','label'),
                Input("overview_graph","tapNodeData"), 
                Input("fake_graph_size","data"),
                Input("filters_dropdown","value"),
                State('selected_genes_store','data'),
                State('mono_graph','elements'),
                State("mono_graph_pos","data"),            
                State('mono_graph','stylesheet'),
                    prevent_initial_call=True
                )
        def display_mono_graph(tapNodeData,fake_graph_size,selected_filter,
                               menu_genes,existing_elements,detail_pos_store,current_stylesheets):
            d = None
            c = None
            s = None
            if tapNodeData is not  None and (ctx.triggered_id=="overview_graph"):
                d = tapNodeData["Cancer"]
                c = tapNodeData['Comparison']
                s = tapNodeData["id"]
            if d is None:
                raise dash.exceptions.PreventUpdate()
            genes_set = set()
            for p in menu_genes["from_pathways"]["genes"]:
                genes_set.update(p)
            genes_set = menu_genes["selected"] + sorted(list(genes_set.difference(menu_genes["selected"])))
            r = detail_graph.display_detail_graph([d],[s],list(genes_set),existing_elements,detail_pos_store if detail_pos_store is not None else dict(),1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],selected_filter,[c],show_pathways=True)
            return *r,s

        # @callback(Output("data_gene_detail_selected","value"),
        #           Input("detail_graph","selectedNodeData"))
        # def data_gene_detail_selected(nodes):
        #     if nodes is None:
        #         return ""
        #     return ";".join([i["id"] for i in nodes])

        @callback(
            Output("box_plots_container", "children"),
            Output("overview_graph","stylesheet"),
            Output("box_plots_to_style","data"),
            Output("box_plots_stats","data"),
            Input("selected_genes_store","data"),  # Only depend on selected_genes_store
            Input("box_categories","data"),
            Input("filters_dropdown","value"),
            Input({"type": "close_boxplot", "gene": ALL}, "n_clicks"),
            Input("overview_graph", "selectedNodeData"),  # Add this line
            State("overview_graph","elements"),
            State("overview_graph","stylesheet"),
            prevent_initial_call=False
        )
        def update_box_plot(
            menu_selected, selected_boxcategories,
            selected_filter, close_clicks, overview_selected_nodes,
            overview_elements, overview_stylesheets):
            ctx = dash.callback_context
            
            # Extract selected cancer(s) from overview graph selection
            selected_cancers = []
            if overview_selected_nodes:
                for node in overview_selected_nodes:
                    if "Cancer" in node:
                        selected_cancers.append(node["Cancer"])
            selected_cancers = list(set(selected_cancers))
                    
            # Initialize box_categories at the start
            box_categories = []
            
            # Handle close button clicks
            if ctx.triggered and "close_boxplot" in str(ctx.triggered_id):
                gene_to_remove = ctx.triggered_id["gene"]
                if menu_selected is not None and "selected" in menu_selected:
                    items = [g for g in menu_selected["selected"] if g != gene_to_remove]
                    if not items:
                        return [], overview_stylesheets, {"categories":[],"genes":[]}, {"stats":[]}
                else:
                    return [], overview_stylesheets, {"categories":[],"genes":[]}, {"stats":[]}

            # Get genes from menu selection
            items = []
            if menu_selected is not None and "selected" in menu_selected:
                items = menu_selected["selected"].copy()
                if "from_pathways" in menu_selected and "genes" in menu_selected["from_pathways"]:
                    genes_set = set()
                    for p in menu_selected["from_pathways"]["genes"]:
                        genes_set.update(p)
                    pathway_genes = sorted(list(genes_set.difference(menu_selected["selected"])))
                    items.extend(pathway_genes)
            
            stylesheets = overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(DataManager.get_instance())
            stylesheets = [s for s in stylesheets if not(s["selector"].startswith("edge#"))]
            
            if len(items) > 0:
                gene_items = [item for item in items if isinstance(item, str) and item.startswith("ENSG")]
                signature_items = [item for item in items if isinstance(item, str) and "_" in item]
                
                diseases_detail, comparisons, signatures, genes_set = get_detail_subset(
                    None, [], signature_items, menu_selected, selected_filter
                )
                
                if gene_items:
                    gene_diseases = DataManager.get_instance().get_diseases_from_genes(gene_items)
                    diseases_detail.extend(gene_diseases)
                    diseases_detail = list(set(diseases_detail))
                
                box_categories_tohighlight = dict({i: set() for i in items})
                gene_inter = DataManager.get_instance().get_genes_intersections(diseases_detail,comparisons,signatures,genes_set,selected_filter)[1]
                if gene_inter is None:
                    return [], stylesheets, {"categories":[],"genes":[]}, {"stats":[]}

                diseases = []
                if overview_elements is not None:
                    for elem in overview_elements:
                        if "data" in elem and "elems" in elem["data"]:
                            for e in elem["data"]["elems"]:
                                if isinstance(e, str) and "_" in e:
                                    disease = e.split("_")[0]
                                    if disease not in diseases:
                                        diseases.append(disease)

                # Get color maps for diseases and comparisons
                disease_cmap = DataManager.get_instance().get_disease_cmap()
                comparison_cmap = DataManager.get_instance().get_comparison_cmap()

                # Create a div for each gene's boxplot
                boxplot_divs = []
                for gene in items:
                    if not isinstance(gene, str) or not gene.startswith("ENSG"):
                        continue
                        
                    selected_patient_and_genes = DataManager.get_instance().get_activations(
                        [gene],
                        selected_cancers if selected_cancers else diseases if len(selected_boxcategories["diseases"])==0 else selected_boxcategories["diseases"],
                        selected_boxcategories["comparisons"]
                    ).sort_values(["box_category"])
                    
                    if selected_patient_and_genes.empty:
                        continue
                        
                    current_box_categories = sorted(pd.unique(selected_patient_and_genes["box_category"]).tolist())
                    box_categories = list(set(box_categories + current_box_categories))
                    symbol = " ".join(DataManager.get_instance().get_symbol(gene))
                    selected_patient_and_genes = selected_patient_and_genes.rename(columns={gene: symbol})
                    
                    # Create the boxplot with color coding
                    fig = go.Figure()
                    for category in current_box_categories:
                        disease, stage = category.split("_")
                        category_data = selected_patient_and_genes[selected_patient_and_genes["box_category"] == category]
                        
                        # Determine color based on disease
                        color = disease_cmap.get(disease, "gray")
                        if isinstance(color, tuple):  # Convert matplotlib RGB tuple to plotly rgba
                            # Convert from 0-1 range to 0-255 range and create rgba string
                            r, g, b = [int(x * 255) for x in color[:3]]
                            rgba_color = f"rgba({r},{g},{b},0.2)"
                            # For the line color, use solid rgb
                            line_color = f"rgb({r},{g},{b})"
                        else:
                            # Default gray colors
                            rgba_color = "rgba(128,128,128,0.2)"
                            line_color = "rgb(128,128,128)"
                        
                        fig.add_trace(go.Box(
                            y=category_data[symbol],
                            x=[category] * len(category_data),
                            name=category,
                            showlegend=False,
                            fillcolor=rgba_color,
                            line=dict(color=line_color)
                        ))
                    
                    fig.update_layout(
                        title=symbol,
                        showlegend=False,
                        margin=dict(t=30, l=10, r=10, b=10),
                        height=200
                    )
                    
                    boxplot_div = html.Div([
                        html.Div([
                            html.Button(
                                "×",
                                id={"type": "close_boxplot", "gene": gene},
                                className="close-button"),
                            dcc.Graph(
                                figure=fig,
                                config={"displayModeBar": False},
                                style={"height": "100%"}
                            )
                        ], className="boxplot-wrapper")
                    ], style={"width": "100%", "marginBottom": "10px"})
                    
                    boxplot_divs.append(boxplot_div)
                
                return boxplot_divs, stylesheets, {"categories": box_categories, "genes": items}, {"stats": []}
            
            return [], stylesheets, {"categories":[],"genes":[]}, {"stats":[]}

        @callback(
            Output("data_gene_menu_selected","value"),
            Input("genes_menu_select","value")
        )
        def update_data_gene_menu_selected(v):
            return v


        clientside_callback(
            ClientsideFunction(
                namespace='clientside',
                function_name='activate_tooltip'
            ),
            Output("detail_graph_tooltip","children"),
            Output("detail_graph_tooltip","show"),
            Output("detail_graph_tooltip","bbox"),
            Output("detail_graph_tooltip","direction"),
            Output("detail_graph_tooltip","className"),
            Input("detail_graph","mouseoverNodeData"),
            Input("detail_graph","mouseoverEdgeData"),
            Input("detail_graph","selectedEdgeData"),
            Input("fake_graph_size","data"),
            State("detail_graph","elements"),
            State("detail_graph","extent"),
            State("detail_graph","stylesheet"),
            State("detail_graph_pos","data"),
            State("detail_graph_tooltip","children"),
            State("detail_graph_tooltip","show"),
            State("detail_graph_tooltip","bbox"),
            State("detail_graph_tooltip","direction"),
            prevent_initial_call=True
            )
        clientside_callback(
            ClientsideFunction(
                namespace='clientside',
                function_name='activate_tooltip'
            ),
            Output("mono_graph_tooltip","children"),
            Output("mono_graph_tooltip","show"),
            Output("mono_graph_tooltip","bbox"),
            Output("mono_graph_tooltip","direction"),
            Output("mono_graph_tooltip","className"),
            Input("mono_graph","mouseoverNodeData"),
            Input("mono_graph","mouseoverEdgeData"),
            Input("mono_graph","selectedEdgeData"),
            Input("fake_graph_size","data"),
            State("mono_graph","elements"),
            State("mono_graph","extent"),
            State("mono_graph","stylesheet"),
            State("mono_graph_pos","data"),
            State("mono_graph_tooltip","children"),
            State("mono_graph_tooltip","show"),
            State("mono_graph_tooltip","bbox"),
            State("mono_graph_tooltip","direction"),
            prevent_initial_call=True
            )
        clientside_callback(
            ClientsideFunction(
                namespace='clientside',
                function_name='activate_tooltip'
            ),
            Output("overview_graph_tooltip","children"),
            Output("overview_graph_tooltip","show"),
            Output("overview_graph_tooltip","bbox"),
            Output("overview_graph_tooltip","direction"),
            Output("overview_graph_tooltip","className"),
            Input("overview_graph","mouseoverNodeData"),
            Input("overview_graph","mouseoverEdgeData"),
            Input("overview_graph","selectedEdgeData"),
            Input("fake_graph_size","data"),
            State("overview_graph","elements"),
            State("overview_graph","extent"),
            State("overview_graph","stylesheet"),
            State("detail_graph_pos","data"),
            State("overview_graph_tooltip","children"),
            State("overview_graph_tooltip","show"),
            State("overview_graph_tooltip","bbox"),
            State("overview_graph_tooltip","direction"),
            prevent_initial_call=True
            )

        # clientside_callback(ClientsideFunction(
        #     namespace="tooltip",
        #     function_name="set_tooltip"),
        #     Output("dummy_div","children"),
        #     Input("detail_graph","mouseoverNodeData"),prevent_initial_call=True
        # )

        clientside_callback(
            """
    function update_width_start(n1,e_width,state,detail_graph_pos,dw,fake_graph_size){
        var graph_width = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientWidth:document.getElementById("mono_graph").clientWidth;
        var graph_height = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientHeight:document.getElementById("mono_graph").clientHeight;

        var e = e_width;
        if(dash_clientside.callback_context.triggered_id===undefined){
            var elem = document.getElementById('detail_graph');
            var span = document.getElementById('detail_resize_span');

            return  Object.assign({},fake_graph_size,{'width':graph_width,'height':graph_height,'width_span':span.clientWidth,"AR":graph_width/graph_height});
        }
        if(e["type"]=="click" && !e["isTrusted"]){
            if(fake_graph_size["width"]===undefined || fake_graph_size["height"]===undefined || Math.abs(fake_graph_size["width"]-graph_width)/graph_width>0.05 ||Math.abs(fake_graph_size["height"]-graph_height)/graph_height>0.05){
                fake_graph_size["width"]=graph_width;
                fake_graph_size["height"]=graph_height;
                fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                fake_graph_size["just_redraw"]=true;
                
                // Update layout directly instead of calling layout_overview
                var overview_graph = document.getElementById("overview_graph");
                if (overview_graph && overview_graph['_cyreg'] && overview_graph['_cyreg']["cy"]) {
                    overview_graph['_cyreg']["cy"].layout({
                        name: "cola",
                        nodeDimensionsIncludeLabels: true,
                        animate: false
                    }).run();
                }

                return fake_graph_size;
            }
            else
                return dash_clientside.no_update;
        }
        if(e["type"]=="mousedown"){
            let in_tooltip = false;
            for(let t of document.querySelectorAll(".dcc-tooltip-bounding_box")){
                in_tooltip = in_tooltip || t.contains(e.target);
            }

            if(false && e.target.classList.contains("resize_span") && !in_tooltip){
                if(dash_clientside.callback_context.triggered_id==="move_in_ov"){
                    if(!state["width"]["is_resizing"]){
                        state["width"]["is_resizing"]=true;
                        state["width"]["original_width"]=dw;
                        dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
                        document.querySelectorAll("#overview_col canvas, #detail_col canvas").forEach(c => c.style.cursor="w-resize");
                        return dash_clientside.no_update;
                    }
                }
            }else{
                if(!in_tooltip){
                    var possibleTargets = ["overview_col","detail_col"];
                    if(false && ((dash_clientside.callback_context.triggered_id==="full_col" && "second_row_div" == e.target.id) || (possibleTargets.includes(e.target.id) && (e.target.clientHeight-e.offsetY<50)&& dash_clientside.callback_context.triggered_id==="move_in_ov"))){
                        if(!state["height"]["is_resizing"]){
                            state["height"]["is_resizing"]=true;
                            dash_clientside.set_props("resize_state", {"data":{"width":{"is_resizing":false},"height":{"is_resizing":true}}});
                            return dash_clientside.no_update;
                        }
                    }else{
                        console.log(possibleTargets.includes(e.target.id),e.target.id,dash_clientside.callback_context.triggered_id==="full_col",dash_clientside.callback_context.triggered_id)
                    }
                }
            }
            return dash_clientside.no_update;
        }
        if(e["type"]=="mouseup"){
            if(state["width"]["is_resizing"]){
                state["width"]["is_resizing"]=false;
                dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
                document.querySelectorAll("#overview_col canvas, #detail_col canvas").forEach(c => c.style.cursor="auto");
                layout_overview();

                fake_graph_size["width"]=graph_width;
                fake_graph_size["height"]=graph_height;
                fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                return fake_graph_size;
            }
            if(state["height"]["is_resizing"]){
                state["height"]["is_resizing"]=false;
                dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
    //        document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
                layout_overview();

                fake_graph_size["width"]=graph_width;
                fake_graph_size["height"]=graph_height;
                fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                return fake_graph_size;
            }
            return dash_clientside.no_update;
        }
        if(e["type"]=="mousemove" && state["width"]["is_resizing"] && e["buttons"]==0){
            console.log("triggered event",document.getElementById("overview_col").dispatchEvent(new MouseEvent("mouseup")));
        }
        if(e["type"]=="mousemove" && state["height"]["is_resizing"] && e["buttons"]==0){
            state["height"]["is_resizing"]=false;
            dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
    //        document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
            layout_overview();

            fake_graph_size["width"]=graph_width;
            fake_graph_size["height"]=graph_height;
            fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
            return fake_graph_size;
        }
        return dash_clientside.no_update;
    }
            """,
            Output('fake_graph_size','data'),
            Input("move_in_ov","n_events"),
            State("move_in_ov","event"),
            State("resize_state","data"),
            State("detail_graph_pos","data"),
            State("detail_col","width"),
            State('fake_graph_size','data'),
            prevent_initial_call=False
                    )

        clientside_callback(
            """
    function update_width(n,e_width,ow,dw,resize_state,fake_graph_size){
        //console.time("update_width");
        var e = e_width;
        var graph_width = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientWidth:document.getElementById("mono_graph").clientWidth;
        var graph_height = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientHeight:document.getElementById("mono_graph").clientHeight;

        function find_row(elem){
            var p =  elem;
            while(!p.classList.contains("row"))
                p = p.parentNode;
            return p;
        }
        
        // Handle window resize
        window.addEventListener('resize', function() {
            if (fake_graph_size === null) fake_graph_size = {};
            fake_graph_size["width"] = graph_width;
            fake_graph_size["height"] = graph_height;
            fake_graph_size["AR"] = graph_width/graph_height;
            fake_graph_size["just_redraw"] = true;
            dash_clientside.set_props("fake_graph_size", fake_graph_size);
            
            // Update layout
            var overview_graph = document.getElementById("overview_graph");
            if (overview_graph && overview_graph['_cyreg'] && overview_graph['_cyreg']["cy"]) {
                overview_graph['_cyreg']["cy"].layout({
                    name: "cola",
                    nodeDimensionsIncludeLabels: true,
                    animate: false
                }).run();
            }
        });

        if(false && resize_state["width"]["is_resizing"]){
            var overviewCol = document.getElementById("overview_col");
            var detailCol = document.getElementById("detail_col");
            if(e["type"]=="mousemove" && e.target && e.target.tagName==="CANVAS" && (overviewCol.contains(e.target) || detailCol.contains(e.target))){
                if(overviewCol.contains(e.target) && ow>2){
                    var w = ow;
                    var width = find_row(e.target).clientWidth;
                    while(w>2 && e.offsetX<(w-1)*width/12){
                        w=w-1;
                    }
                    if(w!=ow){
                        dash_clientside.set_props("overview_col", {width: w});
                        dash_clientside.set_props("detail_col", {width: 12-w});
                        fake_graph_size["width"]=graph_width;
                        fake_graph_size["height"]=graph_height;
                        fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                        fake_graph_size["just_redraw"]=true;
                        dash_clientside.set_props("fake_graph_size", fake_graph_size);
                    }
                }else{
                    if(detailCol.contains(e.target) && dw>2){
                        var width = find_row(e.target).clientWidth;
                        var w = dw;
                        while(w>2 && e.target.clientWidth-e.offsetX<(w-1)*width/12){
                            w=w-1;
                        }
                        if(dw!=w){
                            dash_clientside.set_props("overview_col", {width: 12-w});
                            dash_clientside.set_props("detail_col", {width: w});
                            fake_graph_size["width"]=graph_width;
                            fake_graph_size["width_span"]=document.getElementById("detail_resize_span").clientWidth;
                            fake_graph_size["height"]=graph_height;
                            fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                            fake_graph_size["just_redraw"]=true;
                            dash_clientside.set_props("fake_graph_size", fake_graph_size);
                        }
                    }
                }
            }
        }

        if(false && resize_state["height"]["is_resizing"] && e["type"]=="mousemove"){
            var upHeight = document.getElementById("overview_col").parentNode.clientHeight;
            if(ajust_flex(e.offsetY+upHeight)){
                fake_graph_size["width"]=graph_width;
                fake_graph_size["height"]=graph_height;
                fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                fake_graph_size["just_redraw"]=true;
                dash_clientside.set_props("fake_graph_size", fake_graph_size);
            }
        }

        return dash_clientside.no_update;
    }
                """,
            Output("overview_col","width"),
            Output("detail_col","width"),
            Input("move_in_ov","n_events"),
            State("move_in_ov","event"),
            State("overview_col","width"),
            State("detail_col","width"),
            State("resize_state","data"),
            State("fake_graph_size","data"),
            
            prevent_initial_call=True
                    )

        @callback(
            Output("about_link","style"),
            Output("main_link","style"),
            Input(dash.dash._ID_LOCATION,"pathname"),
            State("about_link","style"),
            State("main_link","style"),
            prevent_initial_call=True
        )
        def update_link_style(pathname,about,main):
            about = dict(about) if about else {}
            main = dict(main) if main else {}
            match pathname:
                case "/":
                    about["display"]="inline"
                    main["display"]="none"
                case "/about":
                    main["display"]="inline"
                    about["display"]="none"
            return about,main

        clientside_callback(
            """function sign_clipboard(n,sign){
                if(n>0){
                    navigator.clipboard.writeText(sign[0]);
                }
                return dash_clientside.no_update;
            }""",
                Output("title","id",allow_duplicate=True),
                    Input({"type":"signature_clipboard","sign":ALL},"n_clicks"),
                    State({"type":"signature_clipboard","sign":ALL},"value"),
                        prevent_initial_call=True

        )

        clientside_callback(    ClientsideFunction(
                namespace='clientside',
                function_name='layout_overview'
            ),

                                    Output('overview_graph','layout'),
                                    Input('overview_graph_layout','data'),
        )

        clientside_callback( ClientsideFunction(
                namespace='clientside',
                function_name='on_box_plot_click'
            ),
                Output("title","id",allow_duplicate=True),
                Input("activation_boxplot","hoverData"),
                Input("activation_boxplot","restyleData"),
                Input("activation_boxplot","figure"),
                State("box_plots_stats","data"),
                prevent_initial_call=True,allow_duplicate=True

        )
        clientside_callback( ClientsideFunction(
                namespace='clientside',
                function_name='display_legend'
            ),
                Output("title","id",allow_duplicate=True),
                Input('multi_legend',"data"),
                Input('mono_legend',"data"),
                        prevent_initial_call=True,allow_duplicate=True

        )
        clientside_callback( """
            function export_images(n_clicks,elem){
                switch(elem){
                    case "overview":
                        download_canvas_image(document.querySelector("#overview_graph canvas:nth-of-type(3)"),"TheBiom_overview.png");
                        break;
                    case "mono_graph":
                        download_canvas_image(document.querySelector("#mono_graph canvas:nth-of-type(3)"),"TheBiom_mono_signature.png",document.getElementById("mono_canvas"));
                        break;
                    case "detail_graph":
                        download_canvas_image(document.querySelector("#detail_graph canvas:nth-of-type(3)"),"TheBiom_multi_signature.png",document.getElementById("multi_canvas"));
                        break;
                    case "box":
                        const boxplotContainer = document.querySelector("#box_plots_container");
                        if (boxplotContainer && boxplotContainer.children.length > 0) {
                            // Get all plotly plots in the container
                            const plots = boxplotContainer.querySelectorAll(".js-plotly-plot");
                            if (plots.length > 0) {
                                // Export all boxplots as a combined image
                                download_all_boxplots(plots, "TheBiom_expression_boxplots");
                            }
                        }
                        break;
                }
                return dash_clientside.no_update;
            }
        """,
                Output("title","id",allow_duplicate=True),
                Input('export_image_btn',"n_clicks"),
                State('exportImage',"value"),
                        prevent_initial_call=True,allow_duplicate=True

        )
        clientside_callback( """
            function export_json(n_clicks,elem,multi_export,mono_export){
                switch(elem){
                    case "overview":
                        export_overview_data();
                        break;
                    case "mono_graph":
                    case "detail_graph":
                        let a = document.createElement("a");
                        let title = document.getElementById("mono_tab").querySelector("a").innerText;
                            
                        a.download= elem == "detail_graph"?"multi_signature_data.json":title+".json";
                        a.href="data:text/json;charset=utf-8,"+encodeURIComponent(JSON.stringify(elem == "detail_graph"?multi_export:mono_export));
                        a.click();

                        break;
                    case "box":
                        break;
                }
                return dash_clientside.no_update;
            }
        """,
                Output("title","id",allow_duplicate=True),
                Input('export_json_btn',"n_clicks"),
                State('exportImage',"value"),
                State('multi_export',"data"),
                State('mono_export',"data"),
                        prevent_initial_call=True
        )
        clientside_callback( """
            function tap_multi_node(tappedNodeData){
                let id = undefined;
                
                if(tappedNodeData !==undefined && tappedNodeData["classes"].includes("pathway")){
                    id = (tappedNodeData["data"]["id"]);
                }
                tapMultiSignPathway(id);
                console.log(id);
                return dash_clientside.no_update;
                                            
            }
        """,
                Output("title","id",allow_duplicate=True),
                Input('detail_graph',"tapNode"),
                        prevent_initial_call=True
        )

        clientside_callback( ClientsideFunction(
            namespace='clientside',
            function_name='set_min_height_box_plot'
            ),
                Output("title","id",allow_duplicate=True),
                Input('box_plots_to_style',"data"),
                        prevent_initial_call=True
        )        

        @callback(
            Output({"type":"signature_checkbox","disease":ALL,"comparison":ALL},"disabled"),
            
            Input("filters_dropdown","value"),
            State({"type":"signature_checkbox","disease":ALL,"comparison":ALL},"disabled"),
            State({"type":"signature_checkbox","disease":ALL,"comparison":ALL},"id")
            )
        def enable_existing_signatures_btn(selected_filter,button_set,ids):
            existing_signatures_list = DataManager.get_instance().get_signatures_id(selected_filter)
            existing_signatures = {}
            for s in existing_signatures_list:
                f = s.split("_")
                d = f[0]
                c = f[1]
                if d not in existing_signatures:
                    existing_signatures[d]=set([c])
                else:
                    existing_signatures[d].add(c)
            for i,disabled_state in  enumerate(button_set):
                d = ids[i]["disease"]
                c = ids[i]["comparison"]
                if(d not in existing_signatures):
                    # print(f"disabling all for {d}")
                    button_set[i]=True
                else:
                    if disabled_state==(c in existing_signatures[d]):
                        button_set[i]= not disabled_state                
            return button_set
        
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="update_selected_signatures_store"
        ),
            Output("selected_signatures","data"),
            Input({"type":"signature_checkbox","disease":ALL},"value"),
            State("overview_graph","selectedNodeData"),
            State({"type":"signature_checkbox","disease":ALL},"id"),
        )
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="update_overview_selection"
        ),
            Output({"type":"signature_checkbox","disease":MATCH,"comparison":MATCH},"id"),
            Input({"type":"signature_checkbox","disease":MATCH,"comparison":MATCH},"n_clicks"),
            Input({"type":"check_diseases","disease":MATCH},"n_clicks"),
            Input({"type":"uncheck_diseases","disease":MATCH},"n_clicks"),
            Input({"type":"check_comparisons","comparison":MATCH},"n_clicks"),
            Input({"type":"uncheck_comparisons","comparison":MATCH},"n_clicks"),
            State({"type":"signature_checkbox","disease":MATCH,"comparison":MATCH},"id"),
            State({"type":"signature_checkbox","disease":MATCH,"comparison":MATCH},"active"),
            State("filters_dropdown","value"),
            prevent_initial_call=True
        )
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="select_new_signatures"
        ),
            Output("title","id",allow_duplicate=True),
            Input("selected_genes_store","data"),
            prevent_initial_call=True
        )
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="select_by_group"
        ),
            Output({"type":"check_diseases","disease":ALL},"id"),
            Output({"type":"uncheck_diseases","disease":ALL},"id"),
            Input({"type":"check_diseases","disease":ALL},"n_clicks"),
            Input({"type":"uncheck_diseases","disease":ALL},"n_clicks"),
            prevent_initial_call=True
        )
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="update_signature_active_state"
        ),
            Output({"type":"signature_checkbox","disease":ALL,"comparison":ALL},"active"),
            Input("overview_graph","selectedNodeData"),
            State({"type":"signature_checkbox","disease":ALL,"comparison":ALL},"id"),
            State("filters_dropdown","value"),
            prevent_initial_call=True
        )
        clientside_callback(ClientsideFunction(
            namespace="clientside",
            function_name="update_signature_accordion_title"
        ),
            Output({"type":"disease_accordion","disease":ALL},"title"),
            Input("overview_graph","selectedNodeData"),
            State({"type":"disease_accordion","disease":ALL},"id"),
            prevent_initial_call=True
        )


        @callback(
            Output({"type":"signature_checkbox","disease":ALL},"value"),
            
            Input("filters_dropdown","value"),
            Input("overview_graph","selectedNodeData"),
            Input("selected_genes_store","data"),
            Input({"type":"check_diseases","disease":ALL},"n_clicks"),
            Input({"type":"uncheck_diseases","disease":ALL},"n_clicks"),
            State({"type":"signature_checkbox","disease":ALL},"id"),
            State({"type":"signature_checkbox","disease":ALL},"value"),
            State({"type":"signature_checkbox","disease":ALL},"options"),                        
            cancel=[
            Input("filters_dropdown","value"),
            Input("overview_graph","selectedNodeData"),
            Input("selected_genes_store","data"),
            Input({"type":"check_diseases","disease":ALL},"n_clicks"),
            Input({"type":"uncheck_diseases","disease":ALL},"n_clicks"),
                                ],    
            )
        def update_checkbox_signatures(selected_filter,overview_selection,menu_genes,check_diseases,uncheck_diseases,ids,cur,options):

            index = {elem["disease"]:i for i,elem in enumerate(ids)}
            match ctx.triggered_id:
                case "filters_dropdown":
                    return [[] for i in ids]
                case "overview_graph":
                    out = [[] for i in ids]
                    for node in overview_selection:
                        f = node["id"].split("_")
                        d= f[0]
                        c = f[1]
                        out[index[d]].append(c)
                    return out
                case "selected_genes_store":
                    diseases_detail,comparisons,signatures,genes_set = get_detail_subset(None, [], [], menu_genes,selected_filter)
                    if(len(genes_set)>0):
                        gene_inter = DataManager.get_instance().get_genes_intersections(diseases_detail,comparisons,signatures,genes_set,selected_filter)[1]
                        # print()
                        for s in gene_inter["id"].explode().unique():
                            f = s.split("_")
                            d= f[0]
                            c = f[1]
                            if c not in cur[index[d]]:
                                cur[index[d]].append(c)
                        return cur
                    else:
                        return cur
                case _:
                    if(isinstance(ctx.triggered_id,dict)):
                        # print("ctx.triggered_id is dict",ctx.triggered_id)
                        if("type" in ctx.triggered_id):
                            if(ctx.triggered_id["type"]=="check_diseases"):
                                d = ctx.triggered_id["disease"]
                                for i,o in enumerate(options[index[d]]):
                                    if not o["disabled"] and o["value"] not in cur[index[d]]:
                                        cur[index[d]].append(o["value"])
                            elif(ctx.triggered_id["type"]=="uncheck_diseases"):
                                d = ctx.triggered_id["disease"]
                                cur[index[d]].clear()
            return cur
        clientside_callback( ClientsideFunction(
        namespace='clientside',
        function_name='enable_check_all_btns'
        ),
            Output({"type":"check_diseases","disease":MATCH},"disabled"),
            Output({"type":"uncheck_diseases","disease":MATCH},"disabled"),
            
            Input("overview_graph","selectedNodeData"),
            State({"type":"check_diseases","disease":MATCH},"id"),
            State({"type":"uncheck_diseases","disease":MATCH},"id"),
                    prevent_initial_call=True)
        clientside_callback( ClientsideFunction(
        namespace='clientside',
        function_name='enable_check_all_btns'
        ),
            Output({"type":"check_comparisons","comparison":MATCH},"disabled"),
            Output({"type":"uncheck_comparisons","comparison":MATCH},"disabled"),
            Input("overview_graph","selectedNodeData"),
            State({"type":"check_comparisons","comparison":MATCH},"id"),
            State({"type":"uncheck_comparisons","comparison":MATCH},"id"),
                    prevent_initial_call=True)
                    

        @callback(
            Output('multi_html_legend', 'children'),
            Input('multi_legend', 'data'),
            Input('legend_tooltip_store', 'data'),
            Input('legend_hover_store', 'data'),
            Input('legend_active_store', 'data'),
            prevent_initial_call=False
        )
        def render_multi_html_legend(legend_data, tooltip_data, hover_data, active_data):
            items = []
            hovered_id = hover_data.get('hovered_id') if hover_data else None
            active_ids = active_data.get('active_ids', []) if active_data else []
            if not legend_data:
                # Return empty list instead of showing "Legend not available"
                pass
            else:
                # Genes
                if 'genes' in legend_data:
                    for k, v in legend_data['genes'].items():
                        item_id = {"type": "legend-item", "category": "genes", "name": k}
                        is_hovered = str(item_id) == hovered_id
                        is_active = str(item_id) in active_ids
                        item_style = {
                            'cursor': 'pointer',
                            'marginBottom': 2,
                            'display': 'flex',
                            'alignItems': 'center',
                            'pointerEvents': 'auto',
                            'background': '#f0f0f0' if is_hovered else ('#d1eaff' if is_active else 'none'),
                            'border': '1.5px solid #007bff' if is_active else ('1.5px solid #888' if is_hovered else '1px solid transparent'),
                            'borderRadius': '4px',
                            'transition': 'background 0.15s, border 0.15s',
                        }
                        class_list = ['legend-item']
                        if is_active:
                            class_list.append('active')
                        items.append(
                            html.Div([
                                html.Span(style={
                                    'backgroundColor': v if isinstance(v, str) else v.get('color', 'black'),
                                    'display': 'inline-block',
                                    'width': 12,
                                    'height': 12,
                                    'marginRight': 3,
                                    'borderRadius': '50%' if 'gene' in k else '0%',
                                    'border': '1px solid #333',
                                    'pointerEvents': 'auto',
                                }),
                                html.Span(k, style={'fontSize': '11px'})
                            ],
                            id=item_id,
                            n_clicks=0,
                            className=' '.join(class_list),
                            style=item_style)
                        )
                # Pathways
                if 'pathway' in legend_data:
                    for k, v in legend_data['pathway'].items():
                        shape = v.get('shape', 'rect')
                        color = v.get('color', 'gray')
                        item_id = {"type": "legend-item", "category": "pathway", "name": k}
                        is_hovered = str(item_id) == hovered_id
                        is_active = str(item_id) in active_ids
                        item_style = {
                            'cursor': 'pointer',
                            'marginBottom': 2,
                            'display': 'flex',
                            'alignItems': 'center',
                            'pointerEvents': 'auto',
                            'background': '#f0f0f0' if is_hovered else ('#d1eaff' if is_active else 'none'),
                            'border': '1.5px solid #007bff' if is_active else ('1.5px solid #888' if is_hovered else '1px solid transparent'),
                            'borderRadius': '4px',
                            'transition': 'background 0.15s, border 0.15s',
                        }
                        shape_elem = html.Span(style={
                            'display': 'inline-block',
                            'width': 0,
                            'height': 0,
                            'borderLeft': '6px solid transparent' if shape == 'triangle' else '',
                            'borderRight': '6px solid transparent' if shape == 'triangle' else '',
                            'borderBottom': f'12px solid {color}' if shape == 'triangle' else '',
                            'backgroundColor': color if shape != 'triangle' else '',
                            'marginRight': 3,
                            'pointerEvents': 'auto',
                        })
                        class_list = ['legend-item']
                        if is_active:
                            class_list.append('active')
                        items.append(
                            html.Div([
                                shape_elem,
                                html.Span(k, style={'fontSize': '11px'})
                            ],
                            id=item_id,
                            n_clicks=0,
                            className=' '.join(class_list),
                            style=item_style)
                        )
                # Diseases
                if 'diseases' in legend_data:
                    for d, color in legend_data['diseases'].items():
                        item_id = {"type": "legend-item", "category": "disease", "name": d}
                        is_hovered = str(item_id) == hovered_id
                        is_active = str(item_id) in active_ids
                        item_style = {
                            'cursor': 'pointer',
                            'marginBottom': 2,
                            'display': 'flex',
                            'alignItems': 'center',
                            'pointerEvents': 'auto',
                            'background': '#f0f0f0' if is_hovered else ('#d1eaff' if is_active else 'none'),
                            'border': '1.5px solid #007bff' if is_active else ('1.5px solid #888' if is_hovered else '1px solid transparent'),
                            'borderRadius': '4px',
                            'transition': 'background 0.15s, border 0.15s',
                        }
                        class_list = ['legend-item']
                        if is_active:
                            class_list.append('active')
                        items.append(
                            html.Div([
                                html.Span(style={
                                    'backgroundColor': color,
                                    'display': 'inline-block',
                                    'width': 12,
                                    'height': 12,
                                    'marginRight': 3,
                                    'borderRadius': '50%',
                                    'border': '1px solid #333',
                                    'pointerEvents': 'auto',
                                }),
                                html.Span(d, style={'fontSize': '11px'})
                            ],
                            id=item_id,
                            n_clicks=0,
                            className=' '.join(class_list),
                            style=item_style)
                        )
                # Comparisons (Stages)
                if 'comparisons' in legend_data:
                    for comp, style in legend_data['comparisons'].items():
                        item_id = {"type": "legend-item", "category": "comparison", "name": comp}
                        is_hovered = str(item_id) == hovered_id
                        is_active = str(item_id) in active_ids
                        item_style = {
                            'cursor': 'pointer',
                            'marginBottom': 2,
                            'display': 'flex',
                            'alignItems': 'center',
                            'pointerEvents': 'auto',
                            'background': '#f0f0f0' if is_hovered else ('#d1eaff' if is_active else 'none'),
                            'border': '1.5px solid #007bff' if is_active else ('1.5px solid #888' if is_hovered else '1px solid transparent'),
                            'borderRadius': '4px',
                            'transition': 'background 0.15s, border 0.15s',
                        }
                        swatch_style = {
                            'display': 'inline-block',
                            'width': 12,
                            'height': 12,
                            'marginRight': 3,
                            'border': '1px solid #333',
                            'pointerEvents': 'auto',
                            'backgroundColor': '#eee',  # fallback color
                        }
                        # Set different stripe directions for each stage in the legend
                        #print(f"[LEGEND DEBUG] comp: {comp}")
                        stage_gradient = 'repeating-linear-gradient(45deg, #888 0 1px, transparent 1px 4px)'
                        stage = None
                        if 'Stage' in comp:
                            parts = comp.split('Stage')
                            if len(parts) > 1:
                                stage = parts[1].strip().split()[0].replace(':','').replace('_','')
                        #print(f"[LEGEND DEBUG] extracted stage: {stage}")
                        if stage == 'I' or stage == '1':
                            stage_gradient = 'repeating-linear-gradient(90deg, #888 0 1px, transparent 1px 4px)'  # vertical
                        elif stage == 'II' or stage == '2':
                            stage_gradient = 'repeating-linear-gradient(135deg, #888 0 1px, transparent 1px 4px)'  # diagonal top right to bottom left
                        elif stage == 'III' or stage == '3':
                            stage_gradient = 'repeating-linear-gradient(0deg, #888 0 1px, transparent 1px 4px)'  # horizontal
                        elif stage == 'IV' or stage == '4':
                            stage_gradient = 'repeating-linear-gradient(45deg, #888 0 1px, transparent 1px 4px)'  # diagonal bottom left to top right
                        #print(f"[LEGEND DEBUG] stage_gradient: {stage_gradient}")
                        swatch_style['backgroundImage'] = stage_gradient
                        swatch_style['backgroundRepeat'] = 'repeat'
                        swatch_style['backgroundSize'] = '6px 6px'
                        swatch = html.Span(style=swatch_style)
                        class_list = ['legend-item']
                        if is_active:
                            class_list.append('active')
                        items.append(
                            html.Div([
                                swatch,
                                html.Span(comp, style={'fontSize': '11px'})
                            ],
                            id=item_id,
                            n_clicks=0,
                            className=' '.join(class_list),
                            style=item_style)
                        )
            # Add custom tooltip div if needed
            tooltip_div = None
            if tooltip_data and tooltip_data.get('open') and tooltip_data.get('target_id'):
                tooltip_div = html.Div(
                    tooltip_data.get('text', ''),
                    style={
                        'position': 'absolute',
                        'left': tooltip_data.get('left', 0),
                        'top': tooltip_data.get('top', 0),
                        'background': 'rgba(0,0,0,0.85)',
                        'color': 'white',
                        'padding': '2px 8px',
                        'borderRadius': '4px',
                        'fontSize': '12px',
                        'zIndex': 9999,
                        'pointerEvents': 'none',
                        'whiteSpace': 'nowrap',
                    }
                )
            return items + ([tooltip_div] if tooltip_div else [])

        # Hover state store
        @callback(
            Output('legend_hover_store', 'data'),
            Input('multi_legend', 'data'),
            [Input({"type": "legend-item", "category": "genes", "name": ALL}, 'n_mouseover')] +
            [Input({"type": "legend-item", "category": "pathway", "name": ALL}, 'n_mouseover')] +
            [Input({"type": "legend-item", "category": "disease", "name": ALL}, 'n_mouseover')] +
            [Input({"type": "legend-item", "category": "comparison", "name": ALL}, 'n_mouseover')],
            [State({"type": "legend-item", "category": "genes", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "pathway", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "disease", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "comparison", "name": ALL}, 'id')],
            prevent_initial_call=True
        )
        def update_legend_hover(legend_data, *args):
            # Extract the actual legend items from legend_data to create dynamic inputs
            all_inputs = args[:len(args)//2]
            all_states = args[len(args)//2:]
            
            for i, input_list in enumerate(all_inputs):
                if input_list:  # Check if the list exists
                    for j, value in enumerate(input_list):
                        if value and value > 0:
                            # Find the corresponding state ID
                            if i < len(all_states) and j < len(all_states[i]):
                                hovered_id = str(all_states[i][j])
                                return {'hovered_id': hovered_id}
            return {'hovered_id': None}

        # Active (clicked) state store
        @callback(
            Output('legend_active_store', 'data', allow_duplicate=True),
            Input('multi_legend', 'data'),
            Input('detail_graph', 'elements'),
            [Input({"type": "legend-item", "category": "genes", "name": ALL}, 'n_clicks')] +
            [Input({"type": "legend-item", "category": "pathway", "name": ALL}, 'n_clicks')] +
            [Input({"type": "legend-item", "category": "disease", "name": ALL}, 'n_clicks')] +
            [Input({"type": "legend-item", "category": "comparison", "name": ALL}, 'n_clicks')],
            [State({"type": "legend-item", "category": "genes", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "pathway", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "disease", "name": ALL}, 'id')] +
            [State({"type": "legend-item", "category": "comparison", "name": ALL}, 'id')],
            State('legend_active_store', 'data'),
            State('filters_dropdown', 'value'),
            prevent_initial_call=True
        )
        def handle_legend_click_print_genes(legend_data, graph_elements, *args):
            # The structure is: [genes_clicks, pathway_clicks, disease_clicks, comparison_clicks, genes_ids, pathway_ids, disease_ids, comparison_ids, active_data, selected_filter]
            if len(args) < 8:
                return dash.no_update
            
            # Extract inputs and states
            genes_clicks = args[0] if args[0] else []
            pathway_clicks = args[1] if args[1] else []
            disease_clicks = args[2] if args[2] else []
            comparison_clicks = args[3] if args[3] else []
            
            genes_ids = args[4] if args[4] else []
            pathway_ids = args[5] if args[5] else []
            disease_ids = args[6] if args[6] else []
            comparison_ids = args[7] if args[7] else []
            
            active_data = args[8] if len(args) > 8 else {}
            selected_filter = args[9] if len(args) > 9 else "Merge"
            
            # Get current active items from the store
            current_active = active_data.get('active_ids', []) if active_data else []
            
            # Handle the toggle logic for selection/deselection
            all_clicks = [genes_clicks, pathway_clicks, disease_clicks, comparison_clicks]
            all_ids = [genes_ids, pathway_ids, disease_ids, comparison_ids]
            
            for i, click_list in enumerate(all_clicks):
                if click_list:  # Check if the list exists
                    for j, value in enumerate(click_list):
                        if value and value > 0:
                            # Find the corresponding state ID
                            if i < len(all_ids) and j < len(all_ids[i]):
                                clicked_id = str(all_ids[i][j])
                                
                                # Toggle the clicked item
                                if clicked_id in current_active:
                                    # Remove from active list
                                    current_active.remove(clicked_id)
                                else:
                                    # Add to active list
                                    current_active.append(clicked_id)
                                
                                return {'active_ids': current_active}
            
            # Don't modify the active state if no changes
            return dash.no_update

        @callback(
            Output('selected_genes_legend', 'children'),
            Output('selected_genes_legend', 'style'),
            Input('legend_active_store', 'data'),
            Input('multi_legend', 'data'),
            prevent_initial_call=False
        )
        def render_selected_genes_legend(active_data, legend_data):
            active_ids = active_data.get('active_ids', []) if active_data else []

            visible_style = {
                "position": "absolute",
                "bottom": "0px",
                "left": "0px",
                "padding": "4px 6px",
                "zIndex": 200,
                "minWidth": "0px",
                "maxWidth": "200px",
                "pointerEvents": "auto",
                "backgroundColor": "rgba(255, 255, 255, 0.95)",
                "border": "1px solid #ccc",
                "borderRadius": "4px",
                "boxShadow": "0 2px 4px rgba(0,0,0,0.1)"
            }
            
            # If no items are selected, return empty children and hide the container
            if not active_ids:
                return [], {'display': 'none'}
            
            # Parse the active IDs to get the selected items
            selected_items = []
            for active_id in active_ids:
                try:
                    # Parse the string representation of the dict
                    import ast
                    item_dict = ast.literal_eval(active_id)
                    selected_items.append(item_dict)
                except:
                    continue
            
            # Create the legend content
            items = []
            
            # Add title
            items.append(
                html.Div(
                    "Selected genes",
                    style={
                        'fontWeight': 'bold',
                        'fontSize': '12px',
                        'marginBottom': '4px',
                        'color': '#333'
                    }
                )
            )
            
            # Add action buttons
            copy_button = html.Div(
                "📋 Copy to clipboard",
                id={"type": "legend-action", "action": "copy"},
                n_clicks=0,
                style={
                    'cursor': 'pointer',
                    'marginBottom': 2,
                    'display': 'flex',
                    'alignItems': 'center',
                    'pointerEvents': 'auto',
                    'background': 'none',
                    'border': '1px solid transparent',
                    'borderRadius': '4px',
                    'transition': 'background 0.15s, border 0.15s',
                    'fontSize': '11px',
                    'padding': '2px 4px'
                },
                className='legend-action-item'
            )
            
            gprofiler_button = html.Div(
                "🔗 View in gProfiler",
                id={"type": "legend-action", "action": "gprofiler"},
                n_clicks=0,
                style={
                    'cursor': 'pointer',
                    'marginBottom': 2,
                    'display': 'flex',
                    'alignItems': 'center',
                    'pointerEvents': 'auto',
                    'background': 'none',
                    'border': '1px solid transparent',
                    'borderRadius': '4px',
                    'transition': 'background 0.15s, border 0.15s',
                    'fontSize': '11px',
                    'padding': '2px 4px'
                },
                className='legend-action-item'
            )
            
            items.extend([copy_button, gprofiler_button])
            
            return items, visible_style

        # Callback to handle copy to clipboard action
        @callback(
            Output('copy_toast_store', 'data'),
            Input({"type": "legend-action", "action": "copy"}, 'n_clicks'),
            State('clipboard_text_store', 'data'),
            prevent_initial_call=True
        )
        def handle_copy_to_clipboard(n_clicks, current_clipboard_data):
            if not n_clicks or n_clicks == 0:
                return dash.no_update
            
            # Count the number of genes in the current clipboard data
            if current_clipboard_data and isinstance(current_clipboard_data, str):
                genes = [g for g in current_clipboard_data.split(';') if g.strip()]
                n = len(genes)
                if n > 0:
                    #print(f"Copy button clicked - copying {n} genes to clipboard")
                    return {"is_open": True, "message": f"Successfully copied {n} gene{'s' if n != 1 else ''}"}
                else:
                    #print("Copy button clicked - copying 0 genes to clipboard")
                    return {"is_open": True, "message": "0 genes copied"}
            
            return {"is_open": True, "message": "0 genes copied"}

        # Callback to handle view in gProfiler action
        @callback(
            Output('gprofiler_url_store', 'data'),
            Input({"type": "legend-action", "action": "gprofiler"}, 'n_clicks'),
            Input('legend_active_store', 'data'),
            Input('detail_graph', 'elements'),
            State('filters_dropdown', 'value'),
            prevent_initial_call=True
        )
        def handle_view_in_gprofiler(n_clicks, active_data, graph_elements, selected_filter):
            if not n_clicks or n_clicks == 0:
                return dash.no_update
            # Get selected signatures from the current detail graph
            selected_signatures = []
            for elem in graph_elements:
                if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                    signatures = elem.get("data", {}).get("Signatures", [])
                    selected_signatures.extend(signatures)
            selected_signatures = list(set(selected_signatures))
            # Get the selected items and extract gene information
            active_ids = active_data.get('active_ids', []) if active_data else []
            if not active_ids:
                return ""
            # Extract genes from selected legend items
            all_genes = set()
            dm = DataManager._instance
            for active_id in active_ids:
                try:
                    import ast
                    item_dict = ast.literal_eval(active_id)
                    category = item_dict.get('category')
                    name = item_dict.get('name')
                    if category and name:
                        genes = get_genes_for_legend_item(category, name, graph_elements, selected_filter, selected_signatures=selected_signatures)
                        if genes:
                            all_genes.update(genes)
                except Exception as e:
                    #print(f"Error parsing legend item {active_id}: {e}")
                    pass
            if not all_genes:
                return ""
            # Generate gProfiler URL dynamically
            try:
                # Use Ensembl IDs directly (gProfiler expects Ensembl IDs)
                ensembl_ids = []
                for gene_id in sorted(all_genes):
                    if gene_id.startswith("ENSG"):  # Ensure it's an Ensembl ID
                        ensembl_ids.append(gene_id)
                # Create gProfiler URL with Ensembl IDs (space-separated)
                genes_param = " ".join(ensembl_ids)
                gprofiler_url = f"https://biit.cs.ut.ee/gprofiler/gost?query={genes_param}&organism=hsapiens&sources=GO:BP,GO:MF,GO:CC,KEGG,REAC,WP&user_threshold=0.05&all_results=false&ordered=false&significant=true&no_iea=false&domain_scope=annotated&measure_underrepresentation=false&evcodes=false&as_ranges=false&background=0&domain_size_type=known&term_size_filter_min=3&term_size_filter_max=500&numeric_namespace=ENTREZGENE_ACC&pictograms=false&min_set_size=3&max_set_size=500"
                #print(f"Generated gProfiler URL for {len(ensembl_ids)} Ensembl IDs")
                return gprofiler_url
            except Exception as e:
                #print(f"Error generating gProfiler URL: {e}")
                return ""

        # Separate callback for printing genes that triggers when active store changes
        @callback(
            Output('legend_active_store', 'data', allow_duplicate=True),
            Input('legend_active_store', 'data'),
            Input('detail_graph', 'elements'),
            State('filters_dropdown', 'value'),
            prevent_initial_call=True
        )
        def print_genes_on_active_change(active_data, graph_elements, selected_filter):
            #print("DEBUG: print_genes_on_active_change called")
            # Get current active items from the store
            current_active = active_data.get('active_ids', []) if active_data else []
            #print(f"DEBUG: Current active items: {current_active}")
            # Get selected signatures from the current detail graph
            selected_signatures = []
            for elem in graph_elements:
                if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                    signatures = elem.get("data", {}).get("Signatures", [])
                    selected_signatures.extend(signatures)
            selected_signatures = list(set(selected_signatures))
            # Print all genes from all active legend items
            if current_active:
                #print(f"\n=== Genes from all selected legend items ({len(current_active)} items) ===")
                all_genes = {}
                for active_id in current_active:
                    try:
                        import ast
                        item_dict = ast.literal_eval(active_id)
                        category = item_dict.get('category')
                        name = item_dict.get('name')
                        if category and name:
                            # Pass selected_signatures for legend logic
                            genes = get_genes_for_legend_item(category, name, graph_elements, selected_filter, selected_signatures=selected_signatures)
                            if genes:
                                all_genes[f"{category} - {name}"] = genes
                    except Exception as e:
                        #print(f"Error parsing legend item {active_id}: {e}")
                        pass
                # Print summary and detailed gene lists
                total_unique_genes = set()
                for legend_name, genes in all_genes.items():
                    total_unique_genes.update(genes)
                    #print(f"\n{legend_name} ({len(genes)} genes):")
                    try:
                        symbols = DataManager._instance.get_symbol(genes)
                        for gene_id in sorted(genes):
                            symbol = symbols.get(gene_id, [gene_id])[0] if gene_id in symbols else gene_id
                            #print(f"  {gene_id} ({symbol})")
                    except KeyError as e:
                        #print(f"Error getting symbols: {e}")
                        for gene_id in sorted(genes):
                            #print(f"  {gene_id}")
                            pass
                #print(f"\n=== SUMMARY ===")
                #print(f"Total unique genes across all selections: {len(total_unique_genes)}")
                #print("=" * 50)
            else:
                #print("No legend items selected")
                pass
            # Don't modify the active state
            return dash.no_update

        def get_genes_for_legend_item(category, name, graph_elements, selected_filter, selected_signatures=None):
            """Get gene list for a specific legend item without printing"""
            if not graph_elements:
                return []
            dm = DataManager._instance
            if not dm:
                return []
            genes = []
            if category == "genes":
                if name == "genes":
                    # All genes in detail graph (all black dots)
                    for elem in graph_elements:
                        if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                            gene_id = elem.get("data", {}).get("id")
                            if gene_id and gene_id.startswith("ENSG"):  # Only actual genes, not signature nodes
                                genes.append(gene_id)
                elif name == "genes in multiple signatures":
                    # Genes in multiple signatures (all green dots)
                    try:
                        green_genes = get_last_green_genes()
                        genes = green_genes
                    except Exception as e:
                        genes = []
            elif category == "pathway":
                if name == "pathways linking multiple signatures":
                    # Genes connected with pathways (all dots linked to red triangles)
                    pathway_genes = set()
                    for elem in graph_elements:
                        if "source" in elem.get("data", {}) and elem.get("data", {}).get("is_pathway_edge", False):
                            source = elem.get("data", {}).get("source")
                            target = elem.get("data", {}).get("target")
                            # Check if this is a gene-to-pathway connection
                            if source and target:
                                # Find which one is the gene (not pathway)
                                for gene_elem in graph_elements:
                                    if "source" not in gene_elem.get("data", {}) and not gene_elem.get("data", {}).get("is_pathways", False):
                                        if gene_elem.get("data", {}).get("id") == source:
                                            pathway_genes.add(source)
                                        elif gene_elem.get("data", {}).get("id") == target:
                                            pathway_genes.add(target)
                    genes = list(pathway_genes)
            elif category == "disease":
                # Cancer like LIHC: print every gene in the selected cancer
                cancer_name = name
                try:
                    # Extract signature IDs from graph elements to filter the query
                    signature_ids = set()
                    for elem in graph_elements:
                        if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                            signatures = elem.get("data", {}).get("Signatures", [])
                            signature_ids.update(signatures)
                    # Get genes for only the selected signatures that match the cancer
                    intersections, items = dm.get_genes_intersections(
                        id_filter=list(signature_ids),
                        selected_filter=selected_filter
                    )
                    if items is not None and not items.empty:
                        # Filter to only include genes from the specified cancer
                        cancer_genes = []
                        for _, row in items.iterrows():
                            gene_id = row['gene']
                            signatures = row['id']
                            # Check if any of the gene's signatures contain the cancer name
                            if any(cancer_name in sig for sig in signatures):
                                cancer_genes.append(gene_id)
                        genes = cancer_genes
                except Exception as e:
                    genes = []
            elif category == "comparison":
                # Normal vs Stage N: print every gene across cancer types in stage N
                stage_name = name
                try:
                    # Extract signature IDs from graph elements to filter the query
                    signature_ids = set()
                    for elem in graph_elements:
                        if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                            signatures = elem.get("data", {}).get("Signatures", [])
                            signature_ids.update(signatures)
                    # Get genes for only the selected signatures that match the stage
                    intersections, items = dm.get_genes_intersections(
                        id_filter=list(signature_ids),
                        selected_filter=selected_filter
                    )
                    if items is not None and not items.empty:
                        # Filter to only include genes from the specified stage
                        stage_genes = []
                        # Extract stage identifier from the legend name
                        # "Normal vs StageII" -> "StageII"
                        stage_identifier = stage_name
                        if "vs" in stage_name:
                            stage_identifier = stage_name.split("vs")[1].strip()
                        # Filter signatures to only include those that match the stage
                        matching_signatures = []
                        for sig in signature_ids:
                            if stage_identifier in sig:
                                stage_pos = sig.find(stage_identifier)
                                if stage_pos != -1:
                                    before_stage = sig[:stage_pos]
                                    after_stage = sig[stage_pos + len(stage_identifier):]
                                    valid_before = not before_stage or before_stage.endswith("vs")
                                    valid_after = not after_stage or after_stage.startswith("_")
                                    if valid_before and valid_after:
                                        matching_signatures.append(sig)
                        for _, row in items.iterrows():
                            gene_id = row['gene']
                            gene_signatures = row['id']
                            if any(sig in matching_signatures for sig in gene_signatures):
                                stage_genes.append(gene_id)
                        genes = stage_genes
                except Exception as e:
                    genes = []
            return genes

        # Callback to update clipboard store when legend selection changes
        @callback(
            Output('clipboard_text_store', 'data', allow_duplicate=True),
            Input('legend_active_store', 'data'),
            Input('detail_graph', 'elements'),
            State('filters_dropdown', 'value'),
            prevent_initial_call=True
        )
        def update_clipboard_store_on_legend_change(active_data, graph_elements, selected_filter):
            #print("DEBUG: update_clipboard_store_on_legend_change called")
            # Get selected signatures from the current detail graph
            selected_signatures = []
            for elem in graph_elements:
                if "source" not in elem.get("data", {}) and not elem.get("data", {}).get("is_pathways", False):
                    signatures = elem.get("data", {}).get("Signatures", [])
                    selected_signatures.extend(signatures)
            selected_signatures = list(set(selected_signatures))
            # Collect all genes from selected legend items
            active_ids = active_data.get('active_ids', []) if active_data else []
            #print(f"DEBUG: Active IDs: {active_ids}")
            if not active_ids:
                #print("DEBUG: No active IDs, returning blank space")
                return " "
            all_genes = set()
            for active_id in active_ids:
                try:
                    import ast
                    item_dict = ast.literal_eval(active_id)
                    category = item_dict.get('category')
                    name = item_dict.get('name')
                    if category and name:
                        genes = get_genes_for_legend_item(category, name, graph_elements, selected_filter, selected_signatures=selected_signatures)
                        if genes:
                            all_genes.update(genes)
                except Exception as e:
                    #print(f"Error parsing legend item {active_id}: {e}")
                    pass
            if not all_genes:
                #print("DEBUG: No genes found, returning blank space")
                return " "
            # Return as semicolon-separated string
            gene_str = ";".join(sorted(all_genes))
            #print(f"DEBUG: Returning gene string: '{gene_str}'")
            return gene_str

        @callback(
            Output('detail_graph', 'elements', allow_duplicate=True),
            Output('detail_graph', 'stylesheet', allow_duplicate=True),
            Input('legend_active_store', 'data'),
            State('detail_graph', 'elements'),
            State('detail_graph', 'stylesheet'),
            State('detail_graph', 'elements'),  # for get_genes_for_legend_item
            State('filters_dropdown', 'value'),
            prevent_initial_call=True
        )
        def highlight_nodes_on_legend(legend_active_data, elements, stylesheet, graph_elements, selected_filter):
            import copy
            # Defensive copy
            elements = copy.deepcopy(elements)
            stylesheet = copy.deepcopy(stylesheet)
            highlight_gene_ids = set()
            if legend_active_data and 'active_ids' in legend_active_data:
                for active_id in legend_active_data['active_ids']:
                    import ast
                    item_dict = ast.literal_eval(active_id)
                    category = item_dict.get('category')
                    name = item_dict.get('name')
                    if category and name:
                        genes = get_genes_for_legend_item(category, name, graph_elements, selected_filter)
                        highlight_gene_ids.update(genes)
            # Update node highlight styles
            if elements is not None:
                update_node_highlight_styles(elements, list(highlight_gene_ids))
            # Ensure highlight style is in stylesheet
            if highlight_gene_ids and not any(s.get("selector") == ".legend-highlight" for s in stylesheet):
                stylesheet.append(legend_highlight_stylesheet)
            return elements, stylesheet

        @callback(
            Output('selected_genes_store', 'data', allow_duplicate=True),
            Input('remove_all_boxplots_button', 'n_clicks'),
            State('selected_genes_store', 'data'),
            prevent_initial_call=True
        )
        def remove_all_boxplots(n_clicks, current):
            if n_clicks:
                # Clear all selected genes and pathway genes, as if all boxplots were closed
                return {"selected": [], "from_pathways": {"ids": [], "genes": [], "signatures": []}, "covered_signatures": [], "pathway_signature_map": {}}
            raise dash.exceptions.PreventUpdate()

        @callback(
            Output('remove_all_boxplots_button_container', 'children'),
            Input('selected_genes_store', 'data'),
            prevent_initial_call=False
        )
        def show_remove_all_boxplots_button(selected_genes_store):
            # Check if there are any boxplots to show the button for
            has_boxplots = False
            if selected_genes_store:
                if selected_genes_store.get('selected'):
                    if len(selected_genes_store['selected']) > 0:
                        has_boxplots = True
                # Check pathway genes
                if selected_genes_store.get('from_pathways') and selected_genes_store['from_pathways'].get('genes'):
                    for gene_group in selected_genes_store['from_pathways']['genes']:
                        if gene_group and len(gene_group) > 0:
                            has_boxplots = True
            if has_boxplots:
                # Use dcc.Markdown to inject a <style> block for the hover effect
                return [
                    dcc.Markdown(
                        '    <style>\n.remove-all-boxplots-btn:hover { background: #e0e0e0 !important; border-color: #bbb !important; }\n</style>',
                        dangerously_allow_html=True,
                        style={'display': 'none'}
                    ),
                    html.Button(
                        ['🗑️ ', 'Remove all boxplots'],
                        id='remove_all_boxplots_button',
                        className='remove-all-boxplots-btn',
                        style={
                            'padding': '3px 10px',
                            'fontSize': '0.85em',
                            'background': '#f5f5f5',
                            'color': '#333',
                            'border': '1px solid #ccc',
                            'borderRadius': '4px',
                            'cursor': 'pointer',
                            'boxShadow': 'none',
                            'transition': 'background 0.2s, color 0.2s, border 0.2s',
                            'outline': 'none',
                            'display': 'flex',
                            'alignItems': 'center',
                            'gap': '4px',
                        }
                    )
                ]
            return None

        # @callback(
        #     Output('pathway_signature_warning_store', 'data'),
        #     Input('overview_graph', 'selectedNodeData'),
        #     State('selected_genes_store', 'data'),
        #     State('filters_dropdown', 'value'),
        #     prevent_initial_call=True
        # )
        # def validate_pathway_signature_consistency(selected_signatures, current_store, selected_filter):
        #     """
        #     When signatures are manually deselected, check if pathways still have sufficient signature coverage.
        #     Instead of removing pathways automatically, provide warnings to the user.
        #     """
        #     if not current_store or not current_store.get('pathway_signature_map'):
        #         return dash.no_update
            
        #     # Extract signature IDs from the node data
        #     current_signatures = set()
        #     if selected_signatures:
        #         for node in selected_signatures:
        #             if isinstance(node, dict) and 'id' in node:
        #                 current_signatures.add(node['id'])
            
        #     pathway_map = current_store['pathway_signature_map']
        #     pathways_with_insufficient_coverage = []
            
        #     # Check each pathway's signature requirements
        #     for pathway_id, required_signatures in pathway_map.items():
        #         required_set = set(required_signatures)
        #         available_signatures = current_signatures.intersection(required_set)
                
        #         # If less than threshold of required signatures are available, flag the pathway
        #         coverage_ratio = len(available_signatures) / len(required_set) if len(required_set) > 0 else 0
        #         if coverage_ratio < PATHWAY_SIGNATURE_COVERAGE_THRESHOLD:
        #             pathways_with_insufficient_coverage.append({
        #             'pathway_id': pathway_id,
        #             'coverage_ratio': coverage_ratio,
        #             'available_signatures': len(available_signatures),
        #             'required_signatures': len(required_set)
        #             })
            
        #     # Return warning data if there are pathways with insufficient coverage
        #     if pathways_with_insufficient_coverage:
        #         warning_data = {
        #             'show_warning': True,
        #             'pathways': pathways_with_insufficient_coverage,
        #             'message': f"Some pathways may have insufficient signature coverage after recent signature deselection."
        #         }
        #         return warning_data
            
        #     return {'show_warning': False}
        

        
        # @callback(
        #     Output('pathway_signature_warning_alert', 'children'),
        #     Output('pathway_signature_warning_alert', 'is_open'),
        #     Input('pathway_signature_warning_store', 'data'),
        #     prevent_initial_call=True
        # )
        # def display_pathway_signature_warning(warning_data):
        #     """
        #     Display a warning alert when pathways have insufficient signature coverage.
        #     """
        #     if not warning_data or not warning_data.get('show_warning'):
        #         return [], False
            
        #     pathways = warning_data.get('pathways', [])
        #     if not pathways:
        #         return [], False
            
        #     # Create warning message
        #     pathway_names = []
        #     for pathway_info in pathways:
        #         pathway_id = pathway_info['pathway_id']
        #         coverage_pct = pathway_info['coverage_ratio'] * 100
        #         pathway_names.append(f"{pathway_id} ({coverage_pct:.1f}% coverage)")
            
        #     warning_message = f"Warning: The following pathways may have insufficient signature coverage: {', '.join(pathway_names)}. Consider reselecting the required signatures or removing these pathways manually."
            
        #     alert = dbc.Alert([
        #         html.H6("Pathway Signature Coverage Warning", className="alert-heading"),
        #         html.P(warning_message),
        #         html.Hr(),
        #         html.P("You can manually remove pathways from the Pathways section if needed.", className="mb-0")
        #     ], color="warning", dismissable=True)
            
        #     return alert, True
