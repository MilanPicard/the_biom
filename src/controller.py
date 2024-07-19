import re
import time
# from dash_extensions.enrich import Dash, html, Input, Output, State, callback,ctx,ALL
from dash import Dash, html, Input, Output, State, callback,ctx,ALL
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
import overview as ov
import detail_graph
import utils
import stats
from  data_manager import DataManager
# import cProfile,pstats,io
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
                Input("genes_menu_select","value"),
                Input({"type":"selected_gene_button","gene":ALL},"n_clicks"),
                Input("detail_graph","selectedNodeData"),
                Input("detail_graph","tapNodeData"),
                Input("filters_dropdown","value"),
                Input("pathways_menu_select","value"),
                Input({"type":"selected_pathway_button","pathway":ALL},"n_clicks"),
                Input("mono_graph","selectedNodeData"),
                Input("mono_graph","tapNodeData"),
                State('selected_genes_store','data'),prevent_initial_call=True
        )
        def add_remove_gene(menu_select,button,fromGraph,multip_tap,selected_filter,menu_pathway,pathway_button,fromMonoGraph,monotap,current):
            match ctx.triggered_id:
                case "genes_menu_select":
                    if menu_select is not None and menu_select != "None" and menu_select not in current["selected"]:
                        current["selected"].append(menu_select)
                    else:
                        raise dash.exceptions.PreventUpdate()
                case "detail_graph":
                    
                    print("fromGraph",fromGraph,multip_tap,ctx.triggered_prop_ids)
                    if("detail_graph.tapNodeData" in ctx.triggered_prop_ids):
                        gene = multip_tap
                        if (not "is_pathways" in gene or not gene["is_pathways"]) and  gene["id"] not in current["selected"]:
                            current["selected"].append(gene["id"])
                        else:
                            raise dash.exceptions.PreventUpdate()
                    else:
                        if fromGraph is not None and len(fromGraph)>0:
                            added_any = False
                            for gene in fromGraph:
                                if (not "is_pathways" in gene or not gene["is_pathways"]) and  gene["id"] not in current["selected"]:
                                    current["selected"].append(gene["id"])
                                    added_any = True
                            if not added_any:
                                raise dash.exceptions.PreventUpdate()
                        else:
                            raise dash.exceptions.PreventUpdate()
                case "mono_graph":
                    if fromMonoGraph is not None and len(fromMonoGraph)>0:
                        for gene in fromMonoGraph:
                            if (not "is_pathways" in gene or not gene["is_pathways"]) and  gene["id"] not in current["selected"]:
                                current["selected"].append(gene["id"])
                    else:
                        raise dash.exceptions.PreventUpdate()
                case "pathways_menu_select":
                    if menu_pathway is not None and menu_pathway != "None" and menu_pathway not in current["from_pathways"]:
                        current["from_pathways"]["ids"].append(menu_pathway)
                        current["from_pathways"]["genes"].append(DataManager.get_instance().get_genes(selected_filter,pathway = menu_pathway))
                    else:
                        raise dash.exceptions.PreventUpdate()

                case "filters_dropdown":
                    current["selected"]=[]
                    current["from_pathways"]["ids"]=[]
                    current["from_pathways"]["genes"]=[]
                case _:
                    if "gene" in ctx.triggered_id:
                        current["selected"].remove(ctx.triggered_id["gene"])
                    else:
                        index = current["from_pathways"]["ids"].index(ctx.triggered_id["pathway"])
                        current["from_pathways"]["ids"].remove(ctx.triggered_id["pathway"])
                        del current["from_pathways"]["genes"][index]

            return current        
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
            print(store_data)
            return update_button(store_data["from_pathways"]["ids"], cur_children,attr="pathway",get_symbol=DataManager.get_instance().get_pathway_label) 
        def update_button(store_data, cur_children,attr,get_symbol):
            patched = Patch()
            if cur_children is None:
                cur_children = []
            n=len(cur_children)
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
                            btn_style = {}
                            btn = html.Button(get_symbol(store_data[i]),id={"type":f"selected_{attr}_button",attr:store_data[i]},className="btn",style=btn_style)
                            patched.append(btn)
                            n+=1
                    for i,g in enumerate(store_data):
                        tc = utils.get_text_color(color_scheme[i%len(color_scheme)])
                        patched[i]["props"]["style"]={
                            "--bs-btn-bg":color_scheme[i%len(color_scheme)],
                            "--bs-btn-hover-bg":color_scheme[i%len(color_scheme)],
                            "--bs-btn-hover-color":tc,
                            "--bs-btn-color":tc}
                else:
                    raise dash.exceptions.PreventUpdate()

            return patched         


        @callback(
                Output('overview_graph','elements'),
                Output('overview_graph_layout','data'),
                Input("disease_filter","value"),
                Input("comparisons_filter","value"),
                Input("overview_graph","selectedNodeData"),
                Input("detail_graph","selectedNodeData"),
                Input("selected_genes_store","data"),
                State("filters_dropdown","value"),
                State("overview_graph","elements")
                )
        def update_overview(diseases,
                            comparisons_filter,
                            selectedSign,
                            selected_detail_gene,
                            selected_genes_store,
                            selected_filter,
                            cur_elems):
            overview_graph_layout = dash.no_update
            if any(["Filter" in c["data"] and c["data"]["Filter"]!=selected_filter for c in cur_elems]):
                cur_elems = ov.get_elements(DataManager.get_instance(),selected_filter=selected_filter)
                overview_graph_layout = {"rerun":True}
                print(cur_elems)
            else:
                genes = set(selected_genes_store["selected"])
                for p in selected_genes_store["from_pathways"]["genes"]:
                    genes.update(p)

                classes = ["highlight","half_highlight"]
                
                signs = set([i["id"] for i in selectedSign]) if selectedSign is not None else set()
                # if selectedSign is not None:
                #     for i in selectedSign:
                #         if i["Cancer"] not in diseases:
                #             diseases.append(i["Cancer"])
                #         if i["Comparison"] not in comparisons_filter:
                #             comparisons_filter.append(i["Comparison"])
                print(signs,diseases,comparisons_filter)
                for i in cur_elems:
                    if "fake" in i["data"] and i["data"]["fake"]:
                        continue               
                    if "Cancer" in i["data"]:
                        if "classes" not in i:
                            i["classes"]=" ".join([i["data"]["Cancer"],"highlight"])
                        if i['data']["Cancer"] in diseases and i['data']["Comparison"] in comparisons_filter:
                            c = "highlight"
                            if( selected_genes_store is not None):
                                c = "half_highlight"
                                if len(genes.intersection(i["data"]["Signature"]))!=0 :
                                    c = "highlight"
                            i["classes"]=utils.switch_class(i["classes"],[c],classes)
                        else:
                            if i["data"]["id"] in signs:
                                i["classes"]=utils.switch_class(i["classes"],["highlight"],classes)
                            else:
                                i["classes"]=utils.switch_class(i["classes"],[],classes)
                    else:
                        if "source" not in i["data"]:
                            print(i)
                        if "classes" not in i:
                            i["classes"]="pathway" if i["data"]["type"]=="pathway" else ""

                        src = i["data"]["source"].split("_")
                        tgt = i["data"]["target"].split("_")
                        # src.remove("vs")
                        # tgt.remove("vs")
                        src_dis = src[0]
                        src_comp = "_".join(src[1:])
                        tgt_comp = "_".join(tgt[1:])
                        tgt_dis = tgt[0]
                        if(src_dis in diseases and tgt_dis in diseases and src_comp in comparisons_filter and tgt_comp in comparisons_filter):
                            if "highlight" not in i["classes"]:
                                i["classes"] = " ".join(i["classes"].split(" ")+["highlight"])
                        else:
                            if "highlight" in i["classes"]:
                                classes = i["classes"].split(" ")
                                classes.remove("highlight")
                                i["classes"]="_".join(classes)
            return cur_elems,overview_graph_layout

        @callback(Output('data_overview_selected','value'),
                Input('overview_graph','selectedNodeData'), prevent_initial_call=True
                )
        def update_overview_selected_data(data):
            if data is None:
                return ""
            return ";".join([i["id"] for i in data])
        @callback(Output('data_menu_selected','value'),
                Input('disease_filter','value'),
                    prevent_initial_call=True
                )
        def update_menu_selected_data(data):
            if data is None:
                return ""
            return ";".join(data)
        @callback(Output('disease_filter','value'),
                Input("selected_genes_store","data"),
                State("disease_filter","value"),
                    prevent_initial_call=True
                )
        def update_menu_disease_filter_data(genes,current):
            genes_set = set(genes["selected"])
            for p in genes["from_pathways"]["genes"]:
                genes_set.update(p)
            genes_diseases= DataManager.get_instance().get_diseases_from_genes(genes_set)
            for d in genes_diseases:
                if d not in current:
                    current.append(d)
            return current
        @callback(
                Output('comparisons_filter','value'),
                Input("selected_genes_store","data"),
                State('comparisons_filter','value'),
                    prevent_initial_call=True
                )
        def update_menu_comparison_filter_data(genes,cur_state):
            genes_set = set(genes["selected"])
            for p in genes["from_pathways"]["genes"]:
                genes_set.update(p)
            if len(genes_set)>0:
                for c in  DataManager.get_instance().get_comparisons_from_genes(genes_set):
                    if c not in cur_state:
                        cur_state.append(c)
            return cur_state

        @callback(Output('detail_graph','elements'),
                Output('detail_graph','stylesheet'),
                Output('detail_graph','layout'),
                Output('detail_graph_pos','data'),
                Output('multi_legend',"data"),
                Output('multi_export',"data"),
                #   Input('overview_graph','selectedNodeData'),
                Input('disease_filter','value'),
                Input('comparisons_filter','value'),
                Input('data_overview_selected','value'),
                Input('selected_genes_store','data'),
                Input("fake_graph_size","data"),
                Input("filters_dropdown","value"),
                State('detail_graph','elements'),
                State("detail_graph_pos","data"),            
                State('detail_graph','stylesheet'),
                    prevent_initial_call=True,
                        cancel=[
                Input('disease_filter','value'),
                Input('comparisons_filter','value'),
                Input('data_overview_selected','value'),
                Input('selected_genes_store','data'),
                Input("fake_graph_size","data"),
                Input("filters_dropdown","value"),
                                ],    background=background,
                )
        def display_detail_graph(diseases,comparisons,signatures,menu_genes,fake_graph_size,selected_filter,existing_elements,detail_pos_store,current_stylesheets):
            if ctx.triggered_id=="fake_graph_size" :
                if (fake_graph_size is None or "just_redraw" not in fake_graph_size or not fake_graph_size["just_redraw"]):
                    raise dash.exceptions.PreventUpdate()
                else:
                    return detail_graph.redraw(existing_elements,detail_pos_store,1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],current_stylesheets)
            if(all([len(diseases)==0 or len(comparisons)==0 ,signatures is None or signatures ==""])):
                return [],[],{"name":"preset"},{},dash.no_update,dash.no_update
            if len(diseases)!=0 or len(signatures)!=0:
                diseases, comparisons, signatures, genes_set = get_detail_subset(diseases, comparisons, signatures, menu_genes)
                r = detail_graph.display_detail_graph(diseases,signatures,genes_set,existing_elements,detail_pos_store if detail_pos_store is not None else dict(),1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],selected_filter,comparisons)
                return r
            else:
                return existing_elements,[],{"name":"preset"},{},dash.no_update,dash.no_update

        def get_detail_subset(diseases, comparisons, signatures, menu_genes):
            if diseases is None:
                diseases =""
            if signatures is None:
                signatures = ""
            if len(signatures)>0:
                diseases =[i for i in diseases]
                comparisons =[i for i in comparisons]
                for s in signatures.split(";"):
                    cancer,comp,fil = s.split("_")
                    if cancer not in diseases:
                        diseases.append(cancer)
                    if comp not in comparisons:
                        comparisons.append(comp)
                        
            genes_set = set()
            for p in menu_genes["from_pathways"]["genes"]:
                genes_set.update(p)
            genes_set = menu_genes["selected"] + sorted(list(genes_set.difference(menu_genes["selected"])))

            diseases = list(filter(lambda a: len(a)>0,diseases))
            signatures = list(filter(lambda a: len(a)>0,signatures.split(";")))
            return diseases,comparisons,signatures,genes_set
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
                                Input('disease_filter','value'),
                Input('comparisons_filter','value'),
                Input('selected_genes_store','data'),
                State('mono_graph','elements'),
                State("mono_graph_pos","data"),            
                State('mono_graph','stylesheet'),

                    prevent_initial_call=True
                )
        def display_mono_graph(tapNodeData,fake_graph_size,selected_filter,diseases,comparisons,menu_genes,existing_elements,detail_pos_store,current_stylesheets):
            d = None
            c = None
            s = None
            if(len(diseases)==1 and len(comparisons)==1):
                d = diseases[0]
                c = comparisons[0]
                s = f"{d}_{c}_{selected_filter}"
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
            r = detail_graph.display_detail_graph([d],[s],list(genes_set),existing_elements,detail_pos_store if detail_pos_store is not None else dict(),1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],selected_filter,[c],all_pathway=True)
            return *r,s

        # @callback(Output("data_gene_detail_selected","value"),
        #           Input("detail_graph","selectedNodeData"))
        # def data_gene_detail_selected(nodes):
        #     if nodes is None:
        #         return ""
        #     return ";".join([i["id"] for i in nodes])

        @callback(
            Output("activation_boxplot","figure"),
            Output("activation_boxplot","className"),
            Output("overview_graph","stylesheet"),
            Output("box_plots_to_style","data"),
            Output("box_plots_stats","data"),
            Input("data_menu_selected","value"),
            Input("data_overview_selected","value"),
            # State("overview_graph","selectedNodeData" ),g
            Input("selected_genes_store","data"),
                Input("box_categories","data"),
                Input('comparisons_filter','value'),
                Input('disease_filter','value'),
                Input("filters_dropdown","value"),

                State("overview_graph","elements"),
                State("overview_graph","stylesheet"),
                prevent_initial_call=False
        )
        def update_box_plot(menu_selected_diseases,overview_selected,menu_selected,selected_boxcategories,comparisons,disease_filter,selected_filter,overview_elements,overview_stylesheets):
            items = []
            if menu_selected is not None and (len(menu_selected["selected"])!=0 or len(menu_selected["from_pathways"]["ids"])>0):
                genes_set = set()
                for p in menu_selected["from_pathways"]["genes"]:
                    genes_set.update(p)
                items += menu_selected["selected"] + sorted(list(genes_set.difference(menu_selected["selected"])))
            stylesheets =overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(DataManager.get_instance())
            stylesheets = [s for s in stylesheets if not(s["selector"].startswith("edge#"))]
            if(len(items)>0):
                diseases,comparisons,signatures,genes_set =get_detail_subset(disease_filter, comparisons, overview_selected, menu_selected)
                selected_patient_and_genes = DataManager.get_instance().get_activations(items,diseases if len(selected_boxcategories["diseases"])==0 else selected_boxcategories["diseases"],selected_boxcategories["comparisons"]).sort_values(["box_category"])
                box_categories_tohighlight =dict({i:set() for i in items})
                gene_inter = DataManager.get_instance().get_genes_intersections(diseases,comparisons,signatures,genes_set,selected_filter)[1]
                
                if gene_inter is not None:
                    gene_inter = gene_inter.set_index("gene")["id"].to_dict()
                    for i in items:

                        signature_to_highlight = gene_inter[i]
                        for s in signature_to_highlight:
                            s = s.split("_")
                            d = s[0]
                            for stage in s[1].split("vs"):
                                box_categories_tohighlight[i].add(f"{d}_{stage}")
                box_categories = sorted(pd.unique(selected_patient_and_genes["box_category"]).tolist())
                symbols = list(map(lambda s : " ".join(s),DataManager.get_instance().get_symbol(items).to_list()))
                selected_patient_and_genes =selected_patient_and_genes.rename(columns = dict(zip(items,symbols)))
                for i in overview_elements:
                    if "elems" in i["data"] and  any([j in items for j in i["data"]["elems"]]):
                        stylesheets.append({"selector":'edge#'+i["data"]["id"],"style":{"line-color":"red","width":6}})
                        # if("classes" in i):
                            # i["classes"] = " ".join(set(i["classes"].split(" ") + ["highlight_edge"]))
                        # else:
                            # i["classes"] =  ["highlight_edge"]
                    # else:
                        # if("classes" in i) and "highlight_edge" in i["classes"]:
                            # classes:list = i["classes"].split(" ")
                            # classes.remove("highlight_edge")
                            # i["classes"] = " ".join(classes)
                # if(len(items)==1):
                #     box = px.box(selected_patient_and_genes,x="box_category",y=symbols[0],color_discrete_sequence=detail_graph.get_color_scheme(items),labels={"box_category":""})
                #     if(box.layout.margin.t is not None and box.layout.margin.t>20):
                #         box.layout.margin.t=20
                #     return box,"visible_plot",stylesheets,{"categories":box_categories,"genes":items}

                # else:
                dfs = []
                color_scheme=detail_graph.get_color_scheme(items)
                disease_cmp ={i:f"rgb({int(255*j[0])},{int(255*j[1])},{int(255*j[2])})" for i,j in DataManager.get_instance().get_disease_cmap().items()}
                for i in range(len(items)):
                    df = selected_patient_and_genes.filter(("box_category",symbols[i],"Cancer"))
                    df = df.rename({symbols[i]:"expression"},axis=1)
                    df["gene"] = symbols[i]
                    dfs.append(df)
                df = pd.concat(dfs)
                box = px.box(df,x="box_category",y="expression",color="Cancer",facet_row="gene",color_discrete_map=disease_cmp,labels={"box_category":"","expression":"expression (log2(TPM+1))","LUAD":"a"})
                box = plotly.subplots.make_subplots(rows=len(items),cols=1,shared_xaxes=True,shared_yaxes=True,row_titles=symbols)
                for i in range(len(items)):
                    df = selected_patient_and_genes.filter(("box_category",symbols[i],"Cancer"))
                    df = df.rename({symbols[i]:"expression"},axis=1)
                    df["gene"] = symbols[i]
                    is_highlighted = df["box_category"].isin(box_categories_tohighlight[items[i]])

                    for c in diseases:
                        trace = go.Box(y=df[is_highlighted].loc[df[is_highlighted]["Cancer"]==c]["expression"],x=df[is_highlighted].loc[df[is_highlighted]["Cancer"]==c]["box_category"],name=c,showlegend =i==0,legendgroup=c,fillcolor=re.sub(r"rgb\(([0-9]+,[0-9]+,[0-9]+)\)",r"rgba(\1,64)",disease_cmp[c]))
                        trace.line={'color':"black"}
                        box.add_trace(trace,row=i+1,col=1)
                        trace = go.Box(y=df[~is_highlighted].loc[df[~is_highlighted]["Cancer"]==c]["expression"],x=df[~is_highlighted].loc[df[~is_highlighted]["Cancer"]==c]["box_category"],name=c,showlegend =i==0,legendgroup=c)
                        trace.line={'color':disease_cmp[c]}
                        box.add_trace(trace,row=i+1,col=1)

                # box = go.Figure(traces)
                for axis in box.select_yaxes(col=1):
                    axis.update(title={"text":"expression (log2(TPM+1))"})
                if(box.layout.margin.t is not None and box.layout.margin.t>20):
                    box.layout.margin.t=20
                return box ,"visible_plot",stylesheets,{"categories":box_categories,"genes":items},{"stats":stats.ttest(dfs,[0.1,0.05,0.01])}
                
            else:
                return go.Figure(data=[
                    ]),"hidden_plot",stylesheets,{"categories":[],"genes":[]},{"stats":[]}
        # clientside_callback(
        #     ClientsideFunction(
        #         namespace='clientside',
        #         function_name='box_plots_stats'
        #     ),Input("activation_boxplot","restyleData"),
        #     Input("activation_boxplot","figure"),
        #     State("box_plots_stats","data"),
        #     State("do_box_plots_stats","data"),
        #         prevent_initial_call=False)

        # @callback(
        #     Output("activation_heatmap","figure"),
        #     Output("activation_heatmap","className"),
        #     Input("overview_graph","selectedNodeData" ),
        #     Input("detail_graph","selectedNodeData"),
        #     Input("activation_boxplot","clickData")
        # )
        # def update_heatmap(overview_selected,detail_selected,selected_box):
        #     if(detail_selected is not None and len(detail_selected)==1 and selected_box is not None and len(selected_box)==1):
        #         g = detail_selected[0]['id']
        #         disease,stage = selected_box['points'][0]['x'].split("_")

        #         selected_patient_and_genes = DataManager.get_instance().get_activations(g,[disease]).sort_values(by=["Stage",g])
        #         y_labels = []
        #         for item in selected_patient_and_genes.itertuples():
        #             if len(y_labels)==0 or y_labels[-1]!=item.box_category:
        #                 y_labels.append(item.box_category)
        #             else:
        #                 y_labels.append("")
        #         fig =  px.imshow(np.expand_dims(selected_patient_and_genes[g].to_numpy(),1),y=y_labels,color_continuous_scale="turbo")
        #         fig.update_xaxes(showticklabels=False)

        #         return fig,"visible_plot"
        #     else:
        #         return go.Figure(data=[
        #         ]),"hidden_plot"

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
        #     # Output("detail_graph_tooltip","bbox"),
            Input("detail_graph","mouseoverNodeData"),
            Input("detail_graph","mouseoverEdgeData"),
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
        #     # Output("mono_graph_tooltip","bbox"),
            Input("mono_graph","mouseoverNodeData"),
            Input("mono_graph","mouseoverEdgeData"),
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
        #     # Output("detail_graph_tooltip","bbox"),
            Input("overview_graph","mouseoverNodeData"),
            Input("overview_graph","mouseoverEdgeData"),
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
        function update_width_start(n1,n2,e_width,e_height,state,detail_graph_pos,dw,fake_graph_size){
        var graph_width = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientWidth:document.getElementById("mono_graph").clientWidth;
        var graph_height = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientHeight:document.getElementById("mono_graph").clientHeight;

            var e = dash_clientside.callback_context.triggered_id=="full_col"?e_height:e_width;
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
                //document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true,"animate":false}).run();
                layout_overview();

                return fake_graph_size;
                }
                else
                    return dash_clientside.no_update;

            }
            if(e["type"]=="mousedown"){
                let in_tooltip = false;
                for(let t of document.querySelectorAll(".dcc-tooltip-bounding-box")){
                    in_tooltip = in_tooltip || t.contains(e.target);
                }

                if(e.target.classList.contains("resize_span") && !in_tooltip){
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
                        if((dash_clientside.callback_context.triggered_id==="full_col" && "second_row_div" == e.target.id) || (possibleTargets.includes(e.target.id) && (e.target.clientHeight-e.offsetY<50)&& dash_clientside.callback_context.triggered_id==="move_in_ov")){
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
                    //document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
                    layout_overview();

                    fake_graph_size["width"]=graph_width;
                    fake_graph_size["height"]=graph_height;
                    fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                    return fake_graph_size;
                }
                if(state["height"]["is_resizing"]){
                    state["height"]["is_resizing"]=false;
                    dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
                    //document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
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
            Input("full_col","n_events"),
            State("move_in_ov","event"),
            State("full_col","event"),
            State("resize_state","data"),
            State("detail_graph_pos","data"),
            State("detail_col","width"),
            State('fake_graph_size','data'),
            prevent_initial_call=False
                    )

        clientside_callback(
                """
        function update_width(n,n_height,e_width,e_height,ow,dw,resize_state,fake_graph_size){
        //console.time("update_width");
        var e = dash_clientside.callback_context.triggered_id=="full_col"?e_height:e_width;

        function find_row(elem){
            var p =  elem;
            while(!p.classList.contains("row"))
                p = p.parentNode;
            return p;
        }
        if(fake_graph_size ===null)
            fake_graph_size = {};
        var graph_width = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientWidth:document.getElementById("mono_graph").clientWidth;
        var graph_height = document.getElementById("detail_graph").clientWidth!=0?document.getElementById("detail_graph").clientHeight:document.getElementById("mono_graph").clientHeight;

        if(resize_state["width"]["is_resizing"]){
            var overviewCol = document.getElementById("overview_col");
            var detailCol = document.getElementById("detail_col");
            if(e["type"]=="mousemove" &&e.target.tagName==="CANVAS" && (overviewCol.contains(e.target) || detailCol.contains(e.target))){
                if(overviewCol.contains(e.target) && ow>2){
                    var w = ow;
                    var width = find_row(e.target).clientWidth;
                    while(w>2&&e.offsetX<(w-1)*width/12){
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
                        while(w>2&&e.target.clientWidth-e.offsetX<(w-1)*width/12){
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
        if(e_width==e && resize_state["height"]["is_resizing"]){
            if(e["type"]=="mousemove" ){
                if(ajust_flex(e.offsetY)){
                    fake_graph_size["width"]=graph_width;
                    fake_graph_size["height"]=graph_height;
                    fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                    fake_graph_size["just_redraw"]=true;
                    dash_clientside.set_props("fake_graph_size", fake_graph_size);
                }
            }
        }

        if(e_height==e && resize_state["height"]["is_resizing"]){
            if(e["type"]=="mousemove" ){
                var upHeight = document.getElementById("overview_col").parentNode.clientHeight;
                if(ajust_flex(e.offsetY+upHeight)){
                    fake_graph_size["width"]=graph_width;
                    fake_graph_size["height"]=graph_height;
                    fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
                    fake_graph_size["just_redraw"]=true;
                    dash_clientside.set_props("fake_graph_size", fake_graph_size);
                }
            }
        }
        if(e == e_width &&  !resize_state["width"]["is_resizing"]){
        if(e["type"]=="click"){
        if(document.getElementById("detail_graph").contains(e.target))
        {
            tapMultiSignPathway();

        }

        }
        }
        //console.timeEnd("update_width");
        return dash_clientside.no_update;
        }
                """,
            Output("overview_col","width"),
            Output("detail_col","width"),
            Input("move_in_ov","n_events"),
            Input("full_col","n_events"),
            State("move_in_ov","event"),
            State("full_col","event"),
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
            match pathname:
                case "/":
                    about["display"]="inline"
                    main["display"]="none"
                case "/about":
                    main["display"]="inline"
                    about["display"]="none"
            return about,main

        # clientside_callback(
        #     """
        #     function update_box_style(box_plots_to_style,elements){
        #     console.log(box_plots_to_style);
        #     if(box_plots_to_style["genes"].length>0){
        #         var cy = document.getElementById("detail_graph")['_cyreg']["cy"];
        #         var traces = document.querySelector("#activation_boxplot .boxlayer").children;
        #         for(var i =0;i<box_plots_to_style["genes"].length;i+=1){
        #             var gene = box_plots_to_style["genes"][i];
        #             var trace = traces[i];
        #             const cy_elem = cy.elements("#" + gene).toArray();
        #                 if(cy_elem!==undefined){
        #                 console.log(gene,cy_elem);
        #                 var signatures = cy_elem[0].data("Signatures");
        #                 for(var j=0;j<box_plots_to_style["categories"].length;j+=1){
        #                     var category = box_plots_to_style["categories"][j].split("_");
        #                     var found = false;
        #                     for(var k=0;k<signatures.length&& !found;k+=1){
        #                         var s = signatures[k].split("_");
        #                         if(s[0]==category[0]){
        #                             var stages = s[1].split("vs");
        #                             if(category[1]==stages[0] || category[1]==stages[1]){
        #                                 found=true;
        #                             }
        #                         }
        #                     }
        #                     if(found){
        #                         //trace.children[j].style.strokeWidth=4;
        #                     }else{
        #                         trace.children[j].classList.add("nothighlight");
        #                     }
        #                     console.log(signatures,category);
        #                 }
        #                 console.log(gene,trace,cy_elem);
        #             }
        #         }
        #     }
        #     }
        #     """,
        #     Input("box_plots_to_style","data"),
        #     Input("detail_graph","elements"),
            
        #     prevent_initial_call=True

        #     )


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

        # clientside_callback(    ClientsideFunction(
        #         namespace='clientside',
        #         function_name='highlight_pathway_neighbourhood'
        #     ),

        #         Input("detail_graph","selectedNodeData"),

        # )

        # clientside_callback( ClientsideFunction(
        #         namespace='clientside',
        #         function_name='on_box_plot_click'
        #     ),
        #         Output("title","id",allow_duplicate=True),
        #             Input("activation_boxplot","hoverData"),
        #             Input("activation_boxplot","restyleData"),
        #     Input("activation_boxplot","figure"),
        #     State("box_plots_stats","data"),
        #                 prevent_initial_call=True,allow_duplicate=True

        # )
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
                        download_canvas_image(document.querySelector("#overview_graph canvas:nth-of-type(3)"),"overview.png");
                        break;
                    case "mono_graph":
                        download_canvas_image(document.querySelector("#mono_graph canvas:nth-of-type(3)"),"mono_signature.png",document.getElementById("mono_canvas"));
                        break;
                    case "detail_graph":
                        download_canvas_image(document.querySelector("#detail_graph canvas:nth-of-type(3)"),"multi_signature.png",document.getElementById("multi_canvas"));
                        break;
                    case "box":
                        const plot = document.querySelector("#activation_boxplot div.js-plotly-plot");
                        download_plotly_image(plot,"boxplot.png");
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
        # @callback(
        #     Output("box_plot_filter", "is_open"),
        #     Input("box_plot_filter_button", "n_clicks"),
        #     State("box_plot_filter", "is_open"),prevent_initial_call=True
        # )
        # def toggle_collapse(n_clicks, is_open):
        #     if n_clicks is not None:
        #         return not is_open
        #     return dash.no_update

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