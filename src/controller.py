from dash_extensions.enrich import Dash, html, Input, Output, State, callback,ctx,ALL
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash_extensions.enrich import clientside_callback,ClientsideFunction, Input, Output
from dash import Patch
import dash
import numpy as np
import plotly.graph_objs as go
import pandas as pd
import plotly.express as px
import overview as ov
import nx_layout
import scipy
from scipy.spatial import ConvexHull
import detail_graph
import utils
class Controller(object):
    _instance = None
    def __new__(cls,dm):
        if cls._instance is None:
            cls._instance = super(Controller, cls).__new__(cls)
            cls._instance.dm = dm
        return cls._instance
@callback(
        Output("genes_menu_select","options"),
        Input("filters_dropdown","value")
)
def update_selectable_genes(selected_filter):
    return [{"label":"None","value":"None"}]+[{"label":f"{' '.join(j['GeneSymbolID'])}","title" : f"{j['counts']} signatures","value":i} for i,j in Controller._instance.dm.get_genes(selected_filter).items()]

@callback(
        Output('selected_genes_store','data'),
        Input("genes_menu_select","value"),
        Input({"type":"selected_gene_button","gene":ALL},"n_clicks"),
        Input("detail_graph","selectedNodeData"),
        Input("filters_dropdown","value"),
        State('selected_genes_store','data'),prevent_initial_call=True
)
def add_remove_gene(menu_select,button,fromGraph,selected_filter,current):
    match ctx.triggered_id:
        case "genes_menu_select":
            if menu_select is not None and menu_select != "None" and menu_select not in current["selected"]:
                current["selected"].append(menu_select)
        case "detail_graph":
            if fromGraph is not None:
                for gene in fromGraph:
                    if gene["id"] not in current["selected"]:
                        current["selected"].append(gene["id"])
        case "filters_dropdown":
            current["selected"]=[]
        case _:
            current["selected"].remove(ctx.triggered_id["gene"])
    return current        
@callback(
    Output("selected_genes_div","children"),
    Input("selected_genes_store","data"),
    State("selected_genes_div","children"),
    prevent_initial_call=False
)
def update_genes_buttons(store_data,cur_children):
    patched = Patch()
    if cur_children is None:
        cur_children = []
    n=len(cur_children)
    if len(store_data["selected"])==0 and n>0:
        patched.clear()
    else:
        if n!=len(store_data["selected"]):
            color_scheme = detail_graph.get_color_scheme(store_data["selected"])
            if(len(store_data["selected"])<len(cur_children)):
                found=False
                for i in range(len(store_data["selected"])):
                    if cur_children[i]["props"]["id"]["gene"]!=store_data["selected"][i]:
                        del patched[i]
                        del cur_children[i]
                        found = True
                        break
                if not found:
                    del patched[len(cur_children)-1]
                    del cur_children[-1]

                n-=1
            else:
                btn_style = {}
                btn = html.Button(Controller._instance.dm.get_symbol(store_data["selected"][-1]),id={"type":"selected_gene_button","gene":store_data["selected"][-1]},className="btn",style=btn_style)
                patched.append(btn)
                n+=1
            for i,g in enumerate(store_data["selected"]):
                tc = utils.get_text_color(color_scheme[i%len(color_scheme)])
                patched[i]["props"]["style"]={
                    "--bs-btn-bg":color_scheme[i%len(color_scheme)],
                    "--bs-btn-hover-bg":color_scheme[i%len(color_scheme)],
                    "--bs-btn-hover-color":tc,
                    "--bs-btn-color":tc}
    return patched          


@callback(
        Output('overview_graph','elements'),
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
    print(ctx.triggered_id)
    if cur_elems[0]["data"]["Filter"]!=selected_filter:
        cur_elems = ov.get_elements(Controller._instance.dm,selected_filter=selected_filter)
        print("change")
    else:
        classes = ["highlight","half_highlight"]
        for i in cur_elems:

            if "Cancer" in i["data"]:
                if "classes" not in i:
                    print(i)
                    i["classes"]=" ".join([i["data"]["Cancer"],"highlight"])
                if i['data']["Cancer"] in diseases and i['data']["Comparison"] in comparisons_filter:
                    c = "highlight"
                    if( selected_genes_store is not None):
                        c = "half_highlight"
                        if len(set(selected_genes_store["selected"]).intersection(i["data"]["Signature"]))!=0:
                            c = "highlight"
                    i["classes"]=utils.switch_class(i["classes"],[c],classes)
                else:
                    i["classes"]=utils.switch_class(i["classes"],[],classes)
            else:
                if "classes" not in i:
                    i["classes"]=""

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
    return cur_elems

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
            prevent_initial_call=True
          )
def update_menu_disease_filter_data(genes):
    return Controller._instance.dm.get_diseases_from_genes(genes["selected"])
@callback(
        Output('comparisons_filter','value'),
        Input("selected_genes_store","data"),
        State('comparisons_filter','value'),
            prevent_initial_call=True
          )
def update_menu_comparison_filter_data(genes,cur_state):
    if len(genes["selected"])>0:
        return Controller._instance.dm.get_comparisons_from_genes(genes["selected"])
    return cur_state

@callback(Output('detail_graph','elements'),
        Output('detail_graph','stylesheet'),
        Output('detail_graph','layout'),
        Output('detail_graph_pos','data'),
        #   Input('overview_graph','selectedNodeData'),
        Input('disease_filter','value'),
        Input('comparisons_filter','value'),
        Input('data_overview_selected','value'),
        Input('selected_genes_store','data'),
        Input("fake_graph_size","data"),
        State('detail_graph','elements'),
        State("detail_graph_pos","data"),            
        State('detail_graph','stylesheet'),
            prevent_initial_call=True
          )
def display_detail_graph(diseases,comparisons,signatures,menu_genes,fake_graph_size,existing_elements,detail_pos_store,current_stylesheets):
    if ctx.triggered_id=="fake_graph_size" :
        if (fake_graph_size is None or "just_redraw" not in fake_graph_size or not fake_graph_size["just_redraw"]):
            raise dash.exceptions.PreventUpdate()
        else:
            print("fake_graph_size",fake_graph_size)
            return detail_graph.redraw(existing_elements,detail_pos_store,1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"],current_stylesheets)
    if(all([len(diseases)==0 or len(comparisons)==0 ,signatures is None or signatures ==""])):
        return [],[],{"name":"preset"},{}
    if len(diseases)!=0 or len(signatures)!=0:
        if diseases is None:
            diseases =""
        if signatures is None:
            signatures = ""
        return detail_graph.display_detail_graph(list(filter(lambda a: len(a)>0,diseases)),list(filter(lambda a: len(a)>0,signatures.split(";"))),menu_genes["selected"],existing_elements,detail_pos_store if detail_pos_store is not None else dict(),1 if fake_graph_size is None or "AR" not in fake_graph_size else fake_graph_size["AR"])
    else:
        return existing_elements,[],{"name":"preset"},{}

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
    Input("data_menu_selected","value"),
    Input("data_overview_selected","value"),
    # State("overview_graph","selectedNodeData" ),
    Input("selected_genes_store","data"),
        State("overview_graph","elements"),
        State("overview_graph","stylesheet"),
        prevent_initial_call=False
)
def update_box_plot(menu_selected_diseases,overview_selected,menu_selected,overview_elements,overview_stylesheets):
    items = []
    if menu_selected is not None and len(menu_selected["selected"])!=0:
        items += menu_selected["selected"]
    stylesheets =overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(Controller._instance.dm)
    stylesheets = [s for s in stylesheets if not(s["selector"].startswith("edge#"))]
    if(len(items)>0):
        g = items[0]
        diseases = menu_selected_diseases.split(";")
        if overview_selected is not None:
            overview_selected = overview_selected.split(";")
            diseases = diseases + list(set([i.split("_")[0] for i in overview_selected]))
        selected_patient_and_genes = Controller._instance.dm.get_activations(items,diseases).sort_values(["box_category"])
        box_categories = sorted(pd.unique(selected_patient_and_genes["box_category"]).tolist())
        symbols = list(map(lambda s : " ".join(s),Controller._instance.dm.get_symbol(items).to_list()))
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
        if(len(items)==1):
            print(selected_patient_and_genes,(symbols),items)
            box = px.box(selected_patient_and_genes,x="box_category",y=symbols[0],color_discrete_sequence=detail_graph.get_color_scheme(items),labels={"box_category":""})
            if(box.layout.margin.t is not None and box.layout.margin.t>20):
                box.layout.margin.t=20
            return box,"visible_plot",stylesheets

        else:
            dfs = []
            color_scheme=detail_graph.get_color_scheme(items)
            for i in range(len(items)):
                df = selected_patient_and_genes.filter(("box_category",symbols[i]))
                df = df.rename({symbols[i]:"expression"},axis=1)
                df["gene"] = symbols[i]
                dfs.append(df)
            df = pd.concat(dfs)
            box = px.box(df,x="box_category",y="expression",color="gene",color_discrete_sequence=color_scheme,labels={"box_category":""})
            if(box.layout.margin.t is not None and box.layout.margin.t>20):
                box.layout.margin.t=20
            return box ,"visible_plot",stylesheets
        
    else:
        return go.Figure(data=[
            ]),"hidden_plot",stylesheets
    
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

#         selected_patient_and_genes = Controller._instance.dm.get_activations(g,[disease]).sort_values(by=["Stage",g])
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
        """
function activate_tooltip(mouseoverNodeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction){
    console.log(mouseoverNodeData);
    if(mouseoverNodeData==null){
    
    //TODO keep cur if close else reset
        //return [cur_children,cur_show,cur_bbox,cur_direction];
        return [[],false,{},"right"];
    }else{
        var direction="left";
        var elem = undefined;

        if( mouseoverNodeData["id"] in pos_store_data){
            elem=pos_store_data[mouseoverNodeData["id"]];
        }else{
            for(var i of elements){
                if( i["data"]["id"] == mouseoverNodeData["id"]){
                    elem = i;
                    break;
                }
            }
        }
        var detail_graph = document.getElementById('detail_graph');
        var detail_resize_span = document.getElementById('detail_resize_span');
        //console.log("update_width_start",dash_clientside.callback_context.triggered_id,Object.assign({},fake_graph_size,{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':detail_resize_span.clientWidth,"AR":elem.clientWidth/elem.clientHeight}))
        var x = 0;
        var y=0;
        var width=0;
        var height=0;
        var wratio=detail_graph.clientWidth/(extent["x2"]-extent["x1"]);
        var hratio=detail_graph.clientHeight/(extent["y2"]-extent["y1"]);
        x = (elem["position"]["x"]-extent["x1"]) *wratio+detail_resize_span.clientWidth;
        y = (elem["position"]["y"]-extent["y1"]) *hratio
        direction = (x/detail_graph.clientWidth>0.5)?"left":"right";
        if( "weight" in mouseoverNodeData){
            width = mouseoverNodeData["weight"];
            height = mouseoverNodeData["weight"];
        }else{
            var polygon = [];
            var zeros = [];
            var y_crosses = [];
            for(var i of stylesheets){
                if(i["selector"] == "node#"+mouseoverNodeData["id"]){
                    width = parseFloat(i["style"]["width"]);
                    height = parseFloat(i["style"]["height"]);
                    var pos = i["style"]["shapePolygonPoints"].split(" ").map(parseFloat);
                    for(var p =0;p<pos.length;p+=2){
                        polygon.push([pos[p],pos[p+1]]);
                        if(pos[p+1]==0.0){
                            zeros.push(p/2);
                            y_crosses.push(pos[p+1]);
                        }
                    }

                    break;
                }
            }
            if(y_crosses.length<2){
                for(var i =0; i<polygon.length;i++){
                    if(polygon[i][0]!=0 && polygon[(i+1)%polygon.length][0]!=0){
                        if(polygon[i][1]*polygon[(i+1)%polygon.length][1]<0){
                            y_crosses.push(
                                polygon[i][0]+(Math.abs(polygon[i][1])/(Math.abs(polygon[i][1])+Math.abs(polygon[(i+1)%polygon.length][1])))*(polygon[(i+1)%polygon.length][0]-polygon[i][0])
                            );
                        }
                    }
                }
            }
            if(y_crosses.length==2){
                if(direction=="right"){
                    x+=0.5*y_crosses.reduce((a,b) => Math.max(a,b),-1)*width*wratio;
                }else{
                    x+=0.5*y_crosses.reduce((a,b) => Math.min(a,b),1)*width*wratio;
                }
                width=0;
                height=0;
            }else{
                console.log("y_crosses",y_crosses);
            }
        }
        return [mouseoverNodeData["tooltip_content"],true,{"x0":x-width*wratio/2,"y0":y-height*hratio/2,"x1":x+width*wratio/2,"y1":y+height*hratio/2},direction];
    }
}
""",
    Output("detail_graph_tooltip","children"),
    Output("detail_graph_tooltip","show"),
    Output("detail_graph_tooltip","bbox"),
    Output("detail_graph_tooltip","direction"),
#     # Output("detail_graph_tooltip","bbox"),
    Input("detail_graph","mouseoverNodeData"),
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
    
# clientside_callback(ClientsideFunction(
#     namespace="tooltip",
#     function_name="set_tooltip"),
#     Output("dummy_div","children"),
#     Input("detail_graph","mouseoverNodeData"),prevent_initial_call=True
# )

clientside_callback(
        """
function update_width_start(n1,n2,e_width,e_height,state,detail_graph_pos,dw,fake_graph_size){
    var e = dash_clientside.callback_context.triggered_id=="full_col"?e_height:e_width;
    if(dash_clientside.callback_context.triggered_id===undefined){
        var elem = document.getElementById('detail_graph');
        var span = document.getElementById('detail_resize_span');
        console.log("update_width_start",dash_clientside.callback_context.triggered_id,Object.assign({},fake_graph_size,{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':span.clientWidth,"AR":elem.clientWidth/elem.clientHeight}))

        return  Object.assign({},fake_graph_size,{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':span.clientWidth,"AR":elem.clientWidth/elem.clientHeight});
    }
    if(e["type"]=="click" && !e["isTrusted"]){
        if(fake_graph_size["width"]===undefined || fake_graph_size["height"]===undefined || Math.abs(fake_graph_size["width"]-document.getElementById("detail_graph").clientWidth)/document.getElementById("detail_graph").clientWidth>0.05 ||Math.abs(fake_graph_size["height"]-document.getElementById("detail_graph").clientHeight)/document.getElementById("detail_graph").clientHeight>0.05){
        fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
        fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
        fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
        fake_graph_size["just_redraw"]=true;
        document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true,"animate":false}).run();

        return fake_graph_size;
        }
        else
                            return dash_clientside.no_update;

    }
    if(e["type"]=="mousedown"){
        if(e.target.tagName=="SPAN"){
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
                var possibleTargets = ["overview_col","detail_col"];
                if((dash_clientside.callback_context.triggered_id==="full_col" && "second_row_div" == e.target.id) || (possibleTargets.includes(e.target.id) && dash_clientside.callback_context.triggered_id==="move_in_ov")){
                    if(!state["height"]["is_resizing"]){
                        state["height"]["is_resizing"]=true;
                        dash_clientside.set_props("resize_state", {"data":{"width":{"is_resizing":false},"height":{"is_resizing":true}}});
                        return dash_clientside.no_update;
                    }
                }else{
                    console.log(possibleTargets.includes(e.target.id),e.target.id,dash_clientside.callback_context.triggered_id==="full_col",dash_clientside.callback_context.triggered_id)
                }
        }
        return dash_clientside.no_update;
    }
    if(e["type"]=="mouseup"){
        if(state["width"]["is_resizing"]){
            state["width"]["is_resizing"]=false;
            dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
            document.querySelectorAll("#overview_col canvas, #detail_col canvas").forEach(c => c.style.cursor="auto");
            document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
                    fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
        fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
        fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
            return fake_graph_size;
        }
        if(state["height"]["is_resizing"]){
            state["height"]["is_resizing"]=false;
            dash_clientside.set_props("resize_state", {"data":{"width":state["width"],"height":state["height"]}});
            document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
                    fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
        fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
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
        document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true}).run();
        fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
        fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
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
                fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
                fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
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
                    fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
                    fake_graph_size["width_span"]=document.getElementById("detail_resize_span").clientWidth;
                    fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
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
            fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
            fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
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
            fake_graph_size["width"]=document.getElementById("detail_graph").clientWidth;
            fake_graph_size["height"]=document.getElementById("detail_graph").clientHeight;
            fake_graph_size["AR"]=fake_graph_size["width"]/fake_graph_size["height"];
            fake_graph_size["just_redraw"]=true;
            dash_clientside.set_props("fake_graph_size", fake_graph_size);
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