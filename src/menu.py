# from dash_extensions.enrich import Dash, dcc, html, Input, Output,callback
from dash import Dash, dcc, html, Input, Output,callback
import dash_bootstrap_components as dbc


def menu(diseases,comparisons,diseases_cmap,comparisons_cmap,genes,filters,pathways):
    
    filter_select = dcc.Dropdown(filters,"Merge",id="filters_dropdown")
    disease_accordions = dbc.Accordion([
                  dbc.AccordionItem(
                    # dcc.Checklist(options=[{'label':c,'disabled':True,'value':c} for c in sorted(comparisons)],value=[],id={"type":"signature_checkbox","disease":d},inline=True)
                    [dbc.Button(c,id={"type":"signature_checkbox","disease":d,"comparison":c},disabled=True,active=False,color="signature") for c in sorted(comparisons)]

                        ,title=d,item_id=f"signature_{d}_accordion",id={"type":"disease_accordion","disease":d}) for d in diseases
                        ],
                        start_collapsed=False,
                        always_open=True,
                        active_item=[],class_name="signature_accordion"
                        )
    return [
            dbc.Accordion([
                dbc.AccordionItem(disease_accordions,title="Signatures",item_id="signature_accordion"),
                # old_disease(diseases,diseases_cmap),
                # old_comp(comparisons),
                diseases_buttons(diseases,diseases_cmap),
                comparisons_buttons(sorted(comparisons)),
                dbc.AccordionItem([
                    filter_select
                ],title="Filter"),
                dbc.AccordionItem([
                    dcc.Dropdown(options= [{"label":"None","value":"None"}]+[{"label":f"{' '.join(j['GeneSymbolID'])}","title" : f"{j['counts']} signatures","value":i} for i,j in genes.items()],value="None",id="genes_menu_select"),
                    html.Div(id="selected_genes_div",style={"minHeight":"2em","borderStyle":"ridge"}),
                    dcc.Store(id="selected_genes_store",data={"selected":[],"from_pathways":{"ids":[],"genes":[]}})
                ],title="Genes",item_id="gene_accordion"),
                dbc.AccordionItem([
                    dcc.Dropdown(options= [{"label":"None","value":"None"}]+[{"label":j['PathwayDisplayName'],"title" : f"{j['counts']} genes","value":i} for i,j in pathways.items()],value="None",id="pathways_menu_select",optionHeight=70),
                    html.Div(id="selected_pathways_div",style={"minHeight":"2em","borderStyle":"ridge"})
                ],title="Pathways",item_id="pathway_accordion"),
                dbc.AccordionItem([
                    dcc.Dropdown(options= [{"label":"Overview","value":"overview"} ,{"label":"mono","value":"mono_graph"} ,{"label":"multi","value":"detail_graph"},{"label":"box","value":"box"}  ],value="None",id="exportImage",optionHeight=70),
                    html.Button("Export Image",id="export_image_btn"),
                    html.Button("Export Json",id="export_json_btn")
                ],title="Exports",item_id="exports_accordion")
            ],
            start_collapsed=False,
            always_open=True,
            active_item=["signature_accordion"],
            style={"height":"100%","overflowY": "auto"}
            ),
            dcc.Store(id="selected_genes"),
                 dcc.Store(id="selected_signatures",data={}),
            ]
def diseases_buttons(diseases,diseases_cmap):
    colors_str = dict([(d,f"rgba({','.join([str(i*255) for i in diseases_cmap[d][:3]])},{diseases_cmap[d][3]})") for d in diseases])
    return check_all_buttons(diseases,"disease",colors_str)
    # return dbc.AccordionItem([
    #     dbc.ButtonGroup([
    #         html.Div(d,style={
    #         'color':colors_str[d]
    #         }),
    #         dbc.Button("All",id={"type":"check_diseases","disease":d},outline=True,color="primary"),
    #         dbc.Button("None",id={"type":"uncheck_diseases","disease":d},outline=True,color="primary",disabled=True)
    #     ],size="sm",class_name="menu_all_none_btn_group") for d in diseases]
    # ,title="Diseases")
def comparisons_buttons(comparisons):
    return check_all_buttons(comparisons,"comparison")
def check_all_buttons(elems,elem_type:str,colors_str=None):
    title = f"Select by {elem_type.capitalize()}s"
    return dbc.AccordionItem([
        dbc.ButtonGroup([
            html.Div(d,style=({'color':colors_str[d]} if colors_str is not None else None)),
            dbc.Button("All",id={"type":f"check_{elem_type}s",elem_type:d},outline=True,color="primary"),
            dbc.Button("None",id={"type":f"uncheck_{elem_type}s",elem_type:d},outline=True,color="primary",disabled=True)
        ],size="sm",class_name="menu_all_none_btn_group") for d in elems]
    ,title=title)
def old_comp(comparisons):
    comparisons_checklist = dcc.Checklist(sorted(comparisons),comparisons,id="comparisons_filter")
    return dbc.AccordionItem([
                    comparisons_checklist,
                    dbc.ButtonGroup([
                        dbc.Button("All Comparisons",id="check_all_comparisons",outline=True,color="primary",disabled=True),
                        dbc.Button("No Comparisons",id="check_no_comparisons",outline=True,color="primary")]
                    ,size="sm")
                ],title="Comparisons")

def old_disease(diseases,diseases_cmap):
    colors_str = dict([(d,f"rgba({','.join([str(i*255) for i in diseases_cmap[d][:3]])},{diseases_cmap[d][3]})") for d in diseases])
    disease_checklist = dcc.Checklist(
        [{"label":html.Div(d,style={
            'color':colors_str[d]
            }),"value":d,
        #   
            } for d in diseases],
        [],id="disease_filter",labelStyle={"display": "flex", "alignItems": "center"})
    return dbc.AccordionItem([
                    disease_checklist,
                    dbc.ButtonGroup([
                        dbc.Button("All Diseases",id="check_all_diseases",outline=True,color="primary"),
                        dbc.Button("No Diseases",id="check_no_diseases",outline=True,color="primary",disabled=True)]
                    ,size="sm")
                ],title="Diseases",item_id="disease_accordion")