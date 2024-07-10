from dash_extensions.enrich import Dash, dcc, html, Input, Output,callback
import dash_bootstrap_components as dbc


def menu(diseases,comparisons,diseases_cmap,comparisons_cmap,genes,filters,pathways):
    colors_str = dict([(d,f"rgba({','.join([str(i*255) for i in diseases_cmap[d][:3]])},{diseases_cmap[d][3]})") for d in diseases])
    disease_checklist = dcc.Checklist(
        [{"label":html.Div(d,style={
            'color':colors_str[d]
            }),"value":d,
        #   
            } for d in diseases],
        [],id="disease_filter",labelStyle={"display": "flex", "alignItems": "center"})
    comparisons_checklist = dcc.Checklist(comparisons,comparisons,id="comparisons_filter")
    filter_select = dcc.Dropdown(filters,"Merge",id="filters_dropdown")
    return [
            dbc.Accordion([
                dbc.AccordionItem([
                    disease_checklist
                ],title="Diseases",item_id="disease_accordion"),
                dbc.AccordionItem([
                    comparisons_checklist
                ],title="Comparisons"),
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
                    dcc.Dropdown(options= [{"label":"Overview","value":"overview"} ,{"label":"mono","value":"mono"} ,{"label":"multi","value":"multi"},{"label":"box","value":"box"}  ],value="None",id="exportImage",optionHeight=70),
                    html.Button("Export Image",id="export_image_btn")
                ],title="Exports",item_id="exports_accordion")
            ],
            start_collapsed=False,
            always_open=True,
            active_item=["gene_accordion","disease_accordion"],
            style={"height":"100%","overflowY": "scroll"}
            ),
            dcc.Store(id="selected_genes")]