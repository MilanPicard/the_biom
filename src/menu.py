from dash_extensions.enrich import Dash, dcc, html, Input, Output,callback
import dash_bootstrap_components as dbc


def menu(diseases,comparisons,diseases_cmap,comparisons_cmap,genes):
    colors_str = dict([(d,f"rgba({','.join([str(i*255) for i in diseases_cmap[d][:3]])},{diseases_cmap[d][3]})") for d in diseases])
    disease_checklist = dcc.Checklist(
        [{"label":html.Div(d,style={
            'color':colors_str[d]
            }),"value":d,
        #   
            } for d in diseases],
        [],id="disease_filter",labelStyle={"display": "flex", "alignItems": "center"})
    comparisons_checklist = dcc.Checklist(comparisons,comparisons,id="comparisons_filter")
    return [
            dbc.Accordion([
                dbc.AccordionItem([
                    disease_checklist
                ],title="Diseases",item_id="disease_accordion"),
                dbc.AccordionItem([
                    comparisons_checklist
                ],title="Comparisons"),
                dbc.AccordionItem([
                    dcc.Dropdown(options= [{"label":"None","value":"None"}]+[{"label":f"{i}","title" : f"{j} signatures","value":i} for i,j in genes.items()],value="None",id="genes_menu_select"),
                    html.Div(id="selected_genes_div",style={"minHeight":"2em","borderStyle":"ridge"}),
                    dcc.Store(id="selected_genes_store",data={"selected":[]})
                ],title="Genes",item_id="gene_accordion")
            ],
            start_collapsed=False,
            always_open=True,
            active_item=["gene_accordion","disease_accordion"],
            style={"height":"100%","overflowY": "scroll"}
            ),
            dcc.Store(id="selected_genes")]