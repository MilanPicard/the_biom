import json
import time
import numpy as np
import pandas as pd
import pyarrow as pa
from collections import defaultdict
from functools import reduce
from itertools import chain
from matplotlib import cm
class DataManager(object):
    _instance = None
    def __new__(cls,signature_file,expression_file,pathway_file):
        if cls._instance is None:
            cls._instance = super(DataManager, cls).__new__(cls)
            cls._instance.signature_file = signature_file
            cls._instance.expression_file = expression_file
            cls._instance.signatures = pd.read_csv(signature_file,dtype=pd.StringDtype())
            cls._instance.signatures["Signature"] = cls._instance.signatures["Signature"].str.split(";").astype(pd.ArrowDtype(pa.list_(pa.string())))
            merged = cls._instance.signatures.groupby(["Cancer","Comparison"]).agg(Filter=pd.NamedAgg(column="Filter",aggfunc=lambda x: "Merge"),Signature=pd.NamedAgg(column="Signature",aggfunc=lambda x:set([g for sign in x for g in sign])),gProfiler=pd.NamedAgg(column="gProfiler",aggfunc=lambda x: "??"))
            cls._instance.signatures=pd.concat([cls._instance.signatures,merged.reset_index()])
            cls._instance.signatures["id"] = cls._instance.signatures["Cancer"]+"_"+cls._instance.signatures["Comparison"]+"_"+cls._instance.signatures["Filter"]
            # cls._instance.signatures["Comparison"] = cls._instance.signatures["Comparison"].str.split("_").astype(pd.ArrowDtype(pa.list_(pa.string()))).list[0::2]
            cls._instance.exploded = cls._instance.signatures.explode("Signature").rename(columns={"Signature":"EnsemblID"})
            cls._instance.activation_data = pd.read_csv(cls._instance.expression_file,delimiter=";",dtype=defaultdict(lambda :float,Cancer=str,Stage=str))
            cls._instance.activation_data["box_category"] = cls._instance.activation_data["Cancer"]+"_"+cls._instance.activation_data["Stage"]
            cls._instance.pathways = pd.read_csv(pathway_file,dtype={"EnsemblID":str,"UniProtID":str,"PathwayStId":str,"PathwayDisplayName":str,"PathwayReactomeLink":str},converters={"GeneSymbolID":lambda s:list(map(lambda symbol:symbol.replace("'",""),s[1:-1].split(",")))},keep_default_na=False)
            cls._instance.pathways["GeneSymbolID"] = cls._instance.pathways["GeneSymbolID"].astype(pd.ArrowDtype(pa.list_(pa.string())))
            symbols = cls._instance.pathways.filter(["EnsemblID","GeneSymbolID"]).groupby("EnsemblID").agg(lambda a:a.iloc[0])
            value_counts = pd.DataFrame({"counts":cls._instance.exploded["EnsemblID"].value_counts()})
            symbols = value_counts.join(symbols,how="left")
            symbols.reset_index(inplace=True)
            symbols["GeneSymbolID"] = symbols["GeneSymbolID"].fillna(symbols["EnsemblID"].transform(lambda a:[a]))
            # genes.fillna(inplace=True)
            cls._instance.symbols = symbols.set_index("EnsemblID").filter(["GeneSymbolID"])

            # with open(pathway_file,"rt") as f:               
                # cls._instance.pathways= json.load(f)
        return cls._instance
    def get_filters(self):
        return self.signatures["Filter"].unique().tolist()
    def get_diseases(self):
        return self.signatures["Cancer"].unique().tolist()
    # def get_stages(self):
    #     stages:list = self.signatures["Comparison"].list.flatten().unique().tolist()
    #     return stages
    def get_comparisons(self):
        stages:list = self.signatures["Comparison"].unique().tolist()
        return stages
    def get_signatures_id(self):
        return self.signatures["id"].to_list()
    def get_genes(self,selected_filter="Merge"):
        value_counts = pd.DataFrame({"counts":self.exploded[self.exploded["Filter"]==selected_filter]["EnsemblID"].value_counts()})
        genes = value_counts.join(self.symbols,how="left")
        genes.reset_index(inplace=True)
        genes["GeneSymbolID"] = genes["GeneSymbolID"].fillna(genes["EnsemblID"].transform(lambda a:[a]))
        # genes.fillna(inplace=True)
        genes = genes.set_index("EnsemblID")

        # symbols = self.pathways[self.pathways["EnsemblID"].isin(value_counts.keys())].filter(["EnsemblID","GeneSymbolID"]).groupby("EnsemblID").agg(test)
        # print(symbols,self.exploded[self.exploded["Filter"]==selected_filter]["Signature"].value_counts(),symbols.join(self.exploded[self.exploded["Filter"]==selected_filter]["Signature"].value_counts(),how="right").fillna())
        return genes.to_dict("index")
    def get_symbol(self,gene):
        return self.symbols.loc[gene]["GeneSymbolID"]
    
    def get_signatures_intersections(self,disease_filter=[],comparisons_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]
        if(len(disease_filter)>0):
            exploded = exploded[exploded["Cancer"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        
        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        grouped_by_gene = exploded.filter(["EnsemblID","id"]).groupby("EnsemblID").agg(lambda a:a.values)
        grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        grouped_by_gene = grouped_by_gene.reset_index()[grouped_by_gene["id"].list.len()>1]
        
        intersections = dict()
        for item in grouped_by_gene.itertuples():
            genes  = sorted(item.id)
            for src in range(len(genes)-1):
                for tgt in range(src+1,len(genes)):
                    if (genes[src],genes[tgt]) in intersections:
                        intersections[(genes[src],genes[tgt])].append(item.EnsemblID)
                    else:
                        intersections[(genes[src],genes[tgt])]= [item.EnsemblID]
        return intersections,self.signatures[self.signatures["Filter"]==selected_filter].to_dict('records')
    def get_genes_intersections(self,disease_filter=[],comparisons_filter=[],id_filter=[],gene_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]

        if(len(disease_filter)>0):
            exploded = exploded.loc[exploded["Cancer"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        if 'None' in gene_filter:
            gene_filter.remove('None')
        if(len(id_filter)>0 or len(gene_filter)>0):
            genes = []
            if len(gene_filter)>0 :
                genes = gene_filter
            if(len(id_filter)>0):
                genes = genes+ exploded[exploded["id"].isin(id_filter)]["EnsemblID"].unique().tolist()
            id_filter = exploded[exploded["EnsemblID"].isin(genes)]["id"].unique()
            exploded = exploded[exploded["id"].isin(id_filter)]

        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        if(exploded.empty):
            return dict(),None
        exploded = exploded.filter(["EnsemblID","id"])
        # grouped_by_id = exploded.groupby("id").agg(lambda a:a.values)
        # grouped_by_id["EnsemblID"] = grouped_by_id["EnsemblID"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # grouped_by_id_multiple = grouped_by_id.reset_index()[grouped_by_id["EnsemblID"].list.len()>1]
        intersections = dict()
        # for item in grouped_by_id_multiple.itertuples():
        #     for src in range(len(item.EnsemblID)-1):
        #         for tgt in range(src+1,len(item.EnsemblID)):
        #             if (item.EnsemblID[src],item.EnsemblID[tgt]) in intersections:
        #                 intersections[(item.EnsemblID[src],item.EnsemblID[tgt])].append(item.id)
        #             else:
        #                 intersections[(item.EnsemblID[src],item.EnsemblID[tgt])]= [item.id]
        #             assert item.EnsemblID[src] in exploded["EnsemblID"].unique() and item.EnsemblID[tgt] in exploded["EnsemblID"].unique()
        grouped_by_id_nb=exploded.groupby("EnsemblID").agg(lambda a:a.size)
        # grouped_by_id_nb["id"] = grouped_by_id_nb["id"].astype(int)
        grouped_by_id_nb = grouped_by_id_nb.rename(columns={"id":"size"})
        grouped_by_id_nb= grouped_by_id_nb.join(exploded.groupby("EnsemblID").agg(lambda a:a.values))
        grouped_by_id_nb["id"] = grouped_by_id_nb["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        grouped_by_id_nb =grouped_by_id_nb.reset_index()

        # grouped_by_id["id"] = grouped_by_id["id"].list.len()
        return intersections,grouped_by_id_nb.rename(columns={"EnsemblID":"gene","id":"id"})
    def get_genes_intersection_data(self,disease_filter=[],comparisons_filter=[],id_filter=[],gene_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]

        if(len(disease_filter)>0):
            exploded = exploded.loc[exploded["Cancer"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        if 'None' in gene_filter:
            gene_filter.remove('None')
        if(len(id_filter)>0 or len(gene_filter)>0):
            genes = []
            if len(gene_filter)>0 :
                genes = gene_filter
            if(len(id_filter)>0):
                genes = genes+ exploded[exploded["id"].isin(id_filter)]["EnsemblID"].unique().tolist()
            id_filter = exploded[exploded["EnsemblID"].isin(genes)]["id"].unique()
            exploded = exploded[exploded["id"].isin(id_filter)]

        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        return exploded
    
    def get_activations(self,genes,diseases):
        genes  = [genes] if not isinstance(genes,list) else genes
        activations = self.activation_data
        if len(diseases)>0:
            activations = activations[self.activation_data["Cancer"].isin(diseases)]
        activations = activations[genes+["Cancer",'Stage',"box_category"]]
        return activations
    
    def get_pathways(self,genes):
        genes = set(genes)
        # filtered = [p for p in self.pathways if any([g in genes for g in p["genes"]])]
        # pathways = {g:[] for g in genes}
        # for p in self.pathways:
        #     for g in p["genes"]:
        #         if g in genes:
        #             pathways[g].append(p)
        # return pathways
        return self.pathways[self.pathways['EnsemblID'].isin(genes)]
    def get_disease_cmap(self):
        unique_diseases = self.signatures["Cancer"].unique()
        colormap = cm.get_cmap("tab10",len(unique_diseases))
        return dict([(unique_diseases[i],colormap(i)) for i in range(len(unique_diseases))])
    # def get_stage_cmap(self):
    #     unique_stages = self.signatures["Comparison"].explode().unique()
    #     colormap = cm.get_cmap("tab10",len(unique_stages))
    #     return dict([(unique_stages[i],colormap(i)) for i in range(len(unique_stages))])
    def get_comparison_cmap(self):
        unique_comparisons = self.signatures["Comparison"].unique().tolist()
        colormap = cm.get_cmap("tab10",len(unique_comparisons))
        return dict([(unique_comparisons[i],colormap(i)) for i in range(len(unique_comparisons))])
    def get_diseases_from_genes(self,genes):
        exploded = self.exploded[self.exploded["EnsemblID"].isin(genes)]
        return exploded["Cancer"].unique().tolist()
    def get_comparisons_from_genes(self,genes):
        exploded = self.exploded[self.exploded["EnsemblID"].isin(genes)]
        return exploded["Comparison"].unique().tolist()


