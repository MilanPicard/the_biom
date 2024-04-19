import json
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
            cls._instance.signatures["id"] = cls._instance.signatures["Disease"]+"_"+cls._instance.signatures["Comparison"]
            cls._instance.signatures["Signature"] = cls._instance.signatures["Signature"].str.split(";").astype(pd.ArrowDtype(pa.list_(pa.string())))
            cls._instance.signatures["Comparison"] = cls._instance.signatures["Comparison"].str.split("_").astype(pd.ArrowDtype(pa.list_(pa.string()))).list[0::2]
            cls._instance.exploded = cls._instance.signatures.explode("Signature")
            cls._instance.activation_data = pd.read_csv(cls._instance.expression_file,delimiter=";",dtype=defaultdict(lambda :float,Disease=str,Stage=str))
            cls._instance.activation_data["box_category"] = cls._instance.activation_data["Disease"]+"_"+cls._instance.activation_data["Stage"]
            with open(pathway_file,"rt") as f:
                cls._instance.pathways= json.load(f)
        return cls._instance

    def get_diseases(self):
        return self.signatures["Disease"].unique().tolist()
    def get_stages(self):
        stages:list = self.signatures["Comparison"].list.flatten().unique().tolist()
        return stages
    def get_comparisons(self):
        stages:list = self.signatures["Comparison"].transform(lambda a:"_".join(a)).unique().tolist()
        return stages
    def get_signatures_id(self):
        return self.signatures["id"].to_list()
    def get_genes(self):
        return self.exploded["Signature"].value_counts().to_dict()
    def get_signatures_intersections(self,disease_filter=[],comparisons_filter=[]):
        exploded = self.exploded
        if(len(disease_filter)>0):
            exploded = exploded[exploded["Disease"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        
        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        grouped_by_gene = exploded.filter(["Signature","id"]).groupby("Signature").agg(lambda a:a.values)
        grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        grouped_by_gene = grouped_by_gene.reset_index()[grouped_by_gene["id"].list.len()>1]
        
        intersections = dict()
        for item in grouped_by_gene.itertuples():
            genes  = sorted(item.id)
            for src in range(len(genes)-1):
                for tgt in range(src+1,len(genes)):
                    if (genes[src],genes[tgt]) in intersections:
                        intersections[(genes[src],genes[tgt])].append(item.Signature)
                    else:
                        intersections[(genes[src],genes[tgt])]= [item.Signature]
        return intersections,self.signatures.to_dict('records')
    def get_genes_intersections(self,disease_filter=[],comparisons_filter=[],id_filter=[],gene_filter=[]):
        exploded = self.exploded

        if(len(disease_filter)>0):
            exploded = exploded.loc[exploded["Disease"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        if 'None' in gene_filter:
            gene_filter.remove('None')
        if(len(id_filter)>0 or len(gene_filter)>0):
            genes = []
            if len(gene_filter)>0 :
                genes = gene_filter
            if(len(id_filter)>0):
                genes = genes+ exploded[exploded["id"].isin(id_filter)]["Signature"].unique().tolist()
            id_filter = exploded[exploded["Signature"].isin(genes)]["id"].unique()
            exploded = exploded[exploded["id"].isin(id_filter)]

        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        exploded = exploded.filter(["Signature","id"])
        grouped_by_id = exploded.groupby("id").agg(lambda a:a.values)
        grouped_by_id["Signature"] = grouped_by_id["Signature"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        grouped_by_id_multiple = grouped_by_id.reset_index()[grouped_by_id["Signature"].list.len()>1]
        intersections = dict()
        for item in grouped_by_id_multiple.itertuples():
            for src in range(len(item.Signature)-1):
                for tgt in range(src+1,len(item.Signature)):
                    if (item.id[src],item.id[tgt]) in intersections:
                        intersections[(item.Signature[src],item.Signature[tgt])].append(item.id)
                    else:
                        intersections[(item.Signature[src],item.Signature[tgt])]= [item.id]
                    assert item.Signature[src] in exploded["Signature"].unique() and item.Signature[tgt] in exploded["Signature"].unique()
        grouped_by_id_nb=exploded.groupby("Signature").agg(lambda a:a.size)
        # grouped_by_id_nb["id"] = grouped_by_id_nb["id"].astype(int)
        grouped_by_id_nb = grouped_by_id_nb.rename(columns={"id":"size"})
        grouped_by_id_nb= grouped_by_id_nb.join(exploded.groupby("Signature").agg(lambda a:a.values))
        grouped_by_id_nb["id"] = grouped_by_id_nb["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        grouped_by_id_nb =grouped_by_id_nb.reset_index()
        # grouped_by_id["id"] = grouped_by_id["id"].list.len()
        return intersections,grouped_by_id_nb.rename(columns={"Signature":"gene","id":"id"})
    
    def get_activations(self,genes,diseases):
        genes  = [genes] if not isinstance(genes,list) else genes
        activations = self.activation_data
        if len(diseases)>0:
            activations = activations[self.activation_data["Disease"].isin(diseases)]
        activations = activations[genes+['Disease','Stage',"box_category"]]
        return activations
    
    def get_pathways(self,genes):
        return dict([(g,self.pathways[g]) for g in genes])
    def get_disease_cmap(self):
        unique_diseases = self.signatures["Disease"].unique()
        colormap = cm.get_cmap("tab10",len(unique_diseases))
        return dict([(unique_diseases[i],colormap(i)) for i in range(len(unique_diseases))])
    def get_stage_cmap(self):
        unique_stages = self.signatures["Comparison"].explode().unique()
        colormap = cm.get_cmap("tab10",len(unique_stages))
        return dict([(unique_stages[i],colormap(i)) for i in range(len(unique_stages))])
    def get_comparison_cmap(self):
        unique_comparisons = self.signatures["Comparison"].transform(lambda a:"_".join(a)).unique().tolist()
        colormap = cm.get_cmap("tab10",len(unique_comparisons))
        return dict([(unique_comparisons[i],colormap(i)) for i in range(len(unique_comparisons))])
    def get_diseases_from_genes(self,genes):
        exploded = self.exploded[self.exploded["Signature"].isin(genes)]
        return exploded["Disease"].unique().tolist()
    def get_comparisons_from_genes(self,genes):
        exploded = self.exploded[self.exploded["Signature"].isin(genes)]
        return exploded["Comparison"].transform(lambda a:"_".join(a)).unique().tolist()


