# NOTE: The code expects the 'Disease' column in the signatures CSV to be renamed to 'Cancer' for internal consistency.
import json
import os
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
            # Rename 'Disease' column to 'Cancer' for compatibility
            if 'Disease' in cls._instance.signatures.columns:
                cls._instance.signatures = cls._instance.signatures.rename(columns={'Disease': 'Cancer'})
            cls._instance.signatures["Signature"] = cls._instance.signatures["Signature"].str.split(";").astype(pd.ArrowDtype(pa.list_(pa.string())))
            if "Merge" not in cls._instance.signatures["Filter"].unique():
                merged = cls._instance.signatures.groupby(["Cancer","Comparison"]).agg(
                    Filter=pd.NamedAgg(column="Filter",aggfunc=lambda x: "Merge"),
                    Signature=pd.NamedAgg(column="Signature",aggfunc=lambda x:set([g for sign in x for g in sign])),
                    gProfiler=pd.NamedAgg(column="gProfiler",aggfunc=lambda x: "??")
                )
                cls._instance.signatures=pd.concat([cls._instance.signatures,merged.reset_index()])
            cls._instance.signatures["id"] = cls._instance.signatures["Cancer"]+"_"+cls._instance.signatures["Comparison"]+"_"+cls._instance.signatures["Filter"]
            # cls._instance.signatures["Comparison"] = cls._instance.signatures["Comparison"].str.split("_").astype(pd.ArrowDtype(pa.list_(pa.string()))).list[0::2]
            cls._instance.exploded = cls._instance.signatures.explode("Signature").rename(columns={"Signature":"EnsemblID"})
            cls._instance.activation_data = pd.read_csv(cls._instance.expression_file,delimiter=",",dtype=defaultdict(lambda :float,Cancer=str,Stage=str,ID=str))
            cls._instance.activation_data["box_category"] = cls._instance.activation_data["Cancer"]+"_"+cls._instance.activation_data["Stage"]
            cls._instance.pathways = pd.read_csv(pathway_file,dtype={"EnsemblID":pd.ArrowDtype(pa.string()),"UniProtID":pd.ArrowDtype(pa.string()),"PathwayStId":pd.ArrowDtype(pa.string()),"PathwayDisplayName":pd.ArrowDtype(pa.string()),"PathwayReactomeLink":pd.ArrowDtype(pa.string())},converters={"GeneSymbolID":lambda s:s.split(";")},keep_default_na=False)
            cls._instance.pathways["GeneSymbolID"] = cls._instance.pathways["GeneSymbolID"].astype(pd.ArrowDtype(pa.list_(pa.string())))
            symbols = cls._instance.pathways[cls._instance.pathways["GeneSymbolID"]!=""].filter(["EnsemblID","GeneSymbolID"]).groupby("EnsemblID").agg(lambda a:a.iloc[0])
            value_counts = pd.DataFrame({"counts":cls._instance.exploded["EnsemblID"].value_counts()})
            symbols = value_counts.join(symbols,how="left")
            symbols.reset_index(inplace=True)
            symbols["GeneSymbolID"] = symbols["GeneSymbolID"].fillna(symbols["EnsemblID"].transform(lambda a:[a]))
            # genes.fillna(inplace=True)
            cls._instance.pathways= cls._instance.pathways[cls._instance.pathways["PathwayStId"]!=""]
            
            cls._instance.pathway_labels = cls._instance.pathways.filter(["PathwayStId","PathwayDisplayName"]).groupby("PathwayStId").agg(lambda a: a.iloc[0])["PathwayDisplayName"]
            cls._instance.symbols = symbols.set_index("EnsemblID").filter(["GeneSymbolID"])

            # with open(pathway_file,"rt") as f:               
                # cls._instance.pathways= json.load(f)
        return cls._instance
    @staticmethod
    def get_instance():
        if DataManager._instance is None:
            if ("THE_BIOM_SIGNATURES" in os.environ and "THE_BIOM_EXPRESSIONS" in os.environ and "THE_BIOM_PATHWAYS" in os.environ):
                signatures = os.environ["THE_BIOM_SIGNATURES"].strip()
                expressions = os.environ["THE_BIOM_EXPRESSIONS"].strip()
                pathways = os.environ["THE_BIOM_PATHWAYS"].strip()
            else:
                signatures =os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","signatures","THe_Biom_DEV_dataset.csv")
                expressions = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","activations","fake_data.csv")
                pathways = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","pathways","fake_pathways")
            return DataManager(signatures,expressions,pathways)
        return DataManager._instance
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
    def get_signatures_id(self,selected_filter=None):
        if selected_filter is None:
            return self.signatures["id"].to_list()
        else:
            return self.signatures.loc[self.signatures["Filter"]==selected_filter]["id"].to_list()

    def get_genes(self,selected_filter="Merge",pathway=None):
        value_counts = pd.DataFrame({"counts":self.exploded[self.exploded["Filter"]==selected_filter]["EnsemblID"].value_counts()})
        genes = value_counts.join(self.symbols,how="left")
        genes.reset_index(inplace=True)
        genes["GeneSymbolID"] = genes["GeneSymbolID"].fillna(genes["EnsemblID"].transform(lambda a:[a]))
        # genes.fillna(inplace=True)
        genes = genes.set_index("EnsemblID")

        # symbols = self.pathways[self.pathways["EnsemblID"].isin(value_counts.keys())].filter(["EnsemblID","GeneSymbolID"]).groupby("EnsemblID").agg(test)
        # print(symbols,self.exploded[self.exploded["Filter"]==selected_filter]["Signature"].value_counts(),symbols.join(self.exploded[self.exploded["Filter"]==selected_filter]["Signature"].value_counts(),how="right").fillna())
        if pathway is not None:
            genes = genes.loc[self.pathways[self.pathways["PathwayStId"]==pathway]["EnsemblID"]]
        return genes.to_dict("index")
    def get_pathway_label(self,pathway):
        return self.pathways[self.pathways["PathwayStId"]==pathway]["PathwayDisplayName"].iloc[0]
    
    def get_symbol(self,gene):
        if isinstance(gene,str):
            return self.symbols.loc[gene]["GeneSymbolID"]
        return self.symbols.loc[list(gene)]["GeneSymbolID"]

    def get_intersection_gprofiler_link(self, signature1, signature2, common_genes):
        """Generate a gProfiler link for the intersection of two signatures.
        
        Args:
            signature1 (str): First signature ID
            signature2 (str): Second signature ID
            common_genes (list): List of genes in the intersection
            
        Returns:
            str: gProfiler link for the intersection
        """
        # Get the cancer type from either signature (they should be the same)
        cancer = signature1.split("_")[0]
        
        # Create a description for the intersection
        stage1 = signature1.split("_")[1].split("vs")[1]
        stage2 = signature2.split("_")[1].split("vs")[1]
        description = f"{cancer}_Intersection_{stage1}_vs_{stage2}"
        
        # Create the gProfiler link with the common genes
        genes_str = ";".join(common_genes)
        return f"https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query={genes_str}&description={description}"

    def get_signatures_intersections(self,disease_filter=[],comparisons_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]
        if(len(disease_filter)>0):
            exploded = exploded[exploded["Cancer"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        
        grouped_by_gene = exploded.filter(["EnsemblID","id"]).groupby("EnsemblID").agg(lambda a:a.values)
        grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        pathways = self.pathways.filter(["EnsemblID","PathwayStId"])
        pathways = pathways[pathways["EnsemblID"].isin(grouped_by_gene.index)].groupby("EnsemblID").agg(lambda a:a.unique())
        def agg_signature(a):
            pd.set_option("max_colwidth",100)
            pairs = set()
            
            for i in range(len(a)):
                for j in range(i+1,len(a)):
                    sign_i=set(a.iloc[i]).difference(set(a.iloc[j]))
                    sign_j=set(a.iloc[j]).difference(set(a.iloc[i]))

                    for s1 in sign_i:
                        for s2 in sign_j:
                            # if s1!=s2:
                                if s1<s2:
                                    pairs.add("__".join((s1,s2)))
                                else:
                                    pairs.add("__".join((s2,s1)))
            pd.reset_option("max_colwidth")

            return pairs
        pathways = grouped_by_gene.join(pathways,how="inner").explode("PathwayStId").filter(["id","PathwayStId"]).groupby("PathwayStId").agg(lambda a: agg_signature(a))
        pathways = pathways.assign(count=lambda x: x.id.apply(len))
        pathways = pathways[pathways["count"]>0].filter(["id"])
        pathways_edges = pathways.filter(["id"]).explode("id").reset_index().groupby("id").agg(lambda a:a.unique())
        pathways_edges["PathwayStId"]=pathways_edges["PathwayStId"].apply(lambda pids:[f"{self.pathway_labels.loc[pid]}({pid})" for pid in pids])
        pathways_edges = pd.DataFrame.from_dict({"PathwayStId":[]})
        grouped_by_gene = grouped_by_gene.reset_index()[grouped_by_gene["id"].list.len()>1]
        
        intersections = dict()
        intersection_links = dict()  # New dictionary to store gProfiler links for intersections
        
        for item in grouped_by_gene.itertuples():
            genes = sorted(item.id)
            for src in range(len(genes)-1):
                for tgt in range(src+1,len(genes)):
                    key = (genes[src],genes[tgt])
                    if key in intersections:
                        intersections[key].append(item.EnsemblID)
                    else:
                        intersections[key] = [item.EnsemblID]
                        # Generate gProfiler link for this intersection
                        intersection_links[key] = self.get_intersection_gprofiler_link(
                            genes[src], 
                            genes[tgt], 
                            [item.EnsemblID]  # This will be updated as more genes are added
                        )
        
        # Update intersection links with complete gene lists
        for key, genes in intersections.items():
            intersection_links[key] = self.get_intersection_gprofiler_link(
                key[0], 
                key[1], 
                genes
            )
            
        return intersections, self.signatures[self.signatures["Filter"]==selected_filter].to_dict('records'), pathways_edges, intersection_links
    def get_disease_comp_from_genes(self,gene_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]
        id_filter = exploded[exploded["EnsemblID"].isin(gene_filter)]["id"].unique().tolist()
        exploded = exploded[exploded["id"].isin(id_filter)]
        return exploded["Cancer"].unique(),exploded["Comparison"].unique()

    def get_genes_intersections(self,disease_filter=[],comparisons_filter=[],id_filter=[],gene_filter=[],selected_filter="Merge"):
        exploded = self.exploded[self.exploded["Filter"]==selected_filter]

        if(len(disease_filter)>0):
            exploded = exploded.loc[exploded["Cancer"].isin(disease_filter)]
        if(len(comparisons_filter)>0):
            exploded = exploded[exploded["Comparison"].isin(comparisons_filter)]
        if 'None' in gene_filter:
            gene_filter.remove('None')
        if(len(gene_filter)>0):
            # genes = []
            # if(len(id_filter)>0):
                # genes = genes+ exploded[exploded["id"].isin(id_filter)]["EnsemblID"].unique().tolist()
            id_filter += exploded[exploded["EnsemblID"].isin(gene_filter)]["id"].unique().tolist()
        if(len(id_filter)>0 or len(gene_filter)>0):
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
        if(len(gene_filter)>0):
            genes = gene_filter
            # id_filter = exploded[exploded["EnsemblID"].isin(genes)]["id"].unique()

            # if(len(id_filter)>0):
                # genes = genes+ exploded[exploded["id"].isin(id_filter)]["EnsemblID"].unique().tolist()
            id_filter += exploded[exploded["EnsemblID"].isin(genes)]["id"].unique().tolist()
        if(len(id_filter)>0 or len(gene_filter)>0):
            exploded = exploded[exploded["id"].isin(id_filter)]

        # grouped_by_gene = exploded.groupby("Signature").agg(lambda a:a.to_list() if len(a) >1 else a)
        # grouped_by_gene["id"] = grouped_by_gene["id"].astype(pd.ArrowDtype(pa.list_(pa.string())))
        # nbsign = exploded.groupby("Signature").agg(lambda x: x.size)
        # grouped_by_gene = grouped_by_gene[nbsign>1]

        # exploded["id"] = exploded["id"].transform(lambda a:[a])
        return exploded
    
    def get_activations(self,genes,diseases,comparisons):
        genes  = [genes] if not isinstance(genes,list) else genes
        activations = self.activation_data
        if len(diseases)>0:
            activations = activations[self.activation_data["Cancer"].isin(diseases)]
        if len(comparisons)>0:
            activations = activations[self.activation_data["Stage"].isin([c.split("vs")[1]  for c in comparisons])]
        activations = activations[genes+["Cancer",'Stage',"box_category"]]
        return activations

    def get_pathways(self,genes,filter=None):
        if(len(genes)==0):
            filtered = self.exploded[self.exploded["Filter"]==filter]
            filtered2 = filtered.filter(["EnsemblID","id"]).groupby("EnsemblID").agg(lambda a: a.unique()).astype({"id":pd.ArrowDtype(pa.list_(pa.string()))})
            # print(filtered2)
            # print(self.pathways.join(filtered2,on="EnsemblID",how="right").dtypes)
            # print(self.pathways.join(filtered2,on="EnsemblID").filter(["PathwayStId","PathwayDisplayName","EnsemblID","id"]).groupby("PathwayStId").agg({
            filtered3 = self.pathways.join(filtered2,on="EnsemblID",how="right").filter(["PathwayStId","PathwayDisplayName","EnsemblID","id"]).groupby("PathwayStId").agg({
                "PathwayDisplayName":(lambda a:a.iloc[0]),
                "EnsemblID":(lambda a:a.unique()),
                "id":(lambda a: a)
            }).astype({"EnsemblID":pd.ArrowDtype(pa.list_(pa.string()))})#,"id":pd.ArrowDtype(pa.list_(pa.list_(pa.string())))})

            def keep(l):
                for i in range(len(l)):
                    set_i = set(l[i])
                    for j in range(i+1,len(l)):
                        if set_i.isdisjoint(set(l[j])):
                            return True
                return False
            filtered3 = filtered3.assign(to_keep=(lambda x: x.id.reset_index().assign(to_keep=lambda y:y.id.apply(keep)).set_index("PathwayStId")["to_keep"]))
            filtered3 = filtered3.loc[filtered3["to_keep"]]
            filtered3 = filtered3.assign(counts=(lambda x: x.id.reset_index().assign(counts=lambda y : y["id"].map(lambda v : len(set([j for i in v for j in i ])))).set_index("PathwayStId")["counts"]))
            filtered3.sort_values("counts",inplace=True,ascending=False)
            
            return filtered3.filter(["PathwayDisplayName","PathwayStId","counts","EnsemblID"]).to_dict("index")
        genes = set(genes)
        # filtered = [p for p in self.pathways if any([g in genes for g in p["genes"]])]
        # pathways = {g:[] for g in genes}
        # for p in self.pathways:
        #     for g in p["genes"]:
        #         if g in genes:
        #             pathways[g].append(p)
        # return pathways
        p=  self.pathways[self.pathways['EnsemblID'].isin(genes)]
        p=p.transform({
            'EnsemblID':lambda col:col.apply(lambda a: [a]).astype(pd.ArrowDtype(pa.list_(pa.string()))),
            'UniProtID':lambda col:col.apply(lambda a: [a]).astype(pd.ArrowDtype(pa.list_(pa.string()))),
            'PathwayStId':lambda col:col,
            'PathwayDisplayName':lambda col:col,
            'PathwayReactomeLink':lambda col:col,
                       })
        p=p.filter(["EnsemblID", "UniProtID" , "PathwayStId","PathwayDisplayName","PathwayReactomeLink"]).groupby("PathwayStId").agg({
            "EnsemblID":lambda a:  [ i  for i in a.list.flatten()],
            "UniProtID":lambda a: a,
            "PathwayDisplayName":lambda a: a.iloc[0],
            "PathwayReactomeLink":lambda a: a.iloc[0],
        })
        return p 
    def get_disease_cmap(self):
        # Stable color mapping:
        # - Keep fixed colors for the original six cancers in their historical order
        #   UCEC, BRCA, LIHC, LUAD, KIRC, HNSC mapped to tab10 indices 0..5
        # - Assign additional cancers the next available distinct colors without
        #   changing the colors of already-mapped cancers
        present_diseases = list(self.signatures["Cancer"].unique())

        fixed_order = ["UCEC", "BRCA", "LIHC", "LUAD", "KIRC", "HNSC"]
        tab10 = cm.get_cmap("tab10")

        # Build a pool of candidate colors (RGBA), starting with tab10 then others
        candidate_colors = []
        # tab10 first (10 colors)
        candidate_colors.extend([tab10(i) for i in range(10)])

        # Extend the pool with other qualitative palettes to cover many cancers
        extra_maps = [
            cm.get_cmap("tab20"),
            cm.get_cmap("tab20b"),
            cm.get_cmap("tab20c"),
            cm.get_cmap("Set3"),
            cm.get_cmap("Dark2"),
            cm.get_cmap("Accent"),
            cm.get_cmap("Paired"),
        ]
        for m in extra_maps:
            colors = getattr(m, "colors", None)
            if colors is None:
                # Sample all entries if the colormap is continuous
                colors = [m(i) for i in range(m.N)]
            for c in colors:
                # Ensure RGBA tuples
                if len(c) == 3:
                    candidate_colors.append((c[0], c[1], c[2], 1.0))
                else:
                    candidate_colors.append((c[0], c[1], c[2], c[3]))

        # Deduplicate while preserving order
        seen = set()
        deduped_colors = []
        for c in candidate_colors:
            key = tuple(round(x, 6) for x in c)  # robust float comparison
            if key not in seen:
                seen.add(key)
                deduped_colors.append(c)

        # Assign fixed colors for the original six cancers when present.
        # IMPORTANT: Use a colormap sized to the number of base diseases present
        # so their colors exactly match the previous behavior when only those are shown.
        cmap = {}
        used = set()
        base_present = [d for d in fixed_order if d in present_diseases]
        if len(base_present) > 0:
            base_cmap = cm.get_cmap("tab10", len(base_present))
            for i, disease in enumerate(base_present):
                color = base_cmap(i)
                cmap[disease] = color
                used.add(tuple(round(x, 6) for x in color))

        # Assign remaining diseases from the pool without reusing colors
        color_index = 0
        for disease in present_diseases:
            if disease in cmap:
                continue
            # Find the next unused color
            while color_index < len(deduped_colors):
                candidate = deduped_colors[color_index]
                color_index += 1
                key = tuple(round(x, 6) for x in candidate)
                if key not in used:
                    cmap[disease] = candidate
                    used.add(key)
                    break

        return cmap
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

    def get_pathway_signatures(self, pathway_id, selected_filter="Merge"):
        """Get all signatures associated with genes in a specific pathway"""
        pathway_genes = self.pathways[self.pathways["PathwayStId"] == pathway_id]["EnsemblID"].tolist()
        if not pathway_genes:
            return []
        
        # Get signatures that contain any of these genes
        pathway_signatures = self.exploded[
            (self.exploded["Filter"] == selected_filter) & 
            (self.exploded["EnsemblID"].isin(pathway_genes))
        ]["id"].unique().tolist()
        
        return pathway_signatures


