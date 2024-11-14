import numpy 
import scanpy as sc
import AnnData
import squidpy as sq
import matplotlib as plt
from typing import Union, Optional
from sklearn.preprocessing import StandardScaler


def compute_ratio_score(adata:AnnData,
                        gene_group_1: list = None,
                        gene_group_2: list = None,
                        key: str = None,
                        scale : bool = False,
                        center_score : bool = False,
                        score_type : str = 'mean_ratio',
                        in_place : bool = True,
                        add_name = None,
                        verbose : bool = True) -> AnnData:
    if gene_group_1 is None and gene_group_2 is None:
        raise ValueError('No gene sets provided to compute ratio')
    
    g1_counts = adata.X[gene_group_1]
    g2_counts = adata.X[gene_group_2]
    
    if scale:
        # Option to parse the sklearn scaler directly 
        scaler = StandardScaler()
        g1_counts = scaler.fit_transform(g1_counts)
        g2_counts = scaler.fit_transform(g2_counts)
    
    if key is not None:
        keys = adata.obs[key].index
        values = adata.obs[key].to_numpy()
    else :
        keys = adata.obs[key].index
        values = [0] * len(keys)
    
    g1_score = compute_score(g1_counts, keys, values, score_type)
    g2_score = compute_score(g2_counts, keys, values, score_type)
    
    return adata