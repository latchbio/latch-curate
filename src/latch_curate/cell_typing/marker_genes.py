from anndata import AnnData

marker_genes = {
    "T cell/NK cell": [
        "CD3D",
        "CD3E",
        "CD8A",
        "CD4", 
        "GNLY",
        "NKG7",
        "KLRD1"
        "IL7R",
        "KLRB1"
    ],
    "B cell/plasma cell": [
        "JCHAIN",
        "CD19",  
        "MS4A1", 
        "CD79A"  
    ],
    "mononuclear phagocyte": [
        "CD14",  
        "CD68",  
        "CSF1R", 
        "HLA-DRA",
        "HLA-DRB1",
        "CX3CL1",
        "TMEM119"
        "ITGAX"  
    ],
    "granulocyte": [
        "S100A8",
        "S100A9",
        "ELANE", 
        "MPO",   
        "CCR3",  
        "FCER1A",
        "IL5RA", 
        "SIGLEC8"
    ],
    "mast cell": [
        "TPSAB1",
        "CMA1",  
        "KIT"    
    ],
    "endothelial cell": [
        "CDH5",
        "VWF",
        "PECAM1"
    ],
    "fibroblast": [
        "COL1A1",
        "COL1A2",
        "DCN"
    ],
    "smooth muscle cell": [
        "ACTA2",
        "TAGLN",
        "MYH11"
    ],
    "epithelial cell": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT5",
        "KRT14",
        "CDH1"
    ],
    "neuron": [
        "RBFOX3",
        "MAP2"
    ],
    "astrocyte": [
        "GFAP",
        "AQP4"
    ],
    "oligodendrocyte": [
        "MOG",
        "MBP"
    ],
    "erythrocyte": [
        "HBB",
        "HBA1",
        "HBA2"
    ],
    "platelet": [
        "PPBP"
    ]
}

def remove_absent_symbols(adata: AnnData, marker_dict: dict[str, list[str]]):
    valid_genes = set(adata.var['gene_symbols'])

    new_dict = {}
    for celltype, gene_list in marker_dict.items():
        filtered = [g for g in gene_list if g in valid_genes]
        new_dict[celltype] = filtered

    return new_dict
