import scanpy as sc
import scvi
import anndata
import numpy as np
import pandas as pd

adata = anndata.read("mammal_InN_migration.final.p9.h5ad")
adata.X = adata.X.tocsr()
adata.X = adata.X.astype(np.int64)
adata.layers["counts"] = adata.X.copy()

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="labsite",
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="orig.ident",
    categorical_covariate_keys=["labsite","transAge","lobe","subtype"]
)
model = scvi.model.SCVI(adata)
model.train(max_epochs=100)

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="lineage",
    unlabeled_category="Unknown",
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
