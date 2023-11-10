from pathlib import Path
import sys
import tempfile
import warnings
import random
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import seaborn as sns
from tqdm import tqdm


sc.set_figure_params(figsize=(4, 4), frameon=False)
torch.set_float32_matmul_precision("high")

warnings.simplefilter(action="ignore", category=FutureWarning)

def usage():
    print('Usage: python script.py [input_h5ad] [output_h5ad]')

def train_scvi_reference(adata = None):
    if adata is None:
        raise ValueError
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", layer="counts")
    model_scvi_ref = scvi.model.SCVI(
        adata,n_layers=2, n_latent=30, gene_likelihood="nb"
    )
    # Start training scVI
    model_scvi_ref.train(max_epochs = 500)
    model_scvi_ref.save(adata_path.stem + '_scVI_model', overwrite=True)
    # Start training scANVI
    model_scanvi_ref = scvi.model.SCANVI.from_scvi_model(model_scvi_ref, unlabeled_category="Unknown",labels_key = 'cell_type')
    model_scanvi_ref.train(max_epochs=50)
    model_scanvi_ref.save(adata_path.stem + '_scANVI_model', overwrite=True)
    



if __name__ == '__main__':
    try:
        adata_path = Path(sys.argv[1])      
        ref_adata = sc.read_h5ad(adata_path)
        ref_adata.layers['counts'] = ref_adata.X
        train_scvi_reference(ref_adata)
        
    except ValueError:
        usage()
