## API Reference
### Input data format

1. Data file: a `.csv` file with genes as rows and one sample as column

|<img width=40/> <img width=40/>|<img width=25/>Cell1<img width=25/>|<img width=25/>Cell2<img width=25/>|<img width=25/>Cell3<img width=25/>|<img width=35/>...<img width=34/>|<img width=25/>CellN<img width=25/>| 
| :-----: | :-----: | :-----: | :-----: | :-----: | :-----: | 
| Gene1 | 1 | 2 | 1 | ... | 0 |
| Gene2 | 4 | 1 | 0 | ... | 4 |
| ... | ... | ... | ... | ... | ... |
| GeneN | 0 | 0 | 2 | ... | 0 |

****  
2. Meta file: a `.csv` file with cell/spot ID and celltype/domain annotation columns
   * The column containing cell ID should be named `Cell` 
   * the column containing the labels should be named `Cell_type` 

|<img width=90/> <img width=90/>|<img width=80/>Cell<img width=81/>|<img width=75/>Cell_type<img width=75/>|
| :-----: | :-----: | :-----: |
| Cell1 | Cell1 | T cell |
| Cell2 | Cell2 | B cell |
| ... | ... | ... |
| CellN | CellN | Monocyte |

****

### Parameter description
* Decompose bulk transcriptomics data into single-cell transcriptomics data:
```python
from bulk2space import Bulk2Space
model = Bulk2Space()

# Decompose bulk transcriptomics data into single-cell transcriptomics data
generate_sc_meta, generate_sc_data = model.train_vae_and_generate(
    input_bulk_path,
    input_sc_data_path,
    input_sc_meta_path,
    input_st_data_path,
    input_st_meta_path,
    ratio_num=1,
    top_marker_num=500,
    gpu=0,
    batch_size=512,
    learning_rate=1e-4,
    hidden_size=256,
    epoch_num=5000,
    vae_save_dir='save_model',
    vae_save_name='vae',
    generate_save_dir='output',
    generate_save_name='output')
```

| Parameter | Description | Default Value |
| --- | --- | --- |
| input_bulk_path | Path to bulk-seq data files (.csv) | None |
| input_sc_data_path | Path to scRNA-seq data files (.csv) | None |
| input_sc_meta_path | Path to scRNA-seq annotation files (.csv) | None |
| input_st_data_path | Path to ST data files (.csv) | None |
| input_st_meta_path | Path to ST metadata files (.csv) | None |
| ratio_num | The multiples of the number of cells of generated scRNA-seq data | (int) `1` |
| top_marker_num | The number of marker genes of each celltype used | (int) `500`  |
| gpu | The GPU ID. Use cpu if `--gpu < 0` | (int) `0` |
| batch_size | The batch size for β-VAE model training | (int) `512` |
| learning_rate | The learning rate for β-VAE model training | (float) `0.0001` |
| hidden_size | The hidden size of β-VAE model | (int) `256` |
| epoch_num | The epoch number for β-VAE model training | (int) `5000` |
| vae_save_dir | Path to save the trained β-VAE model | (str) `save_model` |
| vae_save_name | File name of the trained β-VAE model | (str) `vae` |
| generate_save_dir | Path to save the generated scRNA-seq data | (str) `output` |
| generate_save_name | File name of the generated scRNA-seq data | (str) `output` |







