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
#### gene simulation functions:
```python
import scCube
from scCube import scCube
from scCube.visualization import *
from scCube.utils import *

model = scCube()
# pre-processing
sc_adata = model.pre_process(sc_data=sc_data, 
                             sc_meta=sc_meta,
                             is_normalized=False)
```
**Parameters**

**sc_data**: _DataFrame_

&emsp;DataFrame of input data

**sc_meta**: _DataFrame_

&emsp;DataFrame of input meta

**is_normalized**: _bool, default: `False`_

&emsp;Whether the input data is normalized or not. If `is_normalized=False`, the input data will be normalized by scCube first.

****

```python
generate_sc_meta, generate_sc_data = model.train_vae_and_generate_cell(sc_adata=sc_adata,
                                                                       celltype_key='Cell_type',
                                                                       cell_key='Cell',
                                                                       target_num=None,
                                                                       batch_size=512,
                                                                       epoch_num=10000,
                                                                       lr=0.0001,
                                                                       hidden_size=128,
                                                                       save_model=True,
                                                                       save_path=save_path,
                                                                       project_name=model_name,
                                                                       used_device='cuda:0'
                                                                       )
```
**Parameters**

**sc_adata**: _AnnData_

&emsp;AnnData of pre-processed data 

**celltype_key**: _str_

&emsp;The column name of `cell types` or `domain` in meta

**cell_key**: _str_

&emsp;The column name of `cell` in meta

**target_num**: _Optional[dict], default: `None`_

&emsp;Target number of cells to generate, if `target_num=None`, generate cells by the proportion of cell types of the input data.

**batch_size**: _int, default: `512`_

&emsp;Batch size of training

**epoch_num**: _int, default: `3500`_

&emsp;Epoch number of training

**lr**: _float, default: `0.0001`_

&emsp;Learning reta of training

**hidden_size**: _int, default: `128`_

&emsp;Hidden size of VAE model

**save_model**: _bool, default: `True`_

&emsp;Whether save trained VAE model or not

**save_path**: _str_

&emsp;The save path

**project_name**: _str_

&emsp;The name of trained VAE model

**used_device**: str, default: `cuda:0`_

&emsp;Device name, `cpu` or `cuda

****

```python
generate_sc_meta, generate_sc_data = model.load_vae_and_generate_cell(self,
                                                                      sc_adata: AnnData,
                                                                      celltype_key: str,
                                                                      cell_key: str,
                                                                      target_num: Optional[dict] = None,
                                                                      hidden_size: int = 128,
                                                                      load_path: str = '',
                                                                      used_device: str = 'cuda:0'
                                                                      )
```
**Parameters**

**sc_adata**: _AnnData_

&emsp;AnnData of pre-processed data 

**celltype_key**: _str_

&emsp;The column name of `cell types` or `domain` in meta












