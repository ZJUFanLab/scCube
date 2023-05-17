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
#### gene expression simulation functions:
##### pre-processing
```python
import scCube
from scCube import scCube
from scCube.visualization import *
from scCube.utils import *

model = scCube()
sc_adata = model.pre_process(
    sc_data=sc_data, 
    sc_meta=sc_meta,
    is_normalized=False
    )
```
**Parameters**

**sc_data**: _DataFrame_

&emsp;DataFrame of input data

**sc_meta**: _DataFrame_

&emsp;DataFrame of input meta

**is_normalized**: _bool, default: `False`_

&emsp;Whether the input data is normalized or not. If `is_normalized=False`, the input data will be normalized by scCube first.

****

##### train vae model and generate gene expression
```python
generate_sc_meta, generate_sc_data = model.train_vae_and_generate_cell(
    sc_adata=sc_adata,
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

**used_device**: _str, default: `cuda:0`_

&emsp;Device name, `cpu` or `cuda

****

##### load VAE model and generate gene expression
```python
generate_sc_meta, generate_sc_data = model.load_vae_and_generate_cell(
    sc_adata=sc_adata,
    celltype_key='Cell_type',
    cell_key='Cell',
    target_num=None,
    hidden_size=128,
    load_path=load_path,
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

**hidden_size**: _int, default: `128`_

&emsp;Hidden size of VAE model

**load_path**: _str_

&emsp;The load path

**used_device**: _str, default: `cuda:0`_

&emsp;Device name, `cpu` or `cuda

****

#### spatial pattern simulation functions:
##### generate random spatial patterns with reference-free strategy
```python
generate_sc_data, generate_sc_meta, st_data, st_meta, st_index = model.generate_spatial_data_random(
    generate_sc_data=generate_sc_data,
    generate_sc_meta=generate_sc_meta,
    set_seed=False,
    seed=12345,
    spatial_cell_type=None,
    spatial_dim=2,
    spatial_size=30,
    delta=25,
    lamda=0.75,
    is_spot=False,
    platform='ST',
    gene_type='whole',
    min_cell=10,
    n_gene=None,
    n_cell=5,
    is_split=True,
    split_coord='point_z',
    slice_num=5,
    )
```
**Parameters**

**generate_sc_data**: _DataFrame_

&emsp;DataFrame of generated data

**generate_sc_meta**: _DataFrame_

&emsp;DataFrame of generated meta

**set_seed**: _bool, default: `False`_

&emsp;Whether to set seed for reproducible simulation

**seed**: _int, default: `12345`_

&emsp;The seed number

**spatial_cell_type**: _ Optional[list], default: `None`_

&emsp;The selected cell types with spatial patterns, if`spatial_cell_type=None`, all cell types would be assigned spatial patterns

**spatial_dim**: _int, default: `2`_

&emsp;The spatial dimensionality， `2` or `3`

**spatial_size**: _int, default: `30`_

&emsp;The scope for spatial autocorrelation function, the large values will take more running time

**delta**: _float, default: `25`_

&emsp;The larger value will tend to form spatial patterns with greater connectivity 

**lamda**: _float, default: `0.75`_

&emsp;The larger values will tend to form clearer spatial patterns

**is_spot**: _bool, default: `False`_

&emsp;True -- generate spot-based SRT data; False -- generate imaging-based SRT data

**platform**: _str, default: `ST`_

&emsp;Only works when `is_spot=True`, `ST` -- square neighborhood structure; `Visium` -- hexagonal neighborhood structure; `Slide` -- random neighborhood structure

**gene_type**: _str, default: `whole`_

&emsp;The type of genes to generate, `whole` -- the whole genes; `hvg` -- the highly variable genes; `marker` -- the marker genes of each cell type; `random` -- the randomly selected genes

**min_cell**: _int, default: `10`_

&emsp;Filter the genes expressed in fewer than `min_cell` cells before selected genes, only works when `gene_type='random', 'hvg', or 'marker'`

**n_gene**: _Optional[int], default: `None`_

&emsp;The number of genes to select, only works when `gene_type='random', 'hvg', or 'marker'`

**n_cell**: _int, default: `10`_

&emsp;The average number of cells per spot, only works when `is_spot=True`

**is_split**: _bool, default: `True`_

&emsp;Whether to spilt the 3D generated spatial patterns into a series of 2D spatial patterns, only works when `spatial_dim=3`

**split_coord**: _str, default: `point_z`_

&emsp;The name of split coordinate axis, only works when `spatial_dim=3` and `is_split=True`

**slice_num**: _int, default: `5`_

&emsp;The targeted number of 2D slices, only works when `spatial_dim=3` and `is_split=True`

****

##### generate spatial patterns with reference-based strategy
```python
generate_sc_data, generate_sc_meta = model.generate_spatial_data_reference(
    sc_adata=sc_adata,
    generate_sc_data=generate_sc_data,
    generate_sc_meta=generate_sc_meta,
    celltype_key='Cell_type',
    spatial_key=['x', 'y'],
    cost_metric='sqeuclidean'
    )
```
**Parameters**

**sc_adata**: _AnnData_

&emsp;AnnData of reference data

**generate_sc_data**: _DataFrame_

&emsp;DataFrame of generated data

**generate_sc_meta**: _DataFrame_

&emsp;DataFrame of generated meta

**celltype_key**: _str_

&emsp;The column name of `cell types` or `domain` in meta

**spatial_key**: _list_

&emsp;The column name of `spatial coordinates` in meta

**cost_metric**: _str, defalut: `sqeuclidean`_

&emsp;The cost distance between generate_sc_data and real_data, `sqeuclidean` by default. On numpy the function also accepts from the scipy.spatial.distance.cdist function : ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’, ‘yule’.

****

#### visualization functions:
##### Plot scatter plot of cell type spatial pattern
```python
p = plot_spatial_pattern_scatter(
    obj=generate_sc_meta,
    figwidth=8,
    figheight=8,
    dim=2,
    x="point_x",
    y="point_y",
    z=None,
    label=None,
    palette=None,
    colormap='rainbow',
    size=10,
    alpha=1,
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of generated meta

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**dim**: _int, defalut: `2`_

&emsp;Spatial dimensionality

**x**: _str, defalut: `point_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `point_y`_

&emsp;The name of column containing y coordinate

**z**: _Optional[str], default: `None`_

&emsp;The name of column containing z coordinate, only use when 'dim = 3'

**label**: _Optional[str], default: `None`_

&emsp;The name of column containing cell type information, if 'label=None', plot coordinates without cell type information only.

**palette**: _Optional[list], default: `None`_

&emsp;List of colors used, if 'palette=None', plot scatter plot with colormap colors

**colormap**: str, default: `rainbow`_

&emsp;The name of cmap

**size**: _float, default: `10`_

&emsp;The size of point

**alpha**: _float, default: `1`_

&emsp;The transparency of point

****

##### Plot density plot of cell type spatial pattern
```python
p = plot_spatial_pattern_density(
    obj=generate_sc_meta,
    figwidth=8,
    figheight=8,
    x="point_x",
    y="point_y",
    label="Cell_type",
    show_celltype=None,
    colormap='Blues',
    fill=True,
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of generated meta

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**x**: _str, defalut: `point_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `point_y`_

&emsp;The name of column containing y coordinate

**label**: _str, default: `Cell_type`_

&emsp;The name of column containing cell type information, if 'label=None', plot coordinates without cell type information only.

**show_celltype**: _Optional[str], default: `None`_

&emsp;The cell type selected to plot separately, if 'show_celltype=None', plot all cell type together

**colormap**: str, default: `Blues`_

&emsp;The name of cmap

**fill**: _bool, default: `True`_

&emsp;If 'fill=True', fill in the area between bivariate contours

****

##### Plot scatterpie plot of spot-based data
```python
p = plot_spot_scatterpie(
    obj=prop,
    figwidth=8,
    figheight=8,
    x="spot_x",
    y="spot_y",
    palette=None,
    colormap='rainbow',
    res=50,
    direction="+",
    start=0.0,
    size=100,
    edgecolor="none",
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of cell type proportion per spot

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**x**: _str, defalut: `spot_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `spot_y`_

&emsp;The name of column containing y coordinate

**palette**: _Optional[dict], default: `None`_

&emsp;Dict of color of each cell type, if 'palette == None', plot scatterpie plot with colormap colors

**colormap**: str, default: `rainbow`_

&emsp;The name of cmap

**res**: _int, default: `50`_

&emsp;Number of points around the circle

**direction**: _str, default: `+`_

&emsp;'+' for counter-clockwise, or '-' for clockwise

**start**: _flost, default: `0.0`_

&emsp;Starting position in radians

**size**: _float, default: `100`_

&emsp;The size of point

**edgecolor**: _str, default: `none`_

&emsp;The edge color of point

****

##### Plot scatter plot of proportion of selected cell type
```python
p = plot_spot_prop(
    obj=prop,
    figwidth=8,
    figheight=8,
    x="spot_x",
    y="spot_y",
    colormap='viridis',
    show_celltype= "",
    size=100,
    alpha=1,
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of cell type proportion per spot

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**x**: _str, defalut: `spot_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `spot_y`_

&emsp;The name of column containing y coordinate

**colormap**: str, default: `viridis`_

&emsp;The name of cmap

**show_celltype**: _Union[list, str]_

&emsp;The cell type selected to plot

**size**: _float, default: `100`_

&emsp;The size of point

**alpha**: _float, default: `1`_

&emsp;The transparency of point

****

##### Plot scatter plot of spatial expression pattern of selected gene
```python
p = plot_gene_scatter(
    data=generate_sc_data,
    obj=generate_sc_meta_new,
    figwidth=8,
    figheight=8,
    dim=2,
    label='Cell',
    normalize=True,
    x="point_x",
    y="point_y",
    z="point_z",
    colormap='viridis',
    show_gene: str = "",
    size=10,
    alpha=1,
    )
plt.show(p)
```
**Parameters**

**data**: _DataFrame_

&emsp;DataFrame of generate data

**obj**: _DataFrame_

&emsp;DataFrame of generate meta

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**dim**: _int, default: `2`_

&emsp;Spatial dimensionality

**label**: _str, default: `Cell`_

&emsp;The name of column containing cell/spot name

**normalize**: _bool, default: `True`_

&emsp;If 'normalize=True', normalizing expression value to [0, 1]

**x**: _str, defalut: `point_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `point_x`_

&emsp;The name of column containing y coordinate

**z**: _Optional[str], default: `None`_

&emsp;The name of column containing z coordinate, only use when 'dim = 3'

**colormap**: str, default: `viridis`_

&emsp;The name of cmap

**show_gene**: _str_

&emsp;The gene selected to plot

**size**: _float, default: `10`_

&emsp;The size of point

**alpha**: _float, default: `1`_

&emsp;The transparency of point

****

##### Plot histplot of spot-based data to investigate cell number per spot
```python
p = plot_gene_scatter(
    obj=st_index,
    figwidth=8,
    figheight=8,
    label='spot',
    n_bin=20
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of cell-spot index

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**label**: _str, default: `spot`_

&emsp;The name of column containing spot name

**n_bins**: int, default: `20`_

&emsp;The number of equal-width bins in the range

****

##### Plot 2d scatter plot of cell type spatial pattern of each slices from 3d data
```python
p = plot_slice_scatter(
    obj=generate_sc_meta,
    figwidth=8,
    figheight=8,
    x="point_x",
    y="point_y",
    label='Cell_type',
    palette=None,
    colormap='rainbow',
    size=10,
    alpha=1
    )
plt.show(p)
```
**Parameters**

**obj**: _DataFrame_

&emsp;DataFrame of generated meta

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**x**: _str, defalut: `point_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `point_y`_

&emsp;The name of column containing y coordinate

**label**: _str, default: `Cell_type`_

&emsp;The name of column containing cell type information

**palette**: _Optional[list], default: `None`_

&emsp;List of colors used, if 'palette=None', plot scatter plot with colormap colors

**colormap**: str, default: `rainbow`_

&emsp;The name of cmap

**size**: _float, default: `10`_

&emsp;The size of point

**alpha**: _float, default: `1`_

&emsp;The transparency of point

****

##### Plot 2d scatter plot of spatial expression pattern of selected gene of each slices from 3d data
```python
p = plot_slice_gene_scatter(
    data=generate_sc_data,
    obj=generate_sc_meta,
    figwidth=8,
    figheight=8,
    x="point_x",
    y="point_y",
    label='Cell_type',
    normalize=True,
    show_gene="",
    colormap='viridis',
    size=10,
    alpha=1
    )
plt.show(p)
```
**Parameters**

**data**: _DataFrame_

&emsp;DataFrame of generated data

**obj**: _DataFrame_

&emsp;DataFrame of generated meta

**figwidth**: _float, default: `8`_

&emsp;Figure width

**figheight**: _float, default: `8`_

&emsp;Figure height

**x**: _str, defalut: `point_x`_

&emsp;The name of column containing x coordinate

**y**: _str, defalut: `point_y`_

&emsp;The name of column containing y coordinate

**label**: _str, default: `Cell_type`_

&emsp;The name of column containing cell type information

**normalize**: _bool, default: `True`_

&emsp;If 'normalize=True', normalizing expression value to [0, 1]

**show_gene**: _str_

&emsp;The gene selected to plot

**colormap**: str, default: `viridis`_

&emsp;The name of cmap

**size**: _float, default: `10`_

&emsp;The size of point

**alpha**: _float, default: `1`_

&emsp;The transparency of point

****




