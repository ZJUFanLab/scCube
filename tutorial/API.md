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

DataFrame of input data

`sc_meta`: DataFrame of input meta

`is_normalized`: Whether the input data is normalized or not. If `is_normalized=False`, the input data will be normalized by scCube first. (default: `False`)










