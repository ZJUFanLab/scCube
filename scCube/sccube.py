from .utils import train_vae, generate_vae, load_vae, SPatternGenerator, generate_spot_data, generate_image_data, \
    get_min_grid_index, select_min_grid
import numpy as np
import pandas as pd
import torch
import ot
import scanpy as sc
import anndata as ad
from collections import Counter
from typing import Optional
import os
import os.path as osp
from pandas import DataFrame
from anndata import AnnData
import warnings
warnings.filterwarnings("ignore")


class scCube:
    def __init__(self):
        pass

    def pre_process(self,
                    sc_data: DataFrame,
                    sc_meta: DataFrame,
                    is_normalized: bool = True,
                    ):

        """
        pre-process input data, make sure it's normalized
        :param sc_data: DataFrame of input data
        :param sc_meta: DataFrame of input meta
        :param is_normalized: whether normalized
        :return: AnnData of data
        """
        assert sc_meta.shape[0] == sc_data.shape[1], "shape data no match!!!"

        sc_adata = ad.AnnData(sc_data.T)
        sc_adata.obs = sc_meta

        if not is_normalized:
            print("the input is count matrix, normalizing it firstly...")
            sc.pp.normalize_total(sc_adata, target_sum=1e4)
            sc.pp.log1p(sc_adata)

        return sc_adata

    def train_vae_and_generate_cell(self,
                                    sc_adata: AnnData,
                                    celltype_key: str,
                                    cell_key: str,
                                    target_num: Optional[dict] = None,
                                    batch_size: int = 512,
                                    epoch_num: int = 3500,
                                    lr: float = 0.0001,
                                    hidden_size: int = 128,
                                    save_model: bool = True,
                                    save_path: str = '',
                                    project_name: str = '',
                                    used_device: str = 'cuda:0'
                                    ):
        """
        train VAE model to generate cells
        :param sc_adata: AnnData of data
        :param celltype_key: the column name of `cell types` in meta
        :param cell_key: the column name of `cell` in meta
        :param target_num: target number of cells to generate, if `target_num=None`, generate cells by the proportion
        of cell types of the input data
        :param batch_size: batch size of training
        :param epoch_num: epoch number of training
        :param lr: learning reta of training
        :param hidden_size: hidden size of VAE model
        :param save_model: save trained VAE model or not
        :param save_path: save path
        :param project_name: name of trained VAE model
        :param used_device: device name, `cpu` or `cuda`
        :return: DataFrame of generated cell meta and data
        """
        res_list = self.__parpare_generate(sc_adata, celltype_key, cell_key, target_num)
        print('begin vae training...')
        # ********************* training *********************
        net = train_vae(res_list,
                        batch_size=batch_size,
                        epoch_num=epoch_num,
                        lr=lr,
                        hidden_size=hidden_size,
                        used_device=used_device)
        print('vae training done!')

        if save_model:
            print("saving the trained vae model...")
            path_save = os.path.join(save_path, f"{project_name}.pth")
            if not osp.exists(save_path):
                os.makedirs(save_path)
            torch.save(net.state_dict(), path_save)
            print(f"save trained vae in {path_save}.")

        # ********************* generating *********************
        generate_sc_meta, generate_sc_data = generate_vae(net,
                                                          res_list,
                                                          used_device=used_device)

        return generate_sc_meta, generate_sc_data

    def load_vae_and_generate_cell(self,
                                   sc_adata: AnnData,
                                   celltype_key: str,
                                   cell_key: str,
                                   target_num: Optional[dict] = None,
                                   hidden_size: int = 128,
                                   load_path: str = '',
                                   used_device: str = 'cuda:0'
                                   ):
        """
        load trained VAE model to generate cells
        :param sc_adata: AnnData of data
        :param celltype_key: the column name of `cell types` in meta
        :param cell_key: the column name of `cell` in meta
        :param target_num: target number of cells to generate, if `target_num=None`, generate cells by the proportion
        of cell types of the input data
        :param hidden_size: hidden size of VAE model, make sure it's same as the trained VAE model
        :param load_path: load path
        :param used_device: device name, `cpu` or `cuda`
        :return: DataFrame of generated cell meta and data
        """

        res_list = self.__parpare_generate(sc_adata, celltype_key, cell_key, target_num)
        print(f'loading model from {load_path}')
        # ********************* loading *********************
        net = load_vae(res_list,
                       hidden_size=hidden_size,
                       load_path=load_path,
                       used_device=used_device)
        print('vae loading done!')

        # ********************* generating *********************
        generate_sc_meta, generate_sc_data = generate_vae(net,
                                                          res_list,
                                                          used_device=used_device)

        return generate_sc_meta, generate_sc_data

    def generate_spatial_data_random(self,
                                     generate_sc_data: DataFrame,
                                     generate_sc_meta: DataFrame,
                                     set_seed: bool = False,
                                     seed: int = 12345,
                                     spatial_cell_type: Optional[list] = None,
                                     spatial_dim: int = 2,
                                     spatial_size: int = 30,
                                     delta: float = 25,
                                     lamda: float = 0.75,
                                     is_spot: bool = False,
                                     platform: str = 'ST',
                                     gene_type: str = 'whole',
                                     min_cell: int = 10,
                                     n_gene: Optional[int] = None,
                                     n_cell: int = 10,
                                     is_split: bool = True,
                                     split_coord: str = 'point_z',
                                     slice_num: int = 5,
                                     ):
        """
        generate random spatial patterns of cells with default spatial autocorrelation function
        :param generate_sc_data: DataFrame of generated sc data
        :param generate_sc_meta: DataFrame of generated sc meta
        :param set_seed: whether to set seed for reproducible simulation
        :param seed: seed number
        :param spatial_cell_type: the selected cell types with spatial patterns, if`spatial_cell_type=None`,
        all cell types would be assigned spatial patterns
        :param spatial_dim: spatial dimensionality， `2` or `3`
        :param spatial_size: the scope for spatial autocorrelation function, larger value will yield finer spatial
        patterns, but also will take more running time
        :param delta: larger value will yield more continued spatial patterns
        :param lamda: 0-1, larger value will yield more clear spatial patterns
        :param is_spot: True -- generate spot-based ST data; False -- generate image-based ST data
        :param platform: only works when `is_spot=True`, `ST` -- square neighborhood structure;
        `Visium` -- hexagonal neighborhood structure; `Slide` -- random neighborhood structure
        :param gene_type: the type of genes to generate, `whole` -- the whole genes;
        `hvg` -- the highly variable genes; `marker` -- the marker genes of each cell type;
        `random` -- the randomly selected genes
        :param min_cell: filter the genes expressed in fewer than `min_cell` cells before selected genes,
        only works when `gene_type='random', 'hvg', or 'marker'`
        :param n_gene: the number of genes to select, only works when `gene_type='random', 'hvg', or 'marker'`
        :param n_cell: the mean number of cells per spot, only works when `is_spot=True`
        :param is_split: whether to spilt the 3-D generated spatial patterns into a series of 2-D spatial patterns,
        only works when `spatial_dim=3`
        :param split_coord: the split coordinate axis, only works when `spatial_dim=3` and `is_split=True`
        :param slice_num: the targeted number of 2-D slices, only works when `spatial_dim=3` and `is_split=True`
        :return: generate_sc_data, generate_sc_meta_new, st_data, st_meta, st_index
        """
        spattern_generator = SPatternGenerator()
        generate_sc_meta_new = spattern_generator.run(generate_meta=generate_sc_meta,
                                                      set_seed=set_seed,
                                                      seed=seed,
                                                      spatial_cell_type=spatial_cell_type,
                                                      spatial_dim=spatial_dim,
                                                      spatial_size=spatial_size,
                                                      delta=delta,
                                                      lamda=lamda)

        generate_sc_meta_new = generate_sc_meta_new.loc[generate_sc_data.columns, ]

        if is_spot:
            st_data, st_meta, st_index = generate_spot_data(generate_meta=generate_sc_meta_new,
                                                            generate_data=generate_sc_data,
                                                            spatial_dim=spatial_dim,
                                                            platform=platform,
                                                            gene_type=gene_type,
                                                            min_cell=min_cell,
                                                            n_gene=n_gene,
                                                            n_cell=n_cell
                                                            )
        else:
            st_data, st_meta, st_index = generate_image_data(generate_meta=generate_sc_meta_new,
                                                             generate_data=generate_sc_data,
                                                             spatial_dim=spatial_dim,
                                                             gene_type=gene_type,
                                                             min_cell=min_cell,
                                                             n_gene=n_gene
                                                             )
            if spatial_dim == 3 and is_split:
                st_meta = self.__split_3d_coord(st_meta,
                                                split_coord=split_coord,
                                                n_slice=slice_num)

        return generate_sc_data, generate_sc_meta_new, st_data, st_meta, st_index


    def generate_pattern_reference(self,
                                   sc_adata: AnnData,
                                   generate_sc_data: DataFrame,
                                   generate_sc_meta: DataFrame,
                                   celltype_key: list,
                                   spatial_key: list,
                                   cost_metric: str = 'sqeuclidean'):
        """
        generate real spatial patterns of cells based on ST data
        :param sc_adata: AnnData of data
        :param generate_sc_data: DataFrame of generated sc data
        :param generate_sc_meta: DataFrame of generated sc meta
        :param celltype_key: the column name of `cell type` in meta
        :param spatial_key: the column name of `spatial coordinates` in meta
        :param cost_metric: the cost distance between generate_sc_data and real_data, `sqeuclidean` by default.
        On numpy the function also accepts from the scipy.spatial.distance.cdist function : ‘braycurtis’, ‘canberra’, ‘
        chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’,
        ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’,
        ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’, ‘yule’.
        :return: generate_sc_data, generate_sc_meta
        """

        # real_data = sc_adata.X
        # generate_data = np.array(generate_sc_data.T)
        real_meta = sc_adata.obs

        coord = pd.DataFrame()
        for i in list(set(real_meta[celltype_key])):
            sc_adata_tmp = sc_adata[sc_adata.obs[celltype_key].isin([i])]
            real_data_tmp = np.array(sc_adata_tmp.X)
            real_meta_tmp = sc_adata_tmp.obs

            generate_meta_tmp = generate_sc_meta.loc[generate_sc_meta['Cell_type'] == i]
            generate_meta_tmp.set_index('Cell')
            generate_data_tmp = generate_sc_data.iloc[:, generate_meta_tmp.index]
            generate_data_tmp = np.array(generate_data_tmp.T)

            m = real_data_tmp.shape[0]
            n = generate_data_tmp.shape[0]
            a, b = np.ones((m,)) / m, np.ones((n,)) / n  # uniform distribution on samples

            # loss matrix
            M = ot.dist(real_data_tmp, generate_data_tmp, metric=cost_metric)

            G = ot.emd(a, b, M)

            x_new = []
            y_new = []
            for i in range(n):
                x_tmp = real_meta_tmp.loc[real_meta_tmp.index[G[:, i].argmax()]][spatial_key[0]]
                y_tmp = real_meta_tmp.loc[real_meta_tmp.index[G[:, i].argmax()]][spatial_key[1]]
                x_new.append(x_tmp)
                y_new.append(y_tmp)

            coord_new = pd.DataFrame({spatial_key[0]: x_new, spatial_key[1]: y_new}, index=generate_meta_tmp.index)
            coord = coord.append(coord_new)

        coord = coord.loc[generate_sc_meta.index, ]
        generate_sc_meta = pd.concat([generate_sc_meta, coord], axis=1, join='outer')

        return generate_sc_data, generate_sc_meta

    def __parpare_generate(self,
                           sc_adata: AnnData,
                           celltype_key: str = 'Cell_type',
                           cell_key: str = 'Cell',
                           target_num: Optional[dict] = None,
                           ):

        single_cell = sc_adata.X

        index_2_gene = sc_adata.var_names.tolist()
        feature_size = single_cell.shape[1]

        breed = sc_adata.obs[celltype_key]
        breed_np = breed.values
        breed_set = set(breed_np)
        breed_2_list = list(breed_set)
        # sort by the first letter of the cell type to sure that the target_num is consistent with the breed_2_list
        breed_2_list = sorted(breed_2_list)
        dic = {}
        label = []

        for i in range(len(breed_set)):
            dic[breed_2_list[i]] = i

        cell = sc_adata.obs[cell_key].values
        for i in range(cell.shape[0]):
            label.append(dic[breed_np[i]])

        label = np.array(label)

        if target_num is None:
            print("generating by the proportion of cell types of the input scRNA-seq data...")
            cell_target_num = dict(Counter(breed_np))
        else:
            print("generating by the targeted proportion of cell types...")
            cell_target_num = target_num
            input_breed_2_list = sorted(list(cell_target_num.keys()))
            assert input_breed_2_list == breed_2_list, \
                "the cell type in targeted number list isn't same as it in scRNA-seq meta file!!!"

        # label index the data size of corresponding target
        cell_number_target_num = {}
        for k, v in cell_target_num.items():
            cell_number_target_num[dic[k]] = v

        res_list = []
        res_list.extend((feature_size, single_cell, label, breed_2_list, index_2_gene, cell_number_target_num))

        return res_list

    def __split_3d_coord(self,
                         obj: DataFrame,
                         split_coord: str = 'point_x',
                         n_slice: int = 10,
                         ):
        """
        Split 3d coordinates to n 2d-slices
        :param obj: DataFrame of meta
        :param split_coord: The coordinate to split
        :param n_slice: The number of slices to split
        :return: New DataFrame with slice information
        """
        min_range = int(np.floor(obj[split_coord].min()))
        max_range = int(np.ceil(obj[split_coord].max()))
        grid_seq = [i / n_slice for i in range(min_range, int(max_range * (n_slice + 1)), max_range)]

        obj['slice'] = "unassigned"
        for i in range(0, obj.shape[0]):
            x = obj.loc[obj.index[i], split_coord]
            obj.iloc[i, 5] = get_min_grid_index(select_min_grid(x, grid_seq), grid_seq)

        return obj




