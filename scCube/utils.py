import torch
import pandas as pd
import numpy as np
import scanpy as sc
import copy
import math
import random
from tqdm import tqdm
from torch import nn
from torch.optim import AdamW
from .model import VAE
from torch.utils.data import Dataset, DataLoader
from collections import Counter
from pandas import DataFrame
from typing import Optional
import pdb
import warnings

warnings.filterwarnings('ignore')


# scCube generate single cell
class myDataset(Dataset):  # 需要继承data.Dataset
    def __init__(self, single_cell, label):
        # 1. Initialize file path or list of file names.
        self.sc = single_cell
        self.label = label

    def __getitem__(self, idx):
        # TODO
        tmp_x = self.sc[idx]
        tmp_y_tag = self.label[idx]

        return (tmp_x, tmp_y_tag)  # tag 分类

    def __len__(self):
        # You should change 0 to the total size of your dataset.
        return self.label.shape[0]


def train_vae(
        res_list: list,
        batch_size: int = 512,
        epoch_num: int = 3500,
        lr: float = 0.0001,
        hidden_size: int = 128,
        used_device: str = 'cuda:0'
):
    """
    Training VAE model
    :param res_list: the result list prepared by `prepare_generate` function
    :param batch_size: batch size, eg:512
    :param epoch_num: num epochs, eg:3500
    :param lr: learning rate, eg:0.0001
    :param hidden_size: hidden size, eg:128
    # :param early_stopping: The model waits N epochs before stops since the training loss does not decline, eg:100
    # :param not_early_stop: Whether to use the `early_stop` strategy
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: trained VAE model
    """
    # random_seed = args.random_seed

    hidden_list = [2048, 1024, 512]
    mid_hidden_size = hidden_size
    weight_decay = 5e-4
    dataloader = DataLoader(myDataset(single_cell=res_list[1], label=res_list[2]), batch_size=batch_size,
                            shuffle=True, pin_memory=True)
    criterion = nn.MSELoss()

    vae = VAE(res_list[0], hidden_list, mid_hidden_size).to(used_device)
    optimizer = AdamW(vae.parameters(), lr=lr, weight_decay=weight_decay)

    pbar = tqdm(range(epoch_num))
    min_loss = 1000000000000000
    vae.train()
    early_stop = 0

    for epoch in pbar:
        train_loss = 0

        for batch_idx, data in enumerate(dataloader):
            cell_feature, label = data
            cell_feature = torch.tensor(cell_feature, dtype=torch.float32).to(used_device)

            x_hat, kl_div = vae(cell_feature, used_device)
            loss = criterion(x_hat, cell_feature)

            if kl_div is not None:
                loss += kl_div

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        if train_loss < min_loss:
            min_loss = train_loss
            best_vae = copy.deepcopy(vae)
            # early_stop = 0
        else:
            early_stop += 1
        pbar.set_description('Train Epoch: {}'.format(epoch))
        pbar.set_postfix(loss=f"{train_loss:.4f}", min_loss=f"{min_loss:.4f}")
        # if early_stop > early_stopping and not not_early_stop:  # use early_stop
        #     break

    return best_vae


def load_vae(
        res_list: list,
        load_path: str = '',
        hidden_size: int = 128,
        used_device: str = 'cuda:0'
):
    """
    Loading trained VAE model
    :param res_list: the result list prepared by `prepare_generate` function
    :param load_path: the load directory
    :param hidden_size: hidden size, eg:128
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: trained VAE model
    """

    hidden_list = [2048, 1024, 512]
    vae = VAE(res_list[0], hidden_list, hidden_size).to(used_device)
    vae.load_state_dict(torch.load(load_path, map_location=used_device))
    return vae


def prepare_data(
        cell_all_generate,
        label_all_generate,
        breed_2_list,
        index_2_gene
):
    """
    Convert array to DataFrame
    :param cell_all_generate: cells to generate
    :param label_all_generate: cell type labels to generate
    :param breed_2_list: the cell type list
    :param index_2_gene: the features/genes list
    :return: DataFrame of generated single cell expression profiles and metadata
    """

    cell_all_generate = np.array(cell_all_generate)
    label_all_generate = np.array(label_all_generate)

    cell_all_generate_csv = pd.DataFrame(cell_all_generate)
    label_all_generate_csv = pd.DataFrame(label_all_generate)

    ids = label_all_generate_csv[0].tolist()
    breeds = []
    for id in ids:
        breeds.append(breed_2_list[id])
    name = ["C_" + str(i + 1) for i in range(label_all_generate.shape[0])]

    label_all_generate_csv.insert(1, "Cell_type", np.array(breeds))
    label_all_generate_csv.insert(1, "Cell", np.array(name))
    label_all_generate_csv = label_all_generate_csv.drop([0], axis=1)

    cell_all_generate_csv = cell_all_generate_csv.T
    cell_all_generate_csv.columns = name
    cell_all_generate_csv.index = index_2_gene

    return label_all_generate_csv, cell_all_generate_csv


def generate_vae(
        net,
        res_list: list,
        used_device: str = 'cuda:0'
):
    """
    Generate single cell expression profiles
    :param net: trained VAE model
    :param res_list: the result list prepared by `prepare_generate` function
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: generated single cell expression profiles and metadata
    """
    cell_number_target_num = res_list[5]

    # net in cuda now
    for p in net.parameters():  # reset requires_grad
        p.requires_grad = False  # avoid computation

    net.eval()
    net.to(used_device)
    cell_all_generate = []
    label_all_generate = []

    all_to_generate = 0
    for x in cell_number_target_num.values():
        all_to_generate += x

    if cell_number_target_num != None:
        epochs = 10000  # 10000次
        ratio = 1
    else:
        epochs = 1
        ratio = -1

    # fmt = '{}'.format
    cell_feature = torch.from_numpy(res_list[1]).float()
    label = torch.from_numpy(res_list[2])
    ##############

    with torch.no_grad():
        with tqdm(total=all_to_generate) as pbar:
            for epoch in range(epochs):
                key_list = []  # list
                generate_num = 0

                label_list = label.tolist()
                for i in range(len(label_list)):
                    if cell_number_target_num[label_list[i]] <= 0:
                        continue
                    else:
                        cell_number_target_num[label_list[i]] -= 1
                        generate_num += 1
                        key_list.append(i)

                if cell_number_target_num == None or all_to_generate == 0 or len(key_list) == 0:
                    assert all_to_generate == 0 and len(key_list) == 0
                    break

                # 随机打乱
                random.shuffle(key_list)

                label = label.index_select(0, torch.tensor(key_list))
                cell_feature = cell_feature.index_select(0, torch.tensor(key_list))

                dataloader = DataLoader(myDataset(single_cell=cell_feature, label=label), batch_size=300, shuffle=False,
                                        pin_memory=True, num_workers=0)
                for batch_idx, data in enumerate(dataloader):  # 一个batch
                    cell_feature_batch, label_batch = data
                    cell_feature_batch = cell_feature_batch.to(used_device)
                    pbar.set_description('Generate Epoch: {}'.format(epoch))

                    label_batch = label_batch.cpu().numpy()

                    for j in range(ratio):  # 翻倍多少
                        ans_l, _ = net(cell_feature_batch, used_device)
                        ans_l = ans_l.cpu().data.numpy()
                        # for i in range(ans_l.shape[0]):
                        cell_all_generate.extend(ans_l)
                        label_all_generate.extend(label_batch)

                all_to_generate -= generate_num
                pbar.update(generate_num)

    print("generated done!")
    generate_sc_meta, generate_sc_data = prepare_data(cell_all_generate, label_all_generate, res_list[3], res_list[4])
    print("data have been prepared!")
    return generate_sc_meta, generate_sc_data


# generate spatial pattern (single cell)
class SPatternGenerator:
    def __init__(self):
        pass

    def run(self,
            generate_meta: DataFrame,
            set_seed: bool = False,
            seed: int = 12345,
            spatial_cell_type: Optional[list] = None,
            spatial_dim: int = 2,
            spatial_size: int = 30,
            delta: float = 25,
            lamda: float = 0.75
            ):

        """
            Generate spatial pattern of cell types with single cell resolution
            :param generate_meta: generate single cell meta from VAE
            :param set_seed: False--random True--set seed
            :param seed: random seed
            :param spatial_cell_type: select cell types with spatial pattern, if None, all cell types have spatial patterns
            :param spatial_dim: spatial pattern dim
            :param spatial_size: spatial pattern size
            :param delta: spatial pattern delta, large--big spatial pattern
            :param lamda: spatial pattern lamda (0-1), large--clear pattern small--fuzzy pattern
            :return: DataFrame of single cell meta with spatial coordinates
            """

        if set_seed:
            np.random.seed(seed)

        print('generating spatial coordinates of single cells...')

        # Generates a realisation of the Poisson point process
        sim_point = self.__create_sim_pois(spatial_dim=spatial_dim, cell_num=generate_meta.shape[0]) * spatial_size

        if int(spatial_dim) == 2:
            sim_point['grid_x'] = np.floor(np.array(sim_point['point_x'])) + 0.5
            sim_point['grid_y'] = np.floor(np.array(sim_point['point_y'])) + 0.5
            # grid
            grid = ([(x, y) for x in range(spatial_size) for y in range(spatial_size)])

        elif int(spatial_dim) == 3:
            sim_point['grid_x'] = np.floor(np.array(sim_point['point_x'])) + 0.5
            sim_point['grid_y'] = np.floor(np.array(sim_point['point_y'])) + 0.5
            sim_point['grid_z'] = np.floor(np.array(sim_point['point_z'])) + 0.5
            # grid
            grid = ([(x, y, z) for x in range(spatial_size) for y in range(spatial_size) for z in range(spatial_size)])

        # create distance matrix
        distance = pd.DataFrame(np.array(self.__cal_dist(coord=grid)))
        # Generate random variable
        mean = np.repeat(0, distance.shape[0])
        cov = np.exp(-(distance.values / delta) ** 2)
        sample = np.random.multivariate_normal(mean, cov, tol=1e-6)
        sample_dic = {}
        for i in range(len(sample)):
            sample_dic[i] = sample[i]

        grid = pd.DataFrame(grid)
        if int(spatial_dim) == 2:
            grid.columns = ['grid_x', 'grid_y']
        elif int(spatial_dim) == 3:
            grid.columns = ['grid_x', 'grid_y', 'grid_z']

        grid_out = grid + 0.5
        grid_out['Cell_type'] = "unassigned"

        # get proportion of each cell type from generate_meta
        cell_type_num_dict = dict(Counter(generate_meta['Cell_type'].values))
        cell_num = generate_meta.shape[0]

        for k, v in cell_type_num_dict.items():
            prop = v / cell_num
            dic1 = {key: value for key, value in sample_dic.items() if value < np.quantile(sample, prop)}
            grid_out.loc[dic1.keys(), ['Cell_type']] = k
            sample_dic = {key: value for key, value in sample_dic.items() if value >= np.quantile(sample, prop)}
            sample = np.array(np.array(list(sample_dic.values())))
            cell_num = cell_num - v

        # if args.spatial_cell_type is None:
        #     print('generating spatial patterns of totally ' + str(len(cell_type_num_dict)) + ' cell types...')
        #     for k, v in cell_type_num_dict.items():
        #         prop = v / cell_num
        #         dic1 = {key: value for key, value in sample_dic.items() if value < np.quantile(sample, prop)}
        #         grid_out.loc[dic1.keys(), ['Cell_type']] = k
        #         sample_dic = {key: value for key, value in sample_dic.items() if value >= np.quantile(sample, prop)}
        #         sample = np.array(np.array(list(sample_dic.values())))
        #         cell_num = cell_num - v
        # else:
        #     choose_cell_type_list = args.spatial_cell_type
        #     print('generating spatial patterns of selected ' + str(len(choose_cell_type_list)) + ' cell types...')
        #     choose_cell_type_num_dict = dict(Counter(generate_meta.loc[generate_meta['Cell_type'].isin(choose_cell_type_list), 'Cell_type'].values))
        #     for k, v in choose_cell_type_num_dict.items():
        #         prop = v / cell_num
        #         dic1 = {key: value for key, value in sample_dic.items() if value < np.quantile(sample, prop)}
        #         grid_out.loc[dic1.keys(), ['Cell_type']] = k
        #         sample_dic = {key: value for key, value in sample_dic.items() if value >= np.quantile(sample, prop)}
        #         sample = np.array(np.array(list(sample_dic.values())))
        #         cell_num = cell_num - v

        # other cell type without spatial patterns
        # pdb.set_trace()
        # unchoose_cell_type_num_dict = dict(
        #     Counter(generate_meta.loc[~generate_meta['Cell_type'].isin(choose_cell_type_list), 'Cell_type'].values))
        # for k, v in unchoose_cell_type_num_dict.items():
        #     idx = grid_out.loc[grid_out['Cell_type'] == 'unassigned', ].index
        #     grid_num = int(np.round(v / generate_meta.shape[0] * grid_out.shape[0]))
        #     if len(idx) >= grid_num:
        #         idx_choose = np.random.choice(idx, grid_num, replace=False).tolist()
        #     else:
        #         idx_choose = np.random.choice(idx, grid_num, replace=True).tolist()
        #     grid_out.loc[idx_choose, 'Cell_type'] = k

        # handle unassigned
        cell_type_name = np.unique(generate_meta['Cell_type'].values)
        replace_num = len(grid_out.loc[grid_out['Cell_type'] == 'unassigned', ['Cell_type']])
        grid_out.loc[grid_out['Cell_type'] == 'unassigned', ['Cell_type']] = np.random.choice(cell_type_name,
                                                                                              replace_num)

        if spatial_cell_type is None:
            print('generating spatial patterns of totally ' + str(len(cell_type_num_dict)) + ' cell types...')
            spa_pattern = pd.merge(grid_out, sim_point, how='inner')
        else:
            choose_cell_type_list = spatial_cell_type
            print('generating spatial patterns of selected ' + str(len(choose_cell_type_list)) + ' cell types...')
            spa_pattern = pd.merge(grid_out, sim_point, how='inner')
            # random shuffle other cell type coordinates
            idx = spa_pattern.loc[~spa_pattern['Cell_type'].isin(choose_cell_type_list),].index
            idx_cell_type = spa_pattern.loc[
                ~spa_pattern['Cell_type'].isin(choose_cell_type_list), 'Cell_type'].values.tolist()
            random.shuffle(idx_cell_type)
            spa_pattern.loc[idx, ['Cell_type']] = idx_cell_type

        # correct cell type prop
        generate_cell_type_num_dict = dict(Counter(spa_pattern['Cell_type'].values))
        add_dict = {}
        sub_dict = {}
        for i in list(cell_type_name):
            if generate_cell_type_num_dict[i] < cell_type_num_dict[i]:
                add_dict[i] = cell_type_num_dict[i] - generate_cell_type_num_dict[i]
            elif generate_cell_type_num_dict[i] > cell_type_num_dict[i]:
                sub_dict[i] = generate_cell_type_num_dict[i] - cell_type_num_dict[i]

        # pdb.set_trace()

        correct_index = []
        for k, v in sub_dict.items():
            correct_index.extend(
                np.random.choice(spa_pattern[spa_pattern['Cell_type'] == k].index.values, v, replace=False).tolist())

        add_ct_list = []
        for k, v in add_dict.items():
            add_ct_list.extend([k] * v)

        random.shuffle(add_ct_list)
        spa_pattern.loc[correct_index, ['Cell_type']] = add_ct_list

        # pdb.set_trace()

        # lamda (0-1): large--clear pattern small--fuzzy pattern
        lamda_len = round(len(spa_pattern) * (1 - lamda))
        select_idx = np.random.choice(spa_pattern.index.values, lamda_len, replace=False)
        change_idx = np.random.choice(select_idx, len(select_idx), replace=False)
        spa_pattern.loc[select_idx, ['Cell_type']] = np.array(spa_pattern.loc[change_idx, 'Cell_type'])

        # match Cell - Cell_type - coord
        generate_meta.sort_values(by=['Cell_type'], inplace=True)
        spa_pattern.sort_values(by=['Cell_type'], inplace=True)
        generate_meta.set_index(generate_meta['Cell'], inplace=True)
        spa_pattern.set_index(generate_meta['Cell'], inplace=True)
        if int(spatial_dim) == 2:
            generate_meta = pd.concat([generate_meta, spa_pattern.iloc[:, 3:]], axis=1)
        elif int(spatial_dim) == 3:
            generate_meta = pd.concat([generate_meta, spa_pattern.iloc[:, 4:]], axis=1)

        return generate_meta

    def __cal_dist(self, coord):
        """
            Calculate Euclidean distance between points
            """
        dis = []
        for i in range(len(coord)):
            dis.append([])
            for j in range(len(coord)):
                dis[i].append(math.dist(coord[i], coord[j]))
        return dis

    def __create_sim_pois(self,
                        spatial_dim: int,
                        cell_num: int):
        """
        Generates a realisation of the Poisson point process
        :param spatial_dim: spatial pattern dimensionality, 2 or 3
        :param cell_num: the number of points
        :return: DataFrame of spatial coordinates of points
        """
        assert spatial_dim in [2, 3], "error define spatial pattern dim!!"

        # if args.set_seed:
        #     np.random.seed(args.random_seed)

        if int(spatial_dim) == 2:
            xx = np.random.uniform(0, 1, cell_num)
            yy = np.random.uniform(0, 1, cell_num)
            out = pd.DataFrame({'point_x': xx, 'point_y': yy})
        elif int(spatial_dim) == 3:
            xx = np.random.uniform(0, 1, cell_num)
            yy = np.random.uniform(0, 1, cell_num)
            zz = np.random.uniform(0, 1, cell_num)
            out = pd.DataFrame({'point_x': xx, 'point_y': yy, 'point_z': zz})
        return out


def select_gene(
        generate_meta: DataFrame,
        generate_data: DataFrame,
        min_cell: int = 10,
        gene_type: str = 'random',
        n_gene: int = 200
):
    """
    Select targeted genes for imaged-based spatial transcriptomics data
    :param generate_meta: generate single cell meta from VAE
    :param generate_data: generate single cell data from VAE
    :param min_cell: filter the genes expressed in fewer than `min_cell` cells
    :param gene_type: targeted genes type, whole--whole genes; hvg--highly variable genes;
    marker--marker genes; random--random chosen genes
    :param n_gene: targeted genes number
    :return: selected targeted genes
    """

    # filter the genes expressed in fewer than `min_cell` cells
    filter_gene = (generate_data > 0).sum(axis=1) > min_cell
    generate_data = generate_data[filter_gene]

    if gene_type == 'random':
        gene_list = list(generate_data.index)
        gene = random.sample(gene_list, n_gene)
    else:
        adata = sc.AnnData(generate_data.astype(np.float32).T)
        adata.obs_names = generate_data.columns
        adata.var_names = generate_data.index
        adata.obs['Cell_type'] = generate_meta['Cell_type'].values
        if gene_type == 'hvg':
            sc.pp.highly_variable_genes(adata, n_top_genes=n_gene, flavor='seurat_v3')
            adata = adata[:, adata.var.highly_variable]
            gene = list(adata.var.index)
        elif gene_type == 'marker':
            sc.tl.rank_genes_groups(adata, 'Cell_type', method='wilcoxon')
            marker_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(n_gene)
            gene_list = list(np.unique(np.ravel(np.array(marker_df))))
            gene = random.sample(gene_list, n_gene)

    return gene


def generate_spot_data(
        generate_meta: DataFrame,
        generate_data: DataFrame,
        spatial_dim: int = 2,
        platform: str = 'ST',
        gene_type: str = 'whole',
        min_cell: int = 10,
        n_gene: Optional[int] = None,
        n_cell: int = 10,
):
    """
    Generate spot-based spatial transcriptomics data
    :param generate_meta: generate single cell meta from VAE
    :param generate_data: generate single cell data from VAE
    :param spatial_dim: spatial pattern dimensionality
    :param platform: spatial sequencing platform
    :param gene_type: targeted genes type
    :param min_cell: filter the genes expressed in fewer than `min_cell` cells
    :param n_gene: targeted genes number
    :param n_cell: cell number per spot
    :return:
        spot_data: spot-based spatial transcriptomics data
        spot_meta: spot-based spatial transcriptomics meta
        generate_meta: DataFrame of cell-spot index
    """

    # check spatial dim
    assert spatial_dim == 2, "the spatial pattern dim is supposed to 2 for spot-based ST data!!"
    # check spatial sequencing platform
    assert platform in ['ST', 'Visium', 'Slide'], "the spatial sequencing platform is supposed to ST, Visium, or Slide!!"
    # check gene type
    assert gene_type in ['whole', 'hvg', 'marker', 'random'], "the gene type is supposed to whole, hvg, marker, or random!!"

    # process meta firstly
    print('generating spot-based ST data with ' + str(n_cell) + ' cells per spot...')
    gap = np.floor(np.sqrt(generate_meta.shape[0] / n_cell))
    grid_range = int(np.ceil(max(max(generate_meta['point_x'].values), max(generate_meta['point_y'].values))))
    grid_seq = [i / gap for i in range(0, int(grid_range * (gap + 1)), grid_range)]
    # pdb.set_trace()
    generate_meta['spot_x'] = "unassigned"
    generate_meta['spot_y'] = "unassigned"
    for i in range(0, generate_meta.shape[0]):
        # spot_x TODO:
        generate_meta.iloc[i, 4] = select_min_grid(generate_meta.iloc[i, 2], grid_seq)
        # spot_y
        generate_meta.iloc[i, 5] = select_min_grid(generate_meta.iloc[i, 3], grid_seq)

    # assign spot name
    generate_meta['spot'] = "unassigned"
    x_gap = list(np.unique(generate_meta['spot_x'].values))
    x_gap.sort()
    y_gap = list(np.unique(generate_meta['spot_y'].values))
    y_gap.sort()

    for i in range(0, len(x_gap)):
        for j in range(0, len(y_gap)):
            if generate_meta[(generate_meta['spot_x'] == x_gap[i]) &
                             (generate_meta['spot_y'] == y_gap[j])].shape[0] > 0:
                generate_meta.loc[(generate_meta['spot_x'] == x_gap[i]) &
                                  (generate_meta['spot_y'] == y_gap[j]), 'spot'] = ('spot_' + str(i * len(y_gap) + j + 1))

    # platform
    if platform == 'ST':
        print('generating with the spot layout and neighborhood structure of ' + str(platform) + ' (square)...')
    elif platform == 'Visium':
        print('generating with the spot layout and neighborhood structure of ' + str(platform) + ' (hex)...')
        x_gap_select = x_gap[1::2]
        grid_radius = grid_range / (gap * 2)
        generate_meta.loc[generate_meta['spot_x'].isin(x_gap_select), 'spot_y'] = \
            generate_meta.loc[generate_meta['spot_x'].isin(x_gap_select), 'spot_y'] + grid_radius
    elif platform == 'Slide':
        print('generating with the spot layout and neighborhood structure of ' + str(platform) + '-seq (random)...')
        grid_radius = grid_range / (gap * 2)
        spot_name = list(np.unique(generate_meta['spot'].values))
        # sort: [spot_1, spot_2, spot_3, ..., spot_n]
        spot_name = sorted(spot_name, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        length = np.random.uniform(0, grid_radius, len(spot_name))
        angle = np.pi * np.random.uniform(0, 2, len(spot_name))
        for i in range(0, len(spot_name)):
            generate_meta.loc[generate_meta['spot'] == spot_name[i], 'spot_x'] = \
                generate_meta.loc[generate_meta['spot'] == spot_name[i], 'spot_x'] + length[i] * np.cos(angle[i])
            generate_meta.loc[generate_meta['spot'] == spot_name[i], 'spot_y'] = \
                generate_meta.loc[generate_meta['spot'] == spot_name[i], 'spot_y'] + length[i] * np.sin(angle[i])

    # process data secondly
    spot_data = pd.DataFrame()
    spot_meta = pd.DataFrame(columns=['spot', 'spot_x', 'spot_y'])
    coord_list = list(set(list(zip(generate_meta.spot_x, generate_meta.spot_y, generate_meta.spot))))
    # sort: [spot_1, spot_2, spot_3, ..., spot_n]
    coord_list = sorted(coord_list, key=lambda tup: int("".join([i for i in tup[2] if i.isdigit()])))
    for i in range(0, len(coord_list)):
        row = pd.DataFrame([{'spot': coord_list[i][2], 'spot_x': coord_list[i][0], 'spot_y': coord_list[i][1]}])
        spot_meta = pd.concat([spot_meta, row])
        cell_pool = generate_meta.loc[(generate_meta['spot_x'] == coord_list[i][0]) &
                                      (generate_meta['spot_y'] == coord_list[i][1]), 'Cell']
        spot = generate_data[cell_pool].sum(axis=1)
        if spot.sum() > 25000:
            spot *= 20000 / spot.sum()
        # spot = generate_data[cell_pool].mean(axis=1)
        spot_id = coord_list[i][2]
        spot_data.insert(len(spot_data.columns), spot_id, spot)

    # gene type
    if gene_type == 'whole':
        print('generating with whole genes...')
    else:
        # check gene number
        assert n_gene is not None, "please input the exact number of genes!!"
        assert n_gene <= spot_data.shape[0], "the number of selected genes is too large, please reduce and try again!!"

        if gene_type == 'hvg':
            print('generating with ' + str(n_gene) + ' selected HVGs...')
        elif gene_type == 'marker':
            print('generating with ' + str(n_gene) + ' selected marker genes...')
        elif gene_type == 'random':
            print('generating with ' + str(n_gene) + ' randomly selected genes...')

        genes = select_gene(
            generate_meta=generate_meta,
            generate_data=generate_data,
            min_cell=min_cell,
            gene_type=gene_type,
            n_gene=n_gene
        )
        spot_data = spot_data.loc[genes]

    return spot_data, spot_meta, generate_meta


def generate_image_data(
        generate_meta: DataFrame,
        generate_data: DataFrame,
        spatial_dim: int,
        gene_type: str,
        min_cell: int,
        n_gene: int,
):
    """
    Generate image-based spatial transcriptomics data
    :param generate_meta: generate single cell meta from VAE
    :param generate_data: generate single cell data from VAE
    :param spatial_dim: spatial pattern dimensionality
    :param gene_type: targeted genes type
    :param min_cell: filter the genes expressed in fewer than `min_cell` cells
    :param n_gene: targeted genes number
    :return:
        target_data: image-based spatial transcriptomics data
        generate_meta: image-based spatial transcriptomics meta
        generate_data: ground truth of untargeted genes
    """
    print('generating image-based ST data with ' + str(n_gene) + ' targeted genes...')
    # check spatial dim
    assert spatial_dim in [2, 3], "the spatial pattern dim is supposed to 2 or 3 for image-based ST data!!"
    # check gene type
    assert gene_type in ['whole', 'hvg', 'marker',
                         'random'], "the gene type is supposed to whole, hvg, marker, or random!!"
    # check gene number
    assert n_gene <= generate_data.shape[0], "the number of selected genes is too large, please reduce and try again!!"

    # process data only
    if gene_type == 'hvg':
        print('generating image-based data with ' + str(n_gene) + ' targeted HVGs...')
    elif gene_type == 'marker':
        print('generating with ' + str(n_gene) + ' targeted marker genes...')
    elif gene_type == 'random':
        print('generating with ' + str(n_gene) + ' randomly targeted genes...')

    genes = select_gene(
            generate_meta=generate_meta,
            generate_data=generate_data,
            min_cell=min_cell,
            gene_type=gene_type,
            n_gene=n_gene
        )
    target_data = generate_data.loc[genes]

    return target_data, generate_meta, generate_data


def calculate_spot_prop(
        obj: DataFrame,
        cell_id: str = 'Cell',
        label: str = 'Cell_type',
        spot_id: str = 'spot'
):
    """
    Calculate cell type proportion of each spot
    :param obj: DataFrame of cell-spot index
    :param cell_id: The name of column containing cell id
    :param label: The name of column containing cell type information
    :param spot_id: The name of column containing spot id
    :return: DataFrame of cell type proportion per spot
    """

    prop = obj[[cell_id, label, spot_id]].pivot_table(index=[spot_id], columns=[label],
                                                      aggfunc='count', values=cell_id, fill_value=0)
    prop = prop.div(prop.sum(axis=1), axis=0)
    prop.columns = pd.Index(list(prop.columns))
    prop['spot_x'] = np.array(obj[['spot_x', 'spot_y', 'spot']].pivot_table(index=[spot_id])['spot_x'])
    prop['spot_y'] = np.array(obj[['spot_x', 'spot_y', 'spot']].pivot_table(index=[spot_id])['spot_y'])

    return prop


def select_min_grid(
        x,
        grid_list):
    """
    Select minimum grid as spot coordinates
    """
    choose_grid_list = []
    for i in grid_list:
        if x < i:
            choose_grid_list.append(i)

    choose_grid = np.min(choose_grid_list)

    return choose_grid


def get_min_grid_index(
        x,
        grid_list):
    """
    Get index of minimum grid selected in grid_list
    """
    idx = 0
    for i in range(len(grid_list)):
        if x == grid_list[i]:
            idx = i

    return idx

