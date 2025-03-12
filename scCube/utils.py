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
from .model import VAE, CLASSIFIER, LINEAR_LOGSOFTMAX, MLP_G, MLP_CRITIC
from torch.utils.data import Dataset, DataLoader
from torch.autograd import Variable
import torch.autograd as autograd
from collections import Counter
from pandas import DataFrame
from typing import Optional
import torch.optim as optim
from anndata import AnnData
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


class labelDataset(Dataset):  # 需要继承data.Dataset
    def __init__(self, label):
        # 1. Initialize file path or list of file names.
        self.label = label

    def __getitem__(self, idx):
        tmp_y_tag = self.label[idx]

        return tmp_y_tag  # tag 分类

    def __len__(self):
        # You should change 0 to the total size of your dataset.
        return self.label.shape[0]


def train_vae(
        res_list: list,
        batch_size: int = 512,
        epoch_num: int = 5000,
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


def calc_gradient_penalty(netD, real_data, fake_data, input_sem="", used_device: str = 'cuda:0'):
    # print real_data.size()
    alpha = torch.rand(real_data.shape[0], 1)
    alpha = alpha.expand(real_data.size())
    alpha = alpha.to(used_device)
    interpolates = alpha * real_data + ((1 - alpha) * fake_data)
    interpolates = interpolates.to(used_device)
    interpolates = Variable(interpolates, requires_grad=True)

    disc_interpolates = netD(interpolates, Variable(input_sem))

    ones = torch.ones(disc_interpolates.size())
    ones = ones.to(used_device)
    gradients = autograd.grad(outputs=disc_interpolates, inputs=interpolates,
                              grad_outputs=ones,
                              create_graph=True, retain_graph=True, only_inputs=True)[0]
    # args.GP_Weight = 10
    GP_Weight = 10
    gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean() * GP_Weight
    return gradient_penalty


def train_cgan(
        res_list: list,
        batch_size: int = 512,
        epoch_num: int = 5000,
        lr: float = 0.0001,
        SemSize: int = 256,
        NoiseSize: int = 256,
        used_device: str = 'cuda:0'):
    """
    Training CGAN model
    :param res_list: the result list prepared by `prepare_generate` function
    :param batch_size: batch size, eg:512
    :param epoch_num: num epochs, eg:5000
    :param lr: learning rate, eg:0.0001
    :param SemSize: dim of embedding of class, eg:256
    :param NoiseSize: dim of noise, eg:256
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: trained CGAN model
    """
    dataloader = DataLoader(myDataset(single_cell=res_list[1], label=res_list[2]), batch_size=batch_size, shuffle=True)
    netG = MLP_G(SemSize=SemSize, NoiseSize=NoiseSize, NGH=4096, FeaSize=res_list[0]).to(used_device)
    netD = MLP_CRITIC(SemSize=SemSize, NDH=4096, FeaSize=res_list[0]).to(used_device)
    cls_criterion = nn.NLLLoss().to(used_device)  # cross entropy loss
    optimizerD = optim.Adam(netD.parameters(), lr=lr, betas=(0.5, 0.999))
    optimizerG = optim.Adam(netG.parameters(), lr=lr, betas=(0.5, 0.999))
    pretrain_cls = CLASSIFIER(res_list[0], len(set(res_list[2])), res_list[6], SemSize, dataloader, used_device,
                              _lr=0.001, _nepoch=200)

    for p in pretrain_cls.model.parameters():  # set requires_grad to False
        p.requires_grad = False

    pbar = tqdm(range(epoch_num))
    # begain training
    input_fea = torch.FloatTensor(batch_size, res_list[0]).to(used_device)
    input_sem = torch.FloatTensor(batch_size, SemSize).to(used_device)
    noise = torch.FloatTensor(batch_size, NoiseSize).to(used_device)
    # one = torch.FloatTensor([1])
    #
    one = torch.tensor(1, dtype=torch.float)
    mone = one * -1
    one = one.to(used_device)
    mone = mone.to(used_device)
    input_label = torch.LongTensor(batch_size).to(used_device)

    min_loss_d = 10000000000000
    min_loss_g = 10000000000000
    ntrain = 600
    Critic_Iter = 5
    Cls_Weight = 0.01

    for epoch in pbar:
        train_loss_d = []
        train_loss_g = []
        WD = []
        for i in range(0, ntrain, batch_size):
            iter_d = 0
            for p in netD.parameters():
                p.requires_grad = True
            for batch_idx, data in enumerate(dataloader):

                cell_feature, label = data  # dim of cell sample cell ,label should map to embedding, which dim=dim_noise
                noise = torch.FloatTensor(len(label), NoiseSize).to(used_device)  # (64, 500)

                cell_feature = torch.tensor(cell_feature, dtype=torch.float32).to(used_device)
                netD.zero_grad()
                input_feav = Variable(cell_feature)
                input_sem = pretrain_cls.classifier_embed(
                    label).to(used_device)  # get the embedding of a class. input:4 out put 1024
                input_semv = Variable(input_sem)

                if iter_d < Critic_Iter:  # 每次的样本都不一样
                    iter_d += 1
                    # loss of real data

                    criticD_real = netD(input_feav, input_semv)
                    criticD_real = criticD_real.mean()
                    criticD_real.backward(mone)  # D_real is what we want to maximize, so minimize the loss (-D_real)
                    # loss of generated data
                    noise.normal_(0, 1)
                    noisev = Variable(noise)
                    fake = netG(noisev, input_semv)  # generate samples
                    # detach(): return a new variable, do not compute gradient for it
                    criticD_fake = netD(fake.detach(), input_semv)
                    criticD_fake = criticD_fake.mean()
                    criticD_fake.backward(one)

                    # loss with Lipschitz constraint
                    gradient_penalty = calc_gradient_penalty(netD, cell_feature, fake.data, input_sem,
                                                             used_device=used_device)
                    gradient_penalty.backward()

                    Wasserstein_D = criticD_real - criticD_fake
                    # Final Loss of Discriminator
                    D_cost = criticD_fake - criticD_real + gradient_penalty  # 判别器希望criticD_real大，criticD_fake，gradient_penalty小
                    optimizerD.step()
                    train_loss_d.append(D_cost.item())
                    WD.append(Wasserstein_D.item())
                else:
                    break  # 达到cfg.Critic_Iter，结束该次循环

            # ********************************
            # train Generation
            for p in netD.parameters():  # reset requires_grad
                p.requires_grad = False  # avoid computation
            # GENERATOR
            netG.zero_grad()
            input_semv = Variable(input_sem)
            noise.normal_(0, 1)
            noisev = Variable(noise)
            # pdb.set_trace()

            fake = netG(noisev, input_semv)
            criticG_fake = netD(fake, input_semv)
            criticG_fake = criticG_fake.mean()
            G_cost = -criticG_fake
            # # classification loss
            label = torch.LongTensor(label).to(used_device)
            c_errG = cls_criterion(pretrain_cls.model(fake), Variable(label))  # 输入fake编码 与 label（指示下标）
            errG = G_cost + Cls_Weight * c_errG  # 同时对生成器的生成质量与分类质量进行误差判定。生成器的loss由判别器衡量

            train_loss_g.append(errG.item())
            errG.backward()
            optimizerG.step()

        #     # ********************************
        train_loss_d = np.mean(train_loss_d)
        train_loss_g = np.mean(train_loss_g)
        WD = np.mean(WD)

        fmt = '{:.4f}'.format
        pbar.set_description('Train Epoch: {}'.format(epoch))
        pbar.set_postfix(loss=fmt(WD), loss_d=fmt(train_loss_d), loss_g=fmt(train_loss_g))

        if epoch > 40 and (train_loss_d < min_loss_d or train_loss_g < min_loss_g):
            if train_loss_d < min_loss_d:
                min_loss_d = train_loss_d

            if train_loss_g < min_loss_g:
                min_loss_g = train_loss_g

        #     torch.save(netG.state_dict(), g_path_save)
        #     torch.save(netD.state_dict(), d_path_save)
    print(f"min loss D = {min_loss_d}, min loss G = {min_loss_g}")

    return pretrain_cls, netD, netG


def load_cgan(
        res_list: list,
        load_path: str = '',
        SemSize: int = 256,
        NoiseSize: int = 256,
        used_device: str = 'cuda:0'
):
    """
    Loading trained CGAN model
    :param res_list: the result list prepared by `prepare_generate` function
    :param load_path: the load directory
    :param SemSize: dim of embedding of class, eg:256
    :param NoiseSize: dim of noise, eg:256
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: trained CGAN model
    """

    netG = MLP_G(SemSize=SemSize, NoiseSize=NoiseSize, NGH=4096, FeaSize=res_list[0]).to(used_device)
    # netD = MLP_CRITIC(SemSize=SemSize, NDH=4096, FeaSize=res_list[0]).to(used_device)

    netG.load_state_dict(torch.load(load_path, map_location=used_device))
    # netD.load_state_dict(torch.load(load_path, map_location=used_device))

    return netG


def generate_cgan(
        netG,
        pretrain_cls,
        res_list: list,
        NoiseSize: int = 256,
        used_device: str = 'cuda:0'
):
    """
    Generate single cell expression profiles
    :param netG: trained netG model
    :param pretrain_cls: pretain CLASSIFIER class
    :param res_list: the result list prepared by `prepare_generate` function
    :param NoiseSize: dim of noise, eg:256
    :param used_device: the id of used device, 'cuda:0, 1 ...' or 'cpu'
    :return: generated single cell expression profiles and metadata
    """
    cell_number_target_num = res_list[5]

    generate_label_list = []
    for key in cell_number_target_num.keys():
        generate_label_list += [key] * cell_number_target_num[key]

    label = np.array(generate_label_list)

    # netG in cuda now
    for p in netG.parameters():  # reset requires_grad
        p.requires_grad = False  # avoid computation

    dataloader = DataLoader(labelDataset(label=label), batch_size=1024, shuffle=False)
    ratio = 1
    netG.eval()
    netG.to(used_device)

    pretrain_cls.model.to(used_device)
    cell_all_generate = []
    label_all_generate = []

    with torch.no_grad():
        with tqdm(total=len(dataloader.dataset)) as pbar:
            for batch_idx, data in enumerate(dataloader):  # 一个batch
                label = data
                noise = torch.FloatTensor(label.shape[0], NoiseSize).to(used_device)  # (32, 256)

                input_sem = pretrain_cls.classifier_embed(label).to(used_device)
                input_semv = Variable(input_sem)

                label = label.cpu().numpy()
                # cell_feature = torch.tensor(cell_feature, dtype=torch.float32).cuda()
                for j in range(ratio):  # 翻倍多少
                    noise.normal_(0, 1)  # 每次都重新随机noise
                    noisev = Variable(noise)
                    fake = netG(noisev, input_semv)  # 固定输入的语义

                    fake = fake.cpu().data.numpy()
                    cell_all_generate.extend(fake)
                    label_all_generate.extend(label)

                pbar.update(label.shape[0])

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
            celltype_key: str = 'Cell_type',
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
        :param celltype_key: column name of celltype
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
        grid_out[celltype_key] = "unassigned"

        # get proportion of each cell type from generate_meta
        cell_type_num_dict = dict(Counter(generate_meta[celltype_key].values))
        cell_num = generate_meta.shape[0]

        for k, v in cell_type_num_dict.items():
            prop = v / cell_num
            dic1 = {key: value for key, value in sample_dic.items() if value < np.quantile(sample, prop)}
            grid_out.loc[dic1.keys(), [celltype_key]] = k
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
        cell_type_name = np.unique(generate_meta[celltype_key].values)
        replace_num = len(grid_out.loc[grid_out[celltype_key] == 'unassigned', [celltype_key]])
        grid_out.loc[grid_out[celltype_key] == 'unassigned', [celltype_key]] = np.random.choice(cell_type_name,
                                                                                                replace_num)

        if spatial_cell_type is None:
            print('generating spatial patterns of totally ' + str(len(cell_type_num_dict)) + ' cell types...')
            spa_pattern = pd.merge(grid_out, sim_point, how='inner')
        else:
            choose_cell_type_list = spatial_cell_type
            print('generating spatial patterns of selected ' + str(len(choose_cell_type_list)) + ' cell types...')
            spa_pattern = pd.merge(grid_out, sim_point, how='inner')
            # random shuffle other cell type coordinates
            idx = spa_pattern.loc[~spa_pattern[celltype_key].isin(choose_cell_type_list), ].index
            idx_cell_type = spa_pattern.loc[
                ~spa_pattern[celltype_key].isin(choose_cell_type_list), celltype_key].values.tolist()
            random.shuffle(idx_cell_type)
            spa_pattern.loc[idx, [celltype_key]] = idx_cell_type

        # correct cell type prop
        generate_cell_type_num_dict = dict(Counter(spa_pattern[celltype_key].values))
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
                np.random.choice(spa_pattern[spa_pattern[celltype_key] == k].index.values, v, replace=False).tolist())

        add_ct_list = []
        for k, v in add_dict.items():
            add_ct_list.extend([k] * v)

        random.shuffle(add_ct_list)
        spa_pattern.loc[correct_index, [celltype_key]] = add_ct_list

        # pdb.set_trace()

        # lamda (0-1): large--clear pattern small--fuzzy pattern
        lamda_len = round(len(spa_pattern) * (1 - lamda))
        select_idx = np.random.choice(spa_pattern.index.values, lamda_len, replace=False)
        change_idx = np.random.choice(select_idx, len(select_idx), replace=False)
        spa_pattern.loc[select_idx, [celltype_key]] = np.array(spa_pattern.loc[change_idx, celltype_key])

        # match Cell - Cell_type - coord
        generate_meta.sort_values(by=[celltype_key], inplace=True)
        spa_pattern.sort_values(by=[celltype_key], inplace=True)
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


class SPatternGeneratorSubtype:
    def __init__(self):
        pass

    def run(self,
            generate_meta: DataFrame,
            set_seed: bool = False,
            seed: int = 12345,
            celltype_key: str = 'Cell_type',
            select_cell_type: str = '',
            subtype_key: str = '',
            spatial_dim: int = 2,
            delta: float = 25,
            ):

        """
            Generate spatial pattern of cell types with single cell resolution
            :param generate_meta: generate single cell meta from VAE
            :param set_seed: False--random True--set seed
            :param seed: random seed
            :param celltype_key: column name of celltype
            :param select_cell_type: select cell types to generate subtype spatial patterns
            :param spatial_dim: spatial pattern dim
            :param subtype_key: columns of subtype
            :param delta: spatial pattern delta, large--big spatial pattern
            :return: DataFrame of subtype meta with spatial coordinates
            """

        if set_seed:
            np.random.seed(seed)

        print('generating subtype spatial patterns of ' + select_cell_type, '...')
        generate_meta_sub = generate_meta[generate_meta[celltype_key] == select_cell_type]
        if int(spatial_dim) == 2:
            coord = generate_meta_sub[['point_x', 'point_y']]
        elif int(spatial_dim) == 3:
            coord = generate_meta_sub[['point_x', 'point_y', 'point_z']]

        # create distance matrix
        distance = pd.DataFrame(np.array(self.__cal_dist(coord=np.array(coord))))
        # Generate random variable
        mean = np.repeat(0, distance.shape[0])
        cov = np.exp(-(distance.values / delta) ** 2)
        sample = np.random.multivariate_normal(mean, cov, tol=1e-6)
        sample_dic = {}
        for i in range(len(sample)):
            sample_dic[i] = sample[i]

        coord[subtype_key] = 'unassigned'
        coord = coord.reset_index()

        # get proportion of each cell type from generate_meta
        cell_type_num_dict = dict(Counter(generate_meta_sub[subtype_key].values))
        cell_num = generate_meta_sub.shape[0]

        for k, v in cell_type_num_dict.items():
            prop = v / cell_num
            dic1 = {key: value for key, value in sample_dic.items() if value < np.quantile(sample, prop)}
            coord.loc[dic1.keys(), [subtype_key]] = k
            sample_dic = {key: value for key, value in sample_dic.items() if value >= np.quantile(sample, prop)}
            sample = np.array(np.array(list(sample_dic.values())))
            cell_num = cell_num - v

        # handle unassigned
        cell_type_name = np.unique(generate_meta_sub[subtype_key].values)
        replace_num = len(coord.loc[coord[subtype_key] == 'unassigned', [subtype_key]])
        coord.loc[coord[subtype_key] == 'unassigned', [subtype_key]] = np.random.choice(cell_type_name, replace_num)

        # correct cell type prop
        generate_cell_type_num_dict = dict(Counter(coord[subtype_key].values))
        add_dict = {}
        sub_dict = {}
        for i in list(cell_type_name):
            if generate_cell_type_num_dict[i] < cell_type_num_dict[i]:
                add_dict[i] = cell_type_num_dict[i] - generate_cell_type_num_dict[i]
            elif generate_cell_type_num_dict[i] > cell_type_num_dict[i]:
                sub_dict[i] = generate_cell_type_num_dict[i] - cell_type_num_dict[i]

        correct_index = []
        for k, v in sub_dict.items():
            correct_index.extend(
                np.random.choice(coord[coord[subtype_key] == k].index.values, v, replace=False).tolist())

        add_ct_list = []
        for k, v in add_dict.items():
            add_ct_list.extend([k] * v)

        random.shuffle(add_ct_list)
        coord.loc[correct_index, [subtype_key]] = add_ct_list

        # adjust cell order
        spa_pattern = generate_meta_sub
        coord = coord.sort_values(subtype_key)
        spa_pattern = spa_pattern.sort_values(subtype_key)
        spa_pattern['point_x'] = coord['point_x'].values
        spa_pattern['point_y'] = coord['point_y'].values

        return spa_pattern

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
        count_threshold: Optional[float] = None,
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
    :param count_threshold: sparsity calibration, reset the expression values below this threshold to zero
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

    generate_meta_columns = list(generate_meta.columns)
    spot_x_loc = generate_meta_columns.index('spot_x')
    spot_y_loc = generate_meta_columns.index('spot_y')

    for i in range(0, generate_meta.shape[0]):
        # spot_x TODO:
        generate_meta.iloc[i, spot_x_loc] = select_min_grid(generate_meta.iloc[i, 2], grid_seq)
        # spot_y
        generate_meta.iloc[i, spot_y_loc] = select_min_grid(generate_meta.iloc[i, 3], grid_seq)

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

    if count_threshold is not None:
        spot_data = spot_data.apply(lambda x: (np.where(x < count_threshold, 0, x)))

    return spot_data, spot_meta, generate_meta


def generate_image_data(
        generate_meta: DataFrame,
        generate_data: DataFrame,
        gene_type: str,
        min_cell: int,
        n_gene: int,
        count_threshold: Optional[float] = None,
):
    """
    Generate image-based spatial transcriptomics data
    :param generate_meta: generate single cell meta from VAE
    :param generate_data: generate single cell data from VAE
    :param gene_type: targeted genes type
    :param min_cell: filter the genes expressed in fewer than `min_cell` cells
    :param n_gene: targeted genes number
    :param count_threshold: sparsity calibration, reset the expression values below this threshold to zero
    :return:
        target_data: image-based spatial transcriptomics data
        generate_meta: image-based spatial transcriptomics meta
        generate_data: ground truth of untargeted genes
    """
    print('generating image-based ST data with ' + str(n_gene) + ' targeted genes...')
    # check gene type
    assert gene_type in ['whole', 'hvg', 'marker', 'random'], \
        "the gene type is supposed to whole, hvg, marker, or random!!"
    # check gene number
    if n_gene is not None:
        assert n_gene <= generate_data.shape[0], \
            "the number of selected genes is too large, please reduce and try again!!"

    # gene type
    if gene_type == 'whole':
        print('generating with whole genes...')
        target_data = generate_data
    else:
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

    if count_threshold is not None:
        target_data = target_data.apply(lambda x: (np.where(x < count_threshold, 0, x)))

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

def read_or_pass_through(input_data):
    """
    Helper function to either read a CSV file or pass through data if it's already a DataFrame.
    Additionally checks if the DataFrame contains columns with 'x' and 'y', which might be needed for spatial data.

    Example usage:
    df = read_or_pass_through("path/to/your/data.csv")
    df = read_or_pass_through(your_dataframe)
    """
    if isinstance(input_data, str):
        df = pd.read_csv(input_data, index_col=0)
    elif isinstance(input_data, pd.DataFrame):
        df = input_data
    else:
        raise ValueError("Input data must be either a path to a CSV file or a pandas DataFrame.")

    columns = df.columns
    has_x = any('x' in col for col in columns)
    has_y = any('y' in col for col in columns)

    if not has_x or not has_y:
        print("Warning: The input data may be missing spatial coordinates as it lacks columns containing 'x' and/or 'y'.")

    return df

def create_and_save_anndata(data, meta, spatial_key=['x', 'y'], save_path=None):
    """
    Creates an Anndata object from data and meta. Optionally saves to a file.

    Args:
    data: str or pd.DataFrame - Path to the data CSV file or a DataFrame containing the data.
    meta: str or pd.DataFrame - Path to the meta CSV file or a DataFrame containing the metadata.
    spatial_key: list - List of column names in meta that correspond to spatial coordinates.
    save_path: str - Path to save the Anndata object file. If None, the file is not saved.

    Returns:
    adata: Anndata object
    """
    # Read or pass through data and meta
    data_df = read_or_pass_through(data)
    meta_df = read_or_pass_through(meta)

    # Create an AnnData object from the transposed data
    adata = sc.AnnData(data_df.T)

    # Assign meta as the .obs attribute of AnnData
    adata.obs = meta_df

    # Set spatial coordinates if the appropriate columns are available in meta
    if all(key in meta_df for key in spatial_key):
        adata.obsm['spatial'] = np.array(meta_df[spatial_key])
    else:
        raise ValueError(f"Meta data must contain the specified spatial keys: {spatial_key}")

    # Convert matrix to a sparse format if necessary
    if isinstance(adata.X, np.ndarray):
        adata.X = csr_matrix(adata.X)

    # Save the AnnData object if a save path is provided
    if save_path:
        adata.write_h5ad(save_path)

    return adata


def scale_spatial_coordinates(adata=None, spatial_size=None, spatial_key=None):
    """
    Scales the spatial coordinates of an AnnData object to a specified size.

    Args:
    adata: AnnData object - The AnnData object containing spatial coordinates.
    spatial_size: int - The desired size for scaling spatial coordinates.
    spatial_key: list - List of column names in adata.obs that correspond to spatial coordinates.

    Returns:
    adata: AnnData object - The AnnData object with scaled spatial coordinates.

    Example:
    adata = create_and_save_anndata(data, meta, spatial_key=['x', 'y'], save_path=None)
    adata_scaled = sizer(adata, spatial_size=1000)

    Or with manually specified column names:
    adata_scaled = sizer(adata, spatial_size=1000, spatial_key=['x_coord', 'y_coord'])
    """
    if adata is None:
        raise ValueError("Please provide an AnnData object.")

    if spatial_size is None:
        raise ValueError("Please provide a spatial_size for scaling.")

    if 'spatial' not in adata.obsm:
        raise KeyError(
            "Key 'spatial' not found in adata.obsm. Please make sure the AnnData object contains spatial coordinates.")

    points = adata.obsm['spatial']
    min_vals = points.min(axis=0)
    max_vals = points.max(axis=0)
    normalized_points = (points - min_vals) / (max_vals - min_vals)
    scaled_points = np.floor(normalized_points * spatial_size).astype(int)
    adata.obsm['spatial'] = scaled_points

    if spatial_key is None:
        spatial_key = ['x_pixel', 'y_pixel']

    if len(spatial_key) != 2:
        raise ValueError("spatial_key must contain exactly two elements for 'x' and 'y' coordinates.")

    for key in spatial_key:
        if key not in adata.obs.columns:
            adata.obs[key] = scaled_points[:, spatial_key.index(key)]
        else:
            adata.obs[key] = scaled_points[:, spatial_key.index(key)]

    return adata

# generate custom spatial pattern
# (biologically interpretable, such as a spherical tumor tissue surrounded by normal tissues)
class SPatternGeneratorCustom:
    def __init__(self, sc_adata: AnnData, cell_num: int = 5000, celltype_key: str = 'Cell_type', set_seed: bool = False,
                 seed: int = 12345, spatial_size: int = 30, select_celltype: Optional[list] = None):
        self.sc_adata = sc_adata
        self.cell_num = cell_num
        self.celltype_key = celltype_key
        self.set_seed = set_seed
        self.seed = seed
        self.spatial_size = spatial_size
        self.celltype_list = list(set(self.sc_adata.obs[self.celltype_key]))
        self.select_celltype = select_celltype
        if self.select_celltype is None:
            self.select_celltype = self.celltype_list
        else:
            assert len(set(self.select_celltype) - set(self.celltype_list)) == 0, \
                "select cell types must exist in sc_adata!"

    def simulate_mixing(self, prop_list: Optional[list] = None):

        print('generating unstructured mixed spatial patterns...')
        if self.set_seed:
            np.random.seed(self.seed)
            random.seed(self.seed)

        if prop_list is None:
            print('no `prop_list` is provided, each cell type follows an equal proportion...')
            prop_list = [1 / len(self.select_celltype)] * len(self.select_celltype)

        else:
            assert len(prop_list) == len(self.select_celltype), \
                "the length of `prop_list` must be equal to the length of `select_celltype`!"
            prop_list = [x / sum(prop_list) for x in prop_list]

        celltype_num = [np.round(x * self.cell_num).astype(int) for x in prop_list]

        # Generates a realisation of the Poisson point process
        sim_point = self.__create_sim_pois(cell_num=sum(celltype_num)) * self.spatial_size

        # cell type list
        celltype_list_all = [ct for ct, num in zip(self.select_celltype, celltype_num) for _ in range(num)]
        # random shuffle
        random.shuffle(celltype_list_all)

        sim_point[self.celltype_key] = celltype_list_all

        return sim_point

    def simulate_cluster(self,
                         shape_list: list = ['Circle', 'Oval'],
                         cluster_celltype_list: list = [],
                         cluster_purity_list: list = [],
                         infiltration_celltype_list: list = [[]],
                         infiltration_prop_list: list = [[]],
                         background_celltype: list = [],
                         background_prop: Optional[list] = None,
                         center_x_list: list = [20, 10],
                         center_y_list: list = [20, 10],
                         a_list: list = [15, 20],
                         b_list: list = [10, 15],
                         theta_list: list = [np.pi / 4, np.pi / 4],
                         scale_value_list: list = [4.8, 4.8],
                         twist_value_list: list = [0.5, 0.5]):

        print('generating structured cluster spatial patterns...')
        if self.set_seed:
            np.random.seed(self.seed)
            random.seed(self.seed)

        shapes = ['Circle', 'Oval', 'Irregular']
        assert len(set(shape_list) - set(shapes)) == 0, \
            "the shape in `shape_list` must be `Circle`, `Oval`, or `Irregular`!"

        assert len(set(cluster_celltype_list) - set(self.select_celltype)) == 0, \
            "the cluster cell type in `cluster_celltype_list` must exist in `select_celltype`!"

        infiltration_list_set = set()
        for row in infiltration_celltype_list:
            infiltration_list_set.update(row)
        assert len(infiltration_list_set - set(self.select_celltype)) == 0, \
            "the infiltrated cell type in `infiltration_celltype_list` must exist in `select_celltype`!"

        assert len(set(background_celltype) - set(self.select_celltype)) == 0, \
            "`background_celltype` must exist in `select_celltype`!"

        # Generates a realisation of the Poisson point process
        sim_point = self.__create_sim_pois(cell_num=self.cell_num) * self.spatial_size
        sim_point[self.celltype_key] = 'unassigned'

        var_len_list = [len(shape_list), len(cluster_celltype_list), len(cluster_purity_list),
                        len(infiltration_celltype_list), len(infiltration_prop_list), len(center_x_list),
                        len(center_y_list), len(a_list), len(b_list), len(theta_list),
                        len(scale_value_list), len(twist_value_list)]
        assert len(set(var_len_list)) == 1, \
            "`the length of `shape_list`, `cluster_celltype_list`, `cluster_purity_list`, " \
            "`infiltration_celltype_list`, `infiltration_prop_list`, `center_x_list`, `center_y_list`,  " \
            "`a_list`, `b_list`, `theta_list`, `scale_value_list`, and `twist_value_list` must be equal!"

        for n in range(len(shape_list)):
            shape = shape_list[n]
            cluster_celltype = cluster_celltype_list[n]
            cluster_purity = cluster_purity_list[n]
            infiltration_prop = infiltration_prop_list[n]
            infiltration_celltype = infiltration_celltype_list[n]
            center_x = center_x_list[n]
            center_y = center_y_list[n]
            a = a_list[n]
            b = b_list[n]
            theta = theta_list[n]
            scale_value = scale_value_list[n]
            twist_value = twist_value_list[n]
            idx = []

            # Cluster
            for i in range(sim_point.shape[0]):
                x = sim_point.loc[i, 'point_x']
                y = sim_point.loc[i, 'point_y']
                x = x - center_x
                y = y - center_y

                if shape == 'Circle':
                    assert a == b, 'the value of `a` and `b` must be same when `shape == Circle`!'
                    # 使用椭圆方程进行判断
                    value = (x ** 2) / (a ** 2) + (y ** 2) / (b ** 2)
                    if value < 1:
                        sim_point.loc[i, self.celltype_key] = cluster_celltype
                        idx.append(i)

                elif shape == 'Oval':
                    # theta = np.pi / theta
                    # 旋转点坐标
                    x_rotated = x * np.cos(theta) + y * np.sin(theta)
                    y_rotated = -x * np.sin(theta) + y * np.cos(theta)
                    # 使用椭圆方程进行判断
                    value = (x_rotated ** 2) / (a ** 2) + (y_rotated ** 2) / (b ** 2)
                    if value < 1:
                        sim_point.loc[i, self.celltype_key] = cluster_celltype
                        idx.append(i)

                elif shape == 'Irregular':
                    x_rotated = x * np.cos(theta) + y * np.sin(theta)
                    y_rotated = -x * np.sin(theta) + y * np.cos(theta)
                    # t = np.arctan2(x, y)
                    # r = np.sqrt(x ** 2 + y ** 2) / scale_value
                    # equation = (r <= 16 * np.sin(t) ** 3 * twist_value) and (r >= 13 * np.cos(t))
                    # if equation:
                    #     sim_point.loc[i, self.celltype_key] = cluster_celltype
                    t = np.arctan2(x_rotated, y_rotated)
                    r = np.sqrt(x_rotated ** 2 + y_rotated ** 2) / scale_value
                    # equation = (r <= 16 * np.sin(t) ** 3 * twist_value) and (r >= 13 * np.cos(t))
                    equation = (r <= 16 * np.sin(t) ** 3 * twist_value) and (r >= (13 * np.cos(t) - 5 * np.cos(2 * t) -
                                                                                   2 * np.cos(3 * t) - np.cos(4 * t)))
                    if equation:
                        sim_point.loc[i, self.celltype_key] = cluster_celltype
                        idx.append(i)

            # Infiltration
            if cluster_purity is None:
                cluster_purity = 1
            else:
                assert 0 <= cluster_purity <= 1, '`cluster_purity` must in [0, 1]!'

            # avoid duplication in Infiltration
            sim_point_rep = sim_point.loc[idx].copy()
            cluster_num_all = sim_point_rep[sim_point_rep[self.celltype_key] == cluster_celltype].shape[0]
            cluster_num_inner = np.round(cluster_num_all * cluster_purity).astype(int)
            cluster_num_outer = cluster_num_all - cluster_num_inner
            if infiltration_prop is None:
                print('no `infiltration_prop` is provided, each infiltrated cell type follows an equal proportion...')
                infiltration_prop = [1 / len(infiltration_celltype)] * len(infiltration_celltype)
            else:
                assert len(infiltration_prop) == len(infiltration_celltype), \
                    "the length of `infiltration_prop` must be equal to the length of `infiltration_celltype`!"
                infiltration_prop = [x / sum(infiltration_prop) for x in infiltration_prop]

            infiltration_num_list = [np.round(x * cluster_num_outer).astype(int) for x in infiltration_prop]
            for i in range(len(infiltration_num_list)):
                tmp_idx = np.random.choice(sim_point_rep[sim_point_rep[self.celltype_key] == cluster_celltype].index.values,
                                           infiltration_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = infiltration_celltype[i]

        # background
        if background_prop is None:
            print('no `background_prop` is provided, each background cell type follows an equal proportion...')
            background_prop = [1 / len(background_celltype)] * len(background_celltype)
        else:
            assert len(background_prop) == len(background_celltype), \
                "the length of `background_prop` must be equal to the length of `background_celltype`!"
            background_prop = [x / sum(background_prop) for x in background_prop]

        background_num_all = sim_point[sim_point[self.celltype_key] == 'unassigned'].shape[0]
        background_num_list = [np.round(x * background_num_all).astype(int) for x in background_prop]

        for i in range(len(background_num_list)):
            if i < len(background_num_list) - 1:
                tmp_idx = np.random.choice(sim_point[sim_point[self.celltype_key] == 'unassigned'].index.values,
                                           background_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = background_celltype[i]
            else:
                sim_point.loc[sim_point[self.celltype_key] == 'unassigned', self.celltype_key] = background_celltype[i]

        return sim_point

    def simulate_ring(self,
                      shape_list: list = ['Circle', 'Oval'],
                      ring_celltype_list: list = [[]],
                      ring_purity_list: list = [],
                      infiltration_celltype_list: list = [[]],
                      infiltration_prop_list: list = [[]],
                      background_celltype: list = [],
                      background_prop: Optional[list] = None,
                      center_x_list: list = [20, 10],
                      center_y_list: list = [20, 10],
                      ring_width_list: list = [[2, 3], [2]],
                      a_list: list = [15, 20],
                      b_list: list = [10, 15],
                      theta_list: list = [np.pi / 4, np.pi / 4], ):

        print('generating structured immune ring spatial patterns...')
        if self.set_seed:
            np.random.seed(self.seed)
            random.seed(self.seed)

        shapes = ['Circle', 'Oval']
        assert len(set(shape_list) - set(shapes)) == 0, \
            "the shape in `shape_list` must be `Circle`, or `Oval`!"

        ring_list_set = set()
        for row in ring_celltype_list:
            ring_list_set.update(row)
        assert len(set(ring_list_set) - set(self.select_celltype)) == 0, \
            "the ring cell type in `ring_celltype_list` must exist in `select_celltype`!"

        infiltration_list_set = set()
        for row in infiltration_celltype_list:
            infiltration_list_set.update(row)
        assert len(infiltration_list_set - set(self.select_celltype)) == 0, \
            "the infiltrated cell type in `infiltration_celltype_list` must exist in `select_celltype`!"

        assert len(set(background_celltype) - set(self.select_celltype)) == 0, \
            "`background_celltype` must exist in `select_celltype`!"

        # Generates a realisation of the Poisson point process
        sim_point = self.__create_sim_pois(cell_num=self.cell_num) * self.spatial_size
        sim_point[self.celltype_key] = 'unassigned'

        var_len_list = [len(shape_list), len(ring_celltype_list), len(ring_purity_list),
                        len(infiltration_celltype_list), len(infiltration_prop_list), len(center_x_list),
                        len(center_y_list), len(ring_width_list),
                        len(a_list), len(b_list), len(theta_list)]
        assert len(set(var_len_list)) == 1, \
            "`the length of `shape_list`, `ring_celltype_list`, `ring_purity_list`, " \
            "`infiltration_celltype_list`, `infiltration_prop_list`, `center_x_list`, `center_y_list`, " \
            "`ring_width_list`, `a_list`, `b_list`, and `theta_list` must be equal!"

        major_ring_type = []
        for row in ring_celltype_list:
            if len(row) == 2:
                major_ring_type.append(row[0])
            elif len(row) == 3:
                major_ring_type.append(row[0])
                major_ring_type.append(row[1])

        for n in range(len(shape_list)):
            shape = shape_list[n]
            ring_celltype = ring_celltype_list[n]
            ring_purity = ring_purity_list[n]
            infiltration_prop = infiltration_prop_list[n]
            infiltration_celltype = infiltration_celltype_list[n]
            center_x = center_x_list[n]
            center_y = center_y_list[n]
            a = a_list[n]
            b = b_list[n]
            theta = theta_list[n]
            ring_width = ring_width_list[n]
            assert 1 <= len(ring_width) <= 2, \
                'the length of `ring_width` must be `1` for single ring or `2` for double rings!'
            if len(ring_width) == 1:
                assert len(ring_celltype) == 2, "the length of `ring_celltype` must be 2 single ring!"
            elif len(ring_width) == 2:
                assert len(ring_celltype) == 3, "the length of `ring_celltype` must be 3 for double rings!"

            idx = []

            # Ring
            for i in range(sim_point.shape[0]):
                x = sim_point.loc[i, 'point_x']
                y = sim_point.loc[i, 'point_y']
                x = x - center_x
                y = y - center_y

                if shape == 'Circle':
                    # 使用椭圆方程进行判断
                    assert a == b, 'the value of `a` and `b` must be same when `shape == Circle`!'
                    value_inner = (x ** 2) / (a ** 2) + (y ** 2) / (b ** 2)
                    if len(ring_width) == 1:
                        value_outer = (x ** 2) / ((a + ring_width[0]) ** 2) + (y ** 2) / ((b + ring_width[0]) ** 2)
                        if value_inner < 1:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[0]
                            idx.append(i)
                        elif value_outer < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[1]
                            idx.append(i)
                    elif len(ring_width) == 2:
                        value_middle = (x ** 2) / ((a + ring_width[0]) ** 2) + (y ** 2) / ((b + ring_width[0]) ** 2)
                        value_outer = (x ** 2) / ((a + (ring_width[0] + ring_width[1])) ** 2) + \
                                      (y ** 2) / ((b + (ring_width[0] + ring_width[1])) ** 2)
                        if value_inner < 1:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[0]
                            idx.append(i)
                        elif value_middle < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[1]
                            idx.append(i)
                        elif value_outer < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[2]
                            idx.append(i)

                elif shape == 'Oval':
                    # 旋转点坐标
                    x_rotated = x * np.cos(theta) + y * np.sin(theta)
                    y_rotated = -x * np.sin(theta) + y * np.cos(theta)
                    # 使用椭圆方程进行判断
                    value_inner = (x_rotated ** 2) / (a ** 2) + (y_rotated ** 2) / (b ** 2)
                    if len(ring_width) == 1:
                        value_outer = (x_rotated ** 2) / ((a + ring_width[0]) ** 2) + \
                                      (y_rotated ** 2) / ((b + ring_width[0]) ** 2)
                        if value_inner < 1:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[0]
                            idx.append(i)
                        elif value_outer < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[1]
                            idx.append(i)
                    elif len(ring_width) == 2:
                        value_middle = (x_rotated ** 2) / ((a + ring_width[0]) ** 2) + \
                                       (y_rotated ** 2) / ((b + ring_width[0]) ** 2)
                        value_outer = (x_rotated ** 2) / ((a + (ring_width[0] + ring_width[1])) ** 2) + \
                                      (y_rotated ** 2) / ((b + (ring_width[0] + ring_width[1])) ** 2)
                        if value_inner < 1:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[0]
                            idx.append(i)
                        elif value_middle < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[1]
                            idx.append(i)
                        elif value_outer < 1 and sim_point.loc[i, self.celltype_key] not in major_ring_type:
                            sim_point.loc[i, self.celltype_key] = ring_celltype[2]
                            idx.append(i)

            # Infiltration
            if ring_purity is None:
                ring_purity = 1
            else:
                assert 0 <= ring_purity <= 1, '`cluster_purity` must in [0, 1]!'

            # avoid duplication in Infiltration
            sim_point_rep = sim_point.loc[idx].copy()
            ring_num_all = sim_point_rep[sim_point_rep[self.celltype_key].isin(ring_celltype)].shape[0]
            ring_num_inner = np.round(ring_num_all * ring_purity).astype(int)
            ring_num_outer = ring_num_all - ring_num_inner
            if infiltration_prop is None:
                print(
                    'no `infiltration_prop` is provided, each infiltrated cell type follows an equal proportion...')
                infiltration_prop = [1 / len(infiltration_celltype)] * len(infiltration_celltype)
            else:
                assert len(infiltration_prop) == len(infiltration_celltype), \
                    "the length of `infiltration_prop` must be equal to the length of `infiltration_celltype`!"
                infiltration_prop = [x / sum(infiltration_prop) for x in infiltration_prop]

            infiltration_num_list = [np.round(x * ring_num_outer).astype(int) for x in infiltration_prop]
            for i in range(len(infiltration_num_list)):
                tmp_idx = np.random.choice(
                    sim_point_rep[sim_point_rep[self.celltype_key].isin(ring_celltype)].index.values,
                    infiltration_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = infiltration_celltype[i]

        # background
        if background_prop is None:
            print('no `background_prop` is provided, each background cell type follows an equal proportion...')
            background_prop = [1 / len(background_celltype)] * len(background_celltype)
        else:
            assert len(background_prop) == len(background_celltype), \
                "the length of `background_prop` must be equal to the length of `background_celltype`!"
            background_prop = [x / sum(background_prop) for x in background_prop]

        background_num_all = sim_point[sim_point[self.celltype_key] == 'unassigned'].shape[0]
        background_num_list = [np.round(x * background_num_all).astype(int) for x in background_prop]

        for i in range(len(background_num_list)):
            if i < len(background_num_list) - 1:
                tmp_idx = np.random.choice(sim_point[sim_point[self.celltype_key] == 'unassigned'].index.values,
                                           background_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = background_celltype[i]
            else:
                sim_point.loc[sim_point[self.celltype_key] == 'unassigned', self.celltype_key] = background_celltype[i]

        return sim_point

    def simulate_stripes(self,
                         y1_list: list = [None, None],
                         y2_list: list = [None, None],
                         stripe_width_list: list = [2, 3],
                         stripe_purity_list: list = [],
                         stripe_celltype_list: list = [],
                         infiltration_celltype_list: list = [[]],
                         infiltration_prop_list: list = [[]],
                         background_celltype: list = [],
                         background_prop: Optional[list] = None,
                         ):

        print('generating structured stripes spatial patterns...')
        if self.set_seed:
            np.random.seed(self.seed)
            random.seed(self.seed)

        assert len(set(stripe_celltype_list) - set(self.select_celltype)) == 0, \
            "the stripe cell type in `stripe_celltypee_list` must exist in `select_celltype`!"

        infiltration_list_set = set()
        for row in infiltration_celltype_list:
            infiltration_list_set.update(row)
        assert len(infiltration_list_set - set(self.select_celltype)) == 0, \
            "the infiltrated cell type in `infiltration_celltype_list` must exist in `select_celltype`!"

        assert len(set(background_celltype) - set(self.select_celltype)) == 0, \
            "`background_celltype` must exist in `select_celltype`!"

        # Generates a realisation of the Poisson point process
        sim_point = self.__create_sim_pois(cell_num=self.cell_num) * self.spatial_size
        sim_point[self.celltype_key] = 'unassigned'

        var_len_list = [len(stripe_width_list), len(stripe_purity_list), len(stripe_celltype_list),
                        len(infiltration_celltype_list), len(infiltration_prop_list)]
        assert len(set(var_len_list)) == 1, \
            "`the length of `stripe_width_list`, `stripe_purity_list`, `stripe_celltypee_list`, " \
            "`infiltration_celltype_list`, and `infiltration_prop_list` must be equal!"

        for n in range(len(stripe_celltype_list)):
            y1 = y1_list[n]
            y2 = y2_list[n]
            stripe_celltype = stripe_celltype_list[n]
            stripe_purity = stripe_purity_list[n]
            stripe_width = stripe_width_list[n]
            infiltration_prop = infiltration_prop_list[n]
            infiltration_celltype = infiltration_celltype_list[n]
            idx = []

            x1 = min(sim_point['point_x'])
            x2 = max(sim_point['point_x'])
            if y1 is None:
                y1 = random.uniform(min(sim_point['point_y']), max(sim_point['point_y']))
            if y2 is None:
                y2 = random.uniform(min(sim_point['point_y']), max(sim_point['point_y']))
            # 计算线条方向向量
            direction = np.array([x2 - x1, y2 - y1], dtype=float)
            # 计算线条长度
            length = np.linalg.norm(direction)
            # 归一化方向向量
            direction /= length

            # 计算线条的垂直向量
            perpendicular = np.array([-direction[1], direction[0]], dtype=float)

            # 计算线条的四个顶点
            p1 = np.array([x1, y1], dtype=float) + perpendicular * stripe_width / 2
            p2 = np.array([x1, y1], dtype=float) - perpendicular * stripe_width / 2
            # p3 = np.array([x2, y2], dtype=float) + perpendicular * stripe_width / 2
            # p4 = np.array([x2, y2], dtype=float) - perpendicular * stripe_width / 2

            # Stripes
            for i in range(sim_point.shape[0]):
                x = sim_point.loc[i, 'point_x']
                y = sim_point.loc[i, 'point_y']
                point = np.array([x, y], dtype=float)
                u = np.dot(point - p1, direction)
                v = np.dot(point - p2, perpendicular)
                if 0 <= u <= length and np.abs(v) <= stripe_width / 2:
                    sim_point.loc[i, self.celltype_key] = stripe_celltype
                    idx.append(i)

            # Infiltration
            if stripe_purity is None:
                stripe_purity = 1
            else:
                assert 0 <= stripe_purity <= 1, '`cluster_purity` must in [0, 1]!'

            # avoid duplication in Infiltration
            sim_point_rep = sim_point.loc[idx].copy()
            stripe_num_all = sim_point_rep[sim_point_rep[self.celltype_key] == stripe_celltype].shape[0]
            stripe_num_inner = np.round(stripe_num_all * stripe_purity).astype(int)
            stripe_num_outer = stripe_num_all - stripe_num_inner
            if infiltration_prop is None:
                print('no `infiltration_prop` is provided, each infiltrated cell type follows an equal proportion...')
                infiltration_prop = [1 / len(infiltration_celltype)] * len(infiltration_celltype)
            else:
                assert len(infiltration_prop) == len(infiltration_celltype), \
                    "the length of `infiltration_prop` must be equal to the length of `infiltration_celltype`!"
                infiltration_prop = [x / sum(infiltration_prop) for x in infiltration_prop]

            infiltration_num_list = [np.round(x * stripe_num_outer).astype(int) for x in infiltration_prop]
            for i in range(len(infiltration_num_list)):
                tmp_idx = np.random.choice(sim_point_rep[sim_point_rep[self.celltype_key] == stripe_celltype].index.values,
                                           infiltration_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = infiltration_celltype[i]

        # background
        if background_prop is None:
            print('no `background_prop` is provided, each background cell type follows an equal proportion...')
            background_prop = [1 / len(background_celltype)] * len(background_celltype)
        else:
            assert len(background_prop) == len(background_celltype), \
                "the length of `background_prop` must be equal to the length of `background_celltype`!"
            background_prop = [x / sum(background_prop) for x in background_prop]

        background_num_all = sim_point[sim_point[self.celltype_key] == 'unassigned'].shape[0]
        background_num_list = [np.round(x * background_num_all).astype(int) for x in background_prop]

        for i in range(len(background_num_list)):
            if i < len(background_num_list) - 1:
                tmp_idx = np.random.choice(sim_point[sim_point[self.celltype_key] == 'unassigned'].index.values,
                                           background_num_list[i], replace=False).tolist()
                sim_point.loc[tmp_idx, self.celltype_key] = background_celltype[i]
            else:
                sim_point.loc[sim_point[self.celltype_key] == 'unassigned', self.celltype_key] = background_celltype[i]

        return sim_point

    def __create_sim_pois(self,
                          cell_num: int):
        """
        Generates a realisation of the Poisson point process
        :param cell_num: the number of points
        :return: DataFrame of spatial coordinates of points
        """

        xx = np.random.uniform(0, 1, cell_num)
        yy = np.random.uniform(0, 1, cell_num)
        out = pd.DataFrame({'point_x': xx, 'point_y': yy})

        return out


