import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.autograd import Variable
import torch.optim as optim


# ****************************************************************************
# ******************************** VAE ***************************************
# ****************************************************************************
class VAE(nn.Module):
    def __init__(self, embedding_size, hidden_size_list: list, mid_hidden):
        super(VAE, self).__init__()
        self.embedding_size = embedding_size
        self.hidden_size_list = hidden_size_list
        self.mid_hidden = mid_hidden

        self.enc_feature_size_list = [self.embedding_size] + self.hidden_size_list + [self.mid_hidden * 2]
        self.dec_feature_size_list = [self.embedding_size] + self.hidden_size_list + [self.mid_hidden]

        self.encoder = nn.ModuleList(
            [nn.Linear(self.enc_feature_size_list[i], self.enc_feature_size_list[i + 1]) for i in
             range(len(self.enc_feature_size_list) - 1)])
        self.decoder = nn.ModuleList(
            [nn.Linear(self.dec_feature_size_list[i], self.dec_feature_size_list[i - 1]) for i in
             range(len(self.dec_feature_size_list) - 1, 0, -1)])

    def encode(self, x):
        for i, layer in enumerate(self.encoder):
            x = self.encoder[i](x)
            if i != len(self.encoder) - 1:
                x = F.relu(x)
        return x

    def decode(self, x):
        for i, layer in enumerate(self.decoder):
            x = self.decoder[i](x)
            # if i != len(self.decoder) - 1:  # 考虑去除，防止负数出现
            x = F.relu(x)
        return x

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar) # var => std
        eps = torch.randn_like(std)
        return eps * std + mu

    def forward(self, x, used_device):
        x = x.to(used_device)
        encoder_output = self.encode(x)
        mu, sigma = torch.chunk(encoder_output, 2, dim=1)  # mu, log_var
        hidden = self.reparameterize(mu, sigma)
        x_hat = self.decode(hidden)
        kl_div = 0.5 * torch.sum(torch.exp(sigma) + torch.pow(mu, 2) - 1 - sigma) / (x.shape[0] * x.shape[1])
        return x_hat, kl_div


# ****************************************************************************
# ***************************** CLASSIFIER ***********************************
# ****************************************************************************
def weights_init(m):
    classname = m.__class__.__name__
    if classname.find('Linear') != -1:
        m.weight.data.normal_(0.0, 0.02)
        m.bias.data.fill_(0)
    elif classname.find('BatchNorm') != -1:
        m.weight.data.normal_(1.0, 0.02)
        m.bias.data.fill_(0)


class CLASSIFIER:
    # 先用全连接层分类器
    # 使用nn.NLLLoss()作为损失函数，配合logsoftmax实现对于类别的划分（类别概率大->全连接层输出高->softmax接近1->logsoftmax接近0->NLLLose小）
    def __init__(self, _input_dim, _nclass, dic, SemSize, dataloader='', used_device='', _lr=0.001, _nepoch=20):
        self.input_dim = _input_dim
        self.nclass = _nclass
        self.used_device = used_device
        self.model = LINEAR_LOGSOFTMAX(self.input_dim, self.nclass).to(self.used_device)
        self.model.apply(weights_init)
        self.criterion = nn.NLLLoss().to(self.used_device)
        self.dic = dic
        self.nepoch = _nepoch
        self.optimizer = optim.Adam(self.model.parameters(), lr=_lr, betas=(0.5, 0.999))
        self.dataloader = dataloader
        self.SemSize = SemSize

        self.one_hot_embed = torch.zeros(len(self.dic), self.SemSize)
        self._one_hot_construct()

        self._pretrain()

    def _pretrain(self):
        min_loss = 10000
        for epoch in range(self.nepoch):
            loss_item = []
            for batch_idx, data in enumerate(self.dataloader):
                cell_feature, label = data  # dim of cell sample cell ,label should map to embedding, which dim=dim_noise
                cell_feature = torch.tensor(cell_feature, dtype=torch.float32).to(self.used_device)
                label = torch.LongTensor(label).to(self.used_device)
                self.model.zero_grad()
                inputv = Variable(cell_feature)
                labelv = Variable(label)
                output = self.model(inputv)
                loss = self.criterion(output, labelv)  # 哪个是label哪个的logsoftmax就要小
                loss.backward()
                self.optimizer.step()
                loss_item.append(loss.item())

            loss_item = np.mean(loss_item)
            if loss_item < min_loss:
                min_loss = loss_item

        # torch.save(self.model.state_dict(),
        # osp.join(self.args.save, self.args.name + '_pretrain_classification_epoch_' +
        #          str(self.nepoch) + '.pth'))
        print("classification pretrain min loss:", min_loss)

    def _one_hot_construct(self):
        for i in range(len(self.dic)):
            self.one_hot_embed[i][i] = 1

    def classifier_embed(self, label):
        # 先用one-hot编码
        return self.one_hot_embed[label]


class LINEAR_LOGSOFTMAX(nn.Module):  # 使用logsoftmax作为分类结果。越接近0说明越真实
    def __init__(self, input_dim, nclass):
        super(LINEAR_LOGSOFTMAX, self).__init__()
        self.fc = nn.Linear(input_dim, nclass)
        self.logic = nn.LogSoftmax(dim=1)

    def forward(self, x):
        o = self.logic(self.fc(x))
        return o


# ****************************************************************************
# ******************************** MLP_G *************************************
# ****************************************************************************
class MLP_G(nn.Module):
    def __init__(self, SemSize, NoiseSize, NGH, FeaSize):
        super(MLP_G, self).__init__()
        self.fc1 = nn.Linear(SemSize + NoiseSize, NGH)
        self.fc2 = nn.Linear(NGH, FeaSize)
        self.lrelu = nn.LeakyReLU(0.2, True)
        # self.prelu = nn.PReLU()
        self.relu = nn.ReLU(True)

        self.apply(weights_init)

    def forward(self, noise, sem):
        h = torch.cat((noise, sem), 1)
        h = self.lrelu(self.fc1(h))
        h = self.relu(self.fc2(h))
        return h


# ****************************************************************************
# ****************************** MLP_CRITIC **********************************
# ****************************************************************************
class MLP_CRITIC(nn.Module):
    def __init__(self, SemSize, NDH, FeaSize):
        super(MLP_CRITIC, self).__init__()
        self.fc1 = nn.Linear(FeaSize + SemSize, NDH)
        # self.fc2 = nn.Linear(opt.ndh, opt.ndh)
        self.fc2 = nn.Linear(NDH, 1)
        self.lrelu = nn.LeakyReLU(0.2, True)

        self.apply(weights_init)

    def forward(self, x, sem):

        h = torch.cat((x, sem), 1)
        h = self.lrelu(self.fc1(h))
        h = self.fc2(h)
        return h