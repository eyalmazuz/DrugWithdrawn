from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler
import numpy as np
import pandas as pd
import random
import os
import torch
from torch.utils.data import DataLoader

class DrugDataset(Dataset):
    
    def __init__(self, data):
        self.data = data

    def __getitem__(self, index):
        features = self.data[index][:-1]
        target = self.data[index][-1]
        return features, target

    def __len__(self):
        return len(self.data)


def get_class_distribution(obj):
    count_dict = {
        "not withdrawn": 0,
        "withdrawn": 0
    }
    
    for i in obj:
        if i == 0: 
            count_dict['not withdrawn'] += 1
        elif i == 1: 
            count_dict['withdrawn'] += 1             
        else:
            print("Check classes.")  
    return count_dict

def get_dataloader(dataset, X_train, y_train, X_test, y_test, batch_size=1, to_sample=True):

    """Create a dataset and split it into train/dev/test."""

    dataset_train = [dataset[i] for i in list(X_train['index'].values)]
    dataset_test = [dataset[i] for i in list(X_test['index'].values)]
    

    train_data = DrugDataset(dataset_train)
    test_data = DrugDataset(dataset_test)
    
    target_list = []
    for _, t in train_data:
        target_list.append(t.type(torch.LongTensor).item())
    class_count = [i for i in get_class_distribution(y_train).values()]
    class_weights = 1./torch.tensor(class_count, dtype=torch.float) 
    
    if to_sample:
        class_weights_all = class_weights[target_list]
        print('Class Weights : {}'.format(class_weights))
        print('Class Weights All: {}'.format(class_weights_all))
        weighted_sampler = WeightedRandomSampler(
        weights=class_weights_all,
        num_samples=len(class_weights_all),
        replacement=True)
        train_loader = DataLoader(train_data, batch_size = batch_size, sampler=weighted_sampler)
    else:
        train_loader = DataLoader(train_data, batch_size = batch_size, shuffle=True)
    
    test_loader = DataLoader(test_data, batch_size = batch_size, shuffle = True)
    
    return train_loader, test_loader, class_weights

def get_device(device_id):
    torch.cuda.set_device(device_id)
    """CPU or GPU."""
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print('The code uses GPU')
    else:
        device = torch.device('cpu')
        print('The code uses CPU')
    return device


def set_random_seed(seed_value, use_cuda=True):
    np.random.seed(seed_value)  # cpu vars
    torch.manual_seed(seed_value)  # cpu  vars
    random.seed(seed_value)  # Python
    os.environ['PYTHONHASHSEED'] = str(seed_value)  # Python hash buildin
    if use_cuda:
        torch.cuda.manual_seed(seed_value)
        torch.cuda.manual_seed_all(seed_value)  # gpu vars
        torch.backends.cudnn.deterministic = True  # needed
        torch.backends.cudnn.benchmark = False