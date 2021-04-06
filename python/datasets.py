##################################### datasets.py #######################################
# This file contains the class overrides of the Dataset superclass, specifying the
# training and testing datasets, along with the transformations they operate on the data.
#
# This file is distributed under the following license:
# =======================================================================================
# Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
# Written by Rémi Clerc <remi.clerc@idiap.ch>
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
# conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
# of conditions and the following disclaimer in the documentation and/or other materials
# provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# =======================================================================================

import h5py
import numpy as np
import torch
from torch.utils.data import Dataset
from transformations import RndVHFlip, poisson_noise


class TrainDataset(Dataset):
    def __init__(self, h5_file, augment=False, poisson=False, seed=None, photons=150, background=3):
        super(TrainDataset, self).__init__()
        self.h5_file = h5_file
        self.augment = augment
        self.poisson = poisson
        self.flip = RndVHFlip()
        self.photons = photons
        self.background = background
        self.seed = seed

    def __getitem__(self, idx):
        with h5py.File(self.h5_file, 'r') as f:
            
            # add empty dimensions for compatibility with the model
            item = np.expand_dims(f['lr'][idx] / 255., 0), np.expand_dims(f['hr'][idx] / 255., 0)

            # perform random flips
            if self.augment:
                item = self.flip(item)

            # add Poisson noise
            if self.poisson:
                item = poisson_noise(item[0], self.photons, self.background, self.seed), item[1]
            
            item = torch.from_numpy(item[0].copy()), torch.from_numpy(item[1].copy())
            return item

    def __len__(self):
        with h5py.File(self.h5_file, 'r') as f:
            return len(f['lr'])


class EvalDataset(Dataset):
    def __init__(self, h5_file, poisson=False, seed=None, photons=150, background=3):
        super(EvalDataset, self).__init__()
        self.h5_file = h5_file
        self.poisson = poisson
        self.photons = photons
        self.background = background
        self.seed = seed

    def __getitem__(self, idx):
        with h5py.File(self.h5_file, 'r') as f:

            # add empty dimensions for compatibility with the model
            item = np.expand_dims(f['lr'][idx][:, :] / 255., 0), np.expand_dims(f['hr'][idx][:, :] / 255., 0)

            # add poisson noise
            if self.poisson:
                item = poisson_noise(item[0], self.photons, self.background, self.seed), item[1]
                
            return item

    def __len__(self):
        with h5py.File(self.h5_file, 'r') as f:
            return len(f['lr'])
