################################# train.py ##############################################
# This file contains the code to train SRCNN and the Autoencoder, and is adapted
# from https://github.com/Lornatang/SRCNN-PyTorch
#
# The original version of this file is distributed under the following license:
# =======================================================================================
# Copyright 2020 Dakewe Biotech Corporation. All Rights Reserved.
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =======================================================================================
#
# The modifications that were applied to the original code in this file are distributed
# under the following license:
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

import argparse
import math
import os
import copy
import random
import time

import torch.nn as nn
import torch.optim as optim
import torch.utils.data
import torch.utils.data.distributed
import torch.backends.cudnn as cudnn
from torch.utils.data.dataloader import DataLoader

import numpy as np

from AE import AE
from SRCNN import SRCNN
from datasets import TrainDataset, EvalDataset

######################################################################
# Helper taken from stackoverflow to allow for boolean arguments
# https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
######################################################################

##################### Args parsing ##########################
parser = argparse.ArgumentParser()
parser.add_argument('--train-file', type=str, required=True)    #path to h5 formatted train file
parser.add_argument('--eval-file', type=str, required=True)     #path to h5 formatted eval file
parser.add_argument('--outputs-dir', type=str, required=True)   #output directory
parser.add_argument('--model', type=str, required=True)         #model architecture, either "SRCNN" or "AE"
parser.add_argument('--lr', type=float, default=1e-4)           #learning rate
parser.add_argument('--batch-size', type=int, default=16)
parser.add_argument('--num-epochs', type=int, default=400)      #number of training epochs
parser.add_argument('--num-workers', type=int, default=8)
parser.add_argument('--seed', type=int, default=None)            #random seed
parser.add_argument('--augment', type=str2bool, default=False)  #bool: use random flips during training
parser.add_argument('--poisson', type=str2bool, default=False)  #bool: use poisson noise
parser.add_argument('--photons', type=int, default=150)         #poisson noise parameter
parser.add_argument('--background', type=int, default=3)        #poisson noise parameter

#submit a state file to resume training, the state file is the pth generated at each epoch
parser.add_argument('--resume', type=str, default="")           #path to state file

args = parser.parse_args()

##################### Parameter check ##########################
if not os.path.exists(args.outputs_dir):
    os.makedirs(args.outputs_dir)

##################### System setup ##########################
if args.seed is None:
    args.seed = random.randint(1, 10000)
print("Seed: ", args.seed)
torch.manual_seed(args.seed)

cudnn.benchmark = True
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')


##################### Dataset loading ##########################
train_dataset = TrainDataset(args.train_file,
                             augment=args.augment,
                             poisson=args.poisson,
                             seed=args.seed,
                             photons=args.photons, 
                             background=args.background)
train_dataloader = DataLoader(dataset=train_dataset,
                              batch_size=args.batch_size,
                              shuffle=True,
                              num_workers=args.num_workers,
                              pin_memory=True,
                              drop_last=True)
eval_dataset = EvalDataset(args.eval_file,
                           poisson=args.poisson,
                           seed=args.seed,
                           photons=args.photons,
                           background=args.background)
eval_dataloader = DataLoader(dataset=eval_dataset, batch_size=1)


################### Model and optimizer setup ########################
print("Loading model...")
if args.model == "SRCNN":
    model = SRCNN().to(device)

    optimizer = optim.Adam([
        {"params": model.layer1.parameters()},
        {"params": model.layer2.parameters()},
        {"params": model.layer3.parameters(), "lr": args.lr * 0.1}
    ], lr=args.lr)
    
elif args.model == "AE":
    model = AE().to(device)
    optimizer = optim.Adam(model.parameters(), lr=args.lr) 

resume_epoch = 0
best_epoch = 0
best_psnr = 0.0

# If resuming a previous training, load relevant data
if args.resume:
    # load file
    checkpoint = torch.load(args.resume)
        
    model.load_state_dict(checkpoint['model_state_dict'])
    model.to(device)
    
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

    for state in optimizer.state.values():
        for k, v in state.items():
            if torch.is_tensor(v):
                state[k] = v.to(device)
                
    resume_epoch = checkpoint['epoch'] + 1
    best_epoch = checkpoint['best_epoch']
    best_psnr = checkpoint['best_psnr']

# Define loss
criterion = nn.MSELoss().to(device)



####################### Main loop ############################
model = model.float()
for epoch in range(resume_epoch, args.num_epochs):
    
    start_time = time.time()
    print("Starting epoch {}/{}".format(epoch, args.num_epochs - 1));
    
    ##################### Training loop ##########################
    model.train()
    
    train_losses = []
    
    for inputs, targets in train_dataloader:
        # Set model gradients to zero
        optimizer.zero_grad()
        
        inputs = inputs.to(device).float()
        targets = targets.to(device).float()

        preds = model(inputs)
        loss = criterion(preds, targets)

        loss.backward()
        optimizer.step()


        train_losses.append(loss.item())
    
    train_loss = sum(train_losses)/len(train_losses)
    print('Training loss: {}'.format(train_loss))
        
    ##################### Eval loop ##########################
    model.eval()
    
    validation_losses = []
    validation_psnrs  = []

    for inputs, targets in eval_dataloader:

        inputs = inputs.to(device).float()
        targets = targets.to(device).float()

        with torch.no_grad():
            preds = model(inputs).clamp(0.0, 1.0)

        loss = criterion(preds, targets)
        
        eval_psnr = 10 * math.log10((targets.max() ** 2) / loss.item())
        
        validation_losses.append(loss.item())
        validation_psnrs.append(eval_psnr)
    
    validation_loss = np.mean(validation_losses)
    validation_psnr = np.mean(validation_psnrs)
    print('eval psnr: {:.2f}, eval loss: {}'.format(validation_psnr, validation_loss))
    
    if validation_psnr > best_psnr:
        best_epoch = epoch
        best_psnr = validation_psnr
    
    end_time = time.time()
    print('epoch completed in {} seconds'.format(end_time-start_time))
    
    # save epoch's data
    torch.save({
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'training_loss': train_loss,
        'validation_loss': validation_loss,
        'eval_psnr': validation_psnr,
        'best_epoch': best_epoch,
        'best_psnr': best_psnr
        }, os.path.join(args.outputs_dir, 'epoch_{}.pth'.format(epoch)))
    
##################### Training wrap-up ##########################     
print('best epoch: {}, psnr: {:.2f}'.format(best_epoch, best_psnr))

best_checkpoint = torch.load(args.outputs_dir + 'epoch_{}.pth'.format(best_epoch))

torch.save(best_checkpoint, os.path.join(args.outputs_dir, 'best.pth'))
    

