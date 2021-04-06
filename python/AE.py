#################################### AE.py ##############################################
# This file contains the pytorch module for the
# Autoencoder architecture
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

from torch import nn  

# AE: homebrewed symmetrical autoencoder architecture
class AE(nn.Module):
    def __init__(self):
        super(AE, self).__init__()
       
        #Encoder
        self.conv1 = nn.Conv2d(1, 16, kernel_size=19, padding=9)  
        self.conv2 = nn.Conv2d(16, 4, kernel_size=15, padding=7)
        self.pool = nn.MaxPool2d(kernel_size=2, stride=2, return_indices=True)
       
        #Decoder
        self.t_conv1 = nn.ConvTranspose2d(4, 16, kernel_size=15, padding=7)
        self.t_conv2 = nn.ConvTranspose2d(16, 1, kernel_size=19, padding=9)
        self.unpool = nn.MaxUnpool2d(kernel_size=2, stride=2)

        self.relu = nn.ReLU(inplace=True)

    def forward(self, x):
        x        = self.relu(self.conv1(x))
        x, pool1 = self.pool(x)
        x        = self.relu(self.conv2(x))
        x, pool2 = self.pool(x)
        
        x        = self.unpool(x, pool2)
        x        = self.relu(self.t_conv1(x))
        x        = self.unpool(x, pool1)
        x        = self.relu(self.t_conv2(x))
              
        return x
