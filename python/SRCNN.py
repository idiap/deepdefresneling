############################ SRCNN.py ##########################################
# This file contains the pytorch implementation of SRCNN, and is distributed
# under the following License:
# ==============================================================================
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
# ==============================================================================

from torch import nn

# SRCNN: C. Dong & al.: Image Super-Resolution Using DeepConvolutional Networks
# https://github.com/Lornatang/SRCNN-PyTorch
class SRCNN(nn.Module):
    def __init__(self, init_weights=True):
        super(SRCNN, self).__init__()

        self.layer1 = nn.Sequential(
            nn.Conv2d(1, 64, kernel_size=9, padding=9 // 2),
            nn.ReLU(inplace=True))

        self.layer2 = nn.Sequential(
            nn.Conv2d(64, 32, kernel_size=5, padding=5 // 2),
            nn.ReLU(inplace=True))

        self.layer3 = nn.Conv2d(32, 1, kernel_size=5, padding=5 // 2)


    def forward(self, x):
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        return x
