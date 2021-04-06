################################# transformations.py ####################################
# This file contains the transformations that may be applied to samples during training
# and evaluation.
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


from torchvision import transforms, datasets
import torchvision.transforms.functional as ttf
import torch
import random
import numpy.random as nprnd
import numpy as np


# Random horizontal and vertical flip
class RndVHFlip:
    def __init__(self, pv=0.5, ph=0.5):
        self.pv = pv
        self.ph = ph

    # call method:
    # There is a 50% of each flip direction to occur. The same flips are performed on
    # both the input image and the target image.
    #
    # inputs:
    #	- images: tuple containing the input image and the target image
    #
    # returns:
    #	- transformed: tuple containing the new input image and the new target image
    def __call__(self, images):
        transformed = images

        if random.random() < self.pv:
            transformed = np.flip(transformed[0], 1), np.flip(transformed[1], 1)

        if random.random() < self.ph:
            transformed = np.flip(transformed[0], 2), np.flip(transformed[1], 2)

        return transformed



# poisson noise. Source: fmarelli
def poisson_noise(image, photons=150, background=3, seed=None, in_place=False):
    """
    Simulate the acquisition noise of the camera:
        noise = shot noise

    Parameters
    ----------
    image : array [ZXY]
        the clean image
    photons : int, optional
        the max level of photons per pixel, by default 150
    background : int, optional
        the background level of photons per pixel, by default 3
    seed : int, optional
        the seed for rng, by default None

    Returns
    -------
    array [ZXY]
        the noisy image
    """
    
    if not in_place:
        image = image.copy()

    image_norm = image * (photons - background) / np.max(image)
    image_norm += background
    rng = nprnd.default_rng(seed)
    poisson = rng.poisson(image_norm)
    
    poisson = poisson/np.max(poisson)

    return poisson
