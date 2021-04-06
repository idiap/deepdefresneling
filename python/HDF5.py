####################################### HDF5.py #########################################
# This code takes a dataset folder as input, and creates a HDF5 formatted file containing
# the dataset.
#
# The folder given as input must have the following structure:
# folder/
#   - fresnelet/ (fresnel-transformed samples)
#      - xxx.jpg
#   - greyscale/ (originals or targets)
#      - xxx.jpg
#
# A fresnelet-transformed sample and its target must have the same file name.
#
# The file given as output must have the ".h5" file extension, and is the file to be used
# as input of train.py.
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

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.image as mpimg

import argparse
import os
from os import listdir

import h5py

################## Arguments ##################
parser = argparse.ArgumentParser()
parser.add_argument('--input-folder', type=str, required=True) # path to folder containing the dataset
parser.add_argument('--output-file', type=str, required=True)  # path to output file, with extension ".h5"
args = parser.parse_args()
###############################################

fres_path = args.input_folder + 'fresnelet/'
grey_path = args.input_folder + 'greyscale/'

filenames = listdir(fres_path)

if os.path.exists(args.output_file):
    os.remove(args.output_file) # overwrite
    
with h5py.File(args.output_file, 'a') as hdf5_file:
    #instantiate first image for each array
    file = filenames[0]
    fres_img = np.expand_dims(mpimg.imread(fres_path + file), axis=0)
    grey_img = np.expand_dims(mpimg.imread(grey_path + file), axis=0)

    hr_set = grey_img # high resolution, or greyscale original (name is remnant of SRCNN code)
    lr_set = fres_img # low resolution, or fresenelet transform

    for index in range(1, len(filenames)):
        file = filenames[index]
        fres_img = np.expand_dims(mpimg.imread(fres_path + file), axis=0)
        grey_img = np.expand_dims(mpimg.imread(grey_path + file), axis=0)

        hr_set = np.append(hr_set, grey_img, axis=0)
        lr_set = np.append(lr_set, fres_img, axis=0)

    # save arrays in the file
    hr = hdf5_file.create_dataset('hr', hr_set.shape, data=hr_set, dtype='uint8')
    lr = hdf5_file.create_dataset('lr', lr_set.shape, data=lr_set, dtype='uint8')
