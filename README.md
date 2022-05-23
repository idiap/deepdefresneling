# DeepDefresneling: Deep Learning Methods for Digital Holography in an Embedded System

This repository contains the code that was used and developed during my Master's thesis entitled Deep Learning Methods for Digital Holography in an Embedded System. The goal of this thesis was to use Deep Learning to reconstruct diffracted samples that were acquired in an in-line digital holography setup. Two approaches were considered: Super Resolution Convolutional Neural Networks (SRCNN), and Denoising Convolutional Autoencoders (AE). These solutions were developed with lightweightedness in mind, so that the models can easily be run by the NVIDIA Jetson Nano micro-computer.

## Getting Started

The files you will use will depend on your own application:

* Do you only want to deploy a pre-trained model on your system? -> You will need the weights files in the [pth folder](./pth)
* Do you want to retrain the models with new parameters? -> You will need matlab and the [dataset generator](./matlab/dataset_generator.m)

In either case, you will need python to run the models.


### Prerequisites

#### Python environment
You need to install [python 3.7](https://www.python.org/downloads/).

It is recommended that you create a virtual environment specifically for this project. You can run the following commands in your terminal to get the environment up and running.

```
virtualenv -p /usr/bin/python3.7 env
source env/bin/activate
pip install h5py jupyter notebook matplotlib numpy pillow torch torchvision
```

You can then clone this repository to your device
```
git clone https://github.com/idiap/deepdefresneling .
```

Each time you come back to work on this environment, you will need to activate the environment by running
```
source env/bin/activate
```

#### Matlab environment
There are no prerequisites for running the provided matlab code

## Using the code

### Generating a new dataset
Any existing image dataset can be used to create a new dataset to train the models. During development, we used the dataset [Places205](http://places.csail.mit.edu/downloadData.html).

Matlab and the [dataset generator](./matlab/dataset_generator.m) can then be used to generate the new dataset. You will need to modify the following parameters:

Transform parameters
```
lambda=522e-9; % Wavelength in meters
d = 8e-3; % distance in meters
T = 10e-6; % sampling distance between pixels
```

The paths to the different dataset folders, and your output path
```
output_path = 'xxx/';

dataset_path = 'path_to/places205/images256/';
train_folder = 'a/abbey/';
eval_folder  = 'a/airport_terminal/':
test_folder  = 'a/alley/';
```

The sizes of the new dataset's train, eval and test sections
```
train_set_size = 2000;
eval_set_size = 10;
test_set_size = 10;

```

The python code is made to deal with a dataset compressed to the HDF5 format. You can easily compress the dataset by running the python script [HDF5.py](./python/HDF5.py) in the following way:
**Note that the output-file argument needs you to put .h5 as file extension, you will need to run this code for the train set and the eval set, but it is not necessary for the test set**
```
python HDF5.py --input-folder path/to/dataset/folder --output-file path/to/output/file.h5
```

### Training a model
If you have a HDF5-formatted dataset, you can train AE or SRCNN by using the python script [train.py](./python/train.py) in the following way:
```
python train.py [arguments]
```

With the following required arguments:
```
--train-file    #path to h5 formatted train file
--eval-file     #path to h5 formatted eval file
--outputs-dir   #output directory
--model         #model architecture, either "SRCNN" or "AE"
```

And the following optional arguments:
```
--lr            #learning rate, default 1e-4
--batch-size    #default 16
--num-epochs	#number of training epochs, default 400
--num-workers   #default 8
--seed          #random seed
--augment       #bool: use random flips during training, default False
--poisson       #bool: use poisson noise, default False
--photons       #poisson noise parameter, default 150
--background    #poisson noise parameter, default 3
```

If your training was interrupted, or you want to resume a training at any point, you can use the following optional argument:
```
--resume    #path to state file from which to resume, default ""
```
The state files are generated as pth files after each epoch in the output directory

When training ends, a weights file named **best.pth** will be generated in the output directory, it will contain the state of the epoch that achieved the best evaluation PSNR.

### Testing a model
If you have a pth-formatted weights file and a test dataset, you can test your model with the provided jupyter notebook [test.ipynb](./python/test.ipynb). It contains the relevant instructions for you to run the tests.

## Deployment
The electronics of the physical setup consist of a Basler ace acA2000-165uc camera and a NVIDIA Jetson Nano micro-computer that communicate with each other via USB. The setup of the electronics is to be done on the NVIDIA Jetson Nano.

### Getting started with the NVIDIA Jetson Nano
**First, follow [this tutorial](https://developer.nvidia.com/embedded/learn/get-started-jetson-nano-devkit) to setup the OS on the Jetson Nano.**

Python 3.6 is already installed on the provided OS image.

### Installing Pytorch and Torchvision
Since the Jetson Nano's architecture is a bit different from conventional computers, the installation process for Pytorch and Torchvision may seem a little bit convoluted, as you cannot simply use the standard packages.

First, install pip and setup a folder where everything will stay
```
sudo apt-get install python3-pip
mkdir python
cd python
```

Then, install torch 1.7.0 by getting the compatible wheel from NVIDIA:
```
wget https://nvidia.box.com/shared/static/cs3xn3td6sfgtene6jdvsxlr366m2dhq.whl -O torch-1.7.0-cp36-cp36m-linux_aarch64.whl
sudo apt-get install libopenblas-base libopenmpi-dev 
pip3 install Cython
pip3 install numpy torch-1.7.0-cp36-cp36m-linux_aarch64.whl
```

Then, install torchvision 0.8.1:
```
sudo apt-get install libjpeg-dev zlib1g-dev libpython3-dev libavcodec-dev libavformat-dev libswscale-dev
git clone --branch v0.8.1 https://github.com/pytorch/vision torchvision
cd torchvision
export BUILD_VERSION=0.8.1
python3 setup.py install --user
cd ..
```

You can test that the libraries are correctly installed by running:
```
python3
>>> import torch
>>> torch.__version__         # should be 1.7.0
>>> torch.cuda.is_available() # should be True
>>>
>>> import torchvision
>>> torchvision.__version__   # should be 0.8.1
>>> quit()
```

### Installing Pylon and Pypylon
Pylon is the software used to run the Basler camera, and Pypylon is the python library that allows to use the camera with python scripts, both are necessary to get the system up and running.

Download the [Pylon 6.1.3 debian package](https://www.baslerweb.com/en/sales-support/downloads/software-downloads/pylon-6-1-3-linux-arm-64-bit-debian/) in your Jetson nano and run it by simply double-clicking it in your file explorer. This will allow you to install Pylon.

You will need to install Pypylon from source:
```
sudo apt install gcc python-dev swig
git clone https://github.com/basler/pypylon.git pypylon
cd pypylon
pip3 install .
```

When it is done, you should be able to take pictures with the Basler camera from the Jetson Nano. As a test, you can connect the Basler camera to the Jetson Nano via USB, and run:
```
cd pypylon/samples
python3 save_image.py
```

This will save 5 pictures to the current directory, have a look!

### Acquiring samples and running them through the model
Clone this repository, connect the camera to the Jetson Nano, and you can run the test script!
```
git clone https://github.com/idiap/deepdefresneling ddf
python3 ddf/python/nano_grab.py --weights path/to/pth/file.pth --model [AE/SRCNN]
```

## Author
* **Nara Clerc** - [Minauras](https://github.com/Minauras)

## Acknowledgments

* [**Liu Changyu**](https://github.com/Lornatang) - [SRCNN implementation in pytorch](https://github.com/Lornatang/SRCNN-PyTorch)
* **Nikhil Chacko, Michael Liebling** - [Fresnelet Transform implementation in matlab](http://sybil.ece.ucsb.edu/pages/gcvfft/index.html)

This work was supported by the Swiss National Science Foundation (SNSF), Grants:
* 206021_164022 'Platform for Reproducible Acquisition, Processing, and Sharing of Dynamic, Multi-Modal Data' (PLATFORM_MMD)
* 200020_179217 'Computational biomicroscopy: advanced image processing methods to quantify live biological systems' (COMPBIO)

## License
This is free software: you can redistribute it and/or modify it under the terms of the BSD-2-Clause licence.
