################################ nano_grab.py ###########################################
# This file is a script that is designed to run on the NVIDIA Jetson Nano. It is a demo
# cript that takes a picture from the camera, and applies the given model to it.
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

from pypylon import pylon
import numpy as np
import torch
import torch.backends.cudnn as cudnn
import argparse
import matplotlib
from matplotlib import pyplot as plt

from AE import AE
from SRCNN import SRCNN

parser = argparse.ArgumentParser()
parser.add_argument('--weights', type=str, required=True) # Path to pth formatted weights file
parser.add_argument('--model',   type=str, required=True) # Model architecture, either "SRCNN" or "AE"
args = parser.parse_args()


cudnn.benchmark = False
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print("using device: " + str(device))

# Load model
if args.model == "SRCNN":
    model = SRCNN().to(device)

else:
    model = AE().to(device)
    
# Load model weights
checkpoint = torch.load(args.weights, map_location=device)
model.load_state_dict(checkpoint['model_state_dict'])

model.eval()
print("model loaded succesfully")

# Setup camera
camera = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
camera.Open()

numberOfImagesToGrab = 1
camera.StartGrabbingMax(numberOfImagesToGrab)

# Grab pictures
while camera.IsGrabbing():
    print("Start grabbing")
    grabResult = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)

    if grabResult.GrabSucceeded():
        print("Grab succeeded")
        
        img = torch.tensor(grabResult.Array, device=device).unsqueeze(0).unsqueeze(0)/255.
        # Resize image dimensions to multiples of 4 for AE pooling layer compatibility
        img = torch.narrow(img, 2, 0, np.shape(img)[2] - np.shape(img)[2]%4)
        img = torch.narrow(img, 3, 0, np.shape(img)[3] - np.shape(img)[3]%4)
        
        print("Running prediction...")
        # Run prediction
        model = model.double()
        with torch.no_grad():
            pred = model(img.double()).clamp(0.0, 1.0)
            
        # Create figure
        f, axarr = plt.subplots(1,2, figsize=(20,10))
        axarr[0].imshow(img.squeeze(0).squeeze(0), cmap='gray')
        axarr[1].imshow(pred.squeeze(0).squeeze(0), cmap='gray')

        axarr[0].title.set_text('Input')
        axarr[1].title.set_text('Prediction')

        f.savefig("./grab.png")
        
        print("Done, exiting...")
        
    grabResult.Release()
camera.Close()
