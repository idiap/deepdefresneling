% ============================ dataset_generator.m ======================================
% This code creates a train, eval, and test dataset by sampling an already existing
% datasets and by fresnel-transforming the samples.
% 
% Parameters in sections "Transform parameters" and "Dataset parameters" are up to you to
% modify for your use case.
%
% 
% This file is based on the following papers:
% M. Liebling, T. Blu, M. Unser, "Fresnelets: new multiresolution wavelet bases for
% digital holography", IEEE Trans. Image Proc., vol. 12, no. 1, pp. 29-43, January 2003.
%
% N. Chacko, M. Liebling, T. Blu, "Discretization of continuous convolution operators for
% accurate modeling of wave propagation in digital holography",
% J. Opt. Soc. Am. A, vol. 30, no10, pp. 2012-2020, 2013. 
%			
%
% This file is distributed under the following license:
% =======================================================================================
% Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
% Written by Rémi Clerc <remi.clerc@idiap.ch>
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of
% conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list
% of conditions and the following disclaimer in the documentation and/or other materials
% provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
% WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% =======================================================================================

clc;
clear all;
close all;
addpathstogcvfftfiles()

%% Transform parameters

% Physical parameters of the Fresnel transform
lambda=522e-9; % Wavelength in meters
d = 8e-3; % distance in meters
T = 10e-6; % sampling distance between pixels
tau=sqrt(lambda*d); % tau parameter, by definition
deg=0;

% Transform structure
s1.h = 'FrT';
s1.h_arg = tau/T;
s1.phi1 = 'cardinal-bspline';
s1.deg1 = deg;
s1.phi2 = 'delta';


%% Dataset parameters
output_path = 'xxx/';

dataset_path = 'path_to/places205/images256/';
train_folder = 'a/abbey/';
eval_folder  = 'a/airport_terminal/';
test_folder  = 'a/alley/';

train_set_size = 2000;
eval_set_size = 10;
test_set_size = 10;

%% Train dataset

% output folders
grey_path = strcat(output_path, 'train/greyscale/');
fres_path = strcat(output_path, 'train/fresnelet/');

% file index and suffix
index = 0;
suffix = '.jpg';

% get file list
path = strcat(dataset_path, train_folder);
listing = dir(path);

% jpg files start from index 4
start_index = 4;
end_index = start_index + train_set_size - 1;

tic; % start timer
for k = start_index:end_index
    
    % load image and convert as greyscale
    I = im2gray(imread(strcat(path, listing(k).name)));
    
    % pad to 512x512
    N = length(I);
    zpad = (512 - N)/2; % to have a power of two-sized image
    padded = [zeros(zpad,2*zpad+N);zeros(N,zpad) I zeros(N,zpad);zeros(zpad,2*zpad+N)];
    
    % fresnelet transform
    fres = stretchcontrast(gcvfft(padded, s1),0,255);
    
    % take module of wave function
    fres = uint8(abs(fres));
    
    % crop back to 256x256
    cropped = fres(zpad+1:zpad+N, zpad+1:zpad+N);
    
    % save image
    imwrite(I,          strcat(grey_path, int2str(index), suffix));
    imwrite(cropped,    strcat(fres_path, int2str(index), suffix));
    index = index + 1;    
    
end
time = toc; % stop timer
fprintf('Train dataset: converted %d images in %f seconds\n', end_index-start_index + 1, time);

%% Eval dataset

% output folders
grey_path = strcat(output_path, 'eval/greyscale/');
fres_path = strcat(output_path, 'eval/fresnelet/');

% file index and suffix
index = 0;
suffix = '.jpg';


% get file list
path = strcat(dataset_path, eval_folder);
listing = dir(path);

% jpg files start from index 4
start_index = 4;
end_index = start_index + eval_set_size - 1;

tic; % start timer
for k = start_index:end_index
    
    % load image and convert as greyscale
    I = im2gray(imread(strcat(path, listing(k).name)));
    
    % pad to 512x512
    N = length(I);
    zpad = (512 - N)/2; % to have a power of two-sized image
    padded = [zeros(zpad,2*zpad+N);zeros(N,zpad) I zeros(N,zpad);zeros(zpad,2*zpad+N)];
    
    % fresnelet transform
    fres = stretchcontrast(gcvfft(padded, s1),0,255);
    
    % take module of wave function
    fres = uint8(abs(fres));
    
    % crop back to 256x256
    cropped = fres(zpad+1:zpad+N, zpad+1:zpad+N);
    
    % save image
    imwrite(I,          strcat(grey_path, int2str(index), suffix));
    imwrite(cropped,    strcat(fres_path, int2str(index), suffix));
    index = index + 1;    
    
end
time = toc; % stop timer
fprintf('Eval dataset: converted %d images in %f seconds\n', end_index-start_index + 1, time);

%% Test dataset

% output folders
grey_path = strcat(output_path, 'test/greyscale/');
fres_path = strcat(output_path, 'test/fresnelet/');

% file index and suffix
index = 0;
suffix = '.jpg';

% get file list
path = strcat(dataset_path, test_folder);
listing = dir(path);

% jpg files start from index 4
start_index = 4;
end_index = start_index + test_set_size - 1;

tic; % start timer
for k = start_index:end_index
    
    % load image and convert as greyscale
    I = im2gray(imread(strcat(path, listing(k).name)));
    
    % pad to 512x512
    N = length(I);
    zpad = (512 - N)/2; % to have a power of two-sized image
    padded = [zeros(zpad,2*zpad+N);zeros(N,zpad) I zeros(N,zpad);zeros(zpad,2*zpad+N)];
    
    % fresnelet transform
    fres = stretchcontrast(gcvfft(padded, s1), 0, 255);
    
    % take module of wave function
    fres = uint8(abs(fres));
    
    % crop back to 256x256
    cropped = fres(zpad+1:zpad+N, zpad+1:zpad+N);
    
    % save image
    imwrite(I,          strcat(grey_path, int2str(index), suffix));
    imwrite(cropped,    strcat(fres_path, int2str(index), suffix));
    index = index + 1;    
    
end
time = toc; % stop timer
fprintf('Test dataset: converted %d images in %f seconds\n', end_index-start_index + 1, time);


