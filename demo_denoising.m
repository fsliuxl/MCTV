%% HSI Experiment
clc; clear;
close all;
seed = 2023;
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
%% Path Information
addpath(genpath('utils'));
addpath(genpath('Datasets'));
addpath(genpath('My code'));
%% Load Data
X0 = double(imread('219090_rgb2gray','png'));
clean_data = Normalize(X0);
gaussian_level = 0.0;
sparse_level   = 0.3;
noise_data       = GetNoise(clean_data,gaussian_level,sparse_level);
%% MCTV
opts.X        = clean_data;% used for calculating the quality indices
opts.lambda   = 4/sqrt(max(size(clean_data)));
opts.mu       = 1e-2;
opts.rho      = 1.25;
opts.MaxIter  = 200;
opts.tol      = 1e-8;
[L, S]    = mctv_rpca(noise_data,opts);

