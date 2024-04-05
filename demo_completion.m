clc; clear;
close all;
seed = 2023;
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
%% Path Information
addpath(genpath('Datasets'));
addpath(genpath('quality assess'));
addpath(genpath('My code'));
%% Input
X0 = double(imread('3096_rgb2gray','png'));
X = Normalize(X0); % Normalization
[n1, n2] = size(X);
SamRate = 0.3;
%% Sampling
omega = find(rand(n1*n2,1)<SamRate);
y = zeros(n1,n2);
y(omega) = X(omega);
%% modified CTV
opts.X        = X;% used for calculating the quality indices
opts.mu       = 1e-2;
opts.rho      = 1.25;
opts.MaxIter  = 200;
opts.tol      = 1e-8;
xrec_mctv     = mctv_mc(y,omega,opts);


