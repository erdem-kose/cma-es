clear; clc; close all;
warning('off','all');

addpath(genpath('library'));

marker_size=25;
marker_width=2;
font_size=18;

%% Create CMA object
cmaObj=cmaClass;
cmaObj.ChrLength=3; % Chromosome length or dimension length
cmaObj.Lambda=10;
cmaObj.Mu=10;
cmaObj.FunEvl=100000;
cmaObj.Sigma=0.1;
cmaObj.UpLim=[10.00 10.00 10.00]'; 
cmaObj.LwLim=[-10.00 -10.00 -10.00]';

%% Optimize Subtractive Clustering
tic;
cmaObj=cmaObj.optimize(@hart3);%costfun(x) where x[dimensionindex,sampleindex]
toc;
x=cmaObj.solBest;

%% Get Optimized Clusters
hart3(x);
