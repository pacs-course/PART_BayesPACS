clc
clear all
close all

%% USER DEFINED SECTION
% SET DIMENSIONS
d = 1;
N = 5;
M = 10; % 20; %40;

% SET MASTER DIRECTORY
dirMaster = 'C:\Users\Paola\Desktop\PART_BayesPacs-master';

%% END USER DEFINED SECTION

%% Matlab directory
dirMatlab = [dirMaster '\examples\MATLAB code'];

%% Data directory
dirData=[dirMaster '\examples\data\d' num2str(d)];

%% Import true values and fullchain
chdir(dirData);
true = csvread('true.csv');
full = csvread('fullchain.csv');

%% PART's configurations
configNames ={'kdOneNoSmooth.csv' 'kdPairNoSmooth.csv' 'kdOneSmooth.csv' 'kdPairSmooth.csv' 'mlOneNoSmooth.csv' 'mlPairNoSmooth.csv' 'mlOneSmooth.csv' 'mlPairSmooth.csv'};

%% MEASURES

% base name in C
baseNameInC = ['outPART_M', num2str(M), '_N', num2str(N), '_'];

rmse = [];
postConcRatio=[];
for i=1:length(configNames)
   
	% read PART output
    chdir(dirData);
    fileInC = [baseNameInC, configNames{i}];
    aggr_NEW = csvread(fileInC);
        
    % compute rmse, postConcRatio, KL_div
    chdir(dirMatlab);
    rmse = [ rmse; rmse_posterior_mean(aggr_NEW,full)];
    postConcRatio = [ postConcRatio; relative_error(aggr_NEW, full, true)];
	
	clear aggr_NEW;
end

clc
chdir(dirData);

% WRITE ON FILE ACCURACY
diary on
configNames

rmse;
postConcRatio
diary off
