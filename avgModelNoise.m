%% avgModelNoise.m %%
%
%       By: Josh Wilson
%       Created: June 2022
%
%   Subtracts the mean model (1 set of avg scans) from another scan and plots the residuals.
%   At some point, I will do a model comparison of the types of noise.
%
%   The data2 input is the mean model; data1 is the true time series. You can set r2 cutoffs (and other stuff) in the script.
%
%
%   Usage: avgModelNoise('data1=s0401mc12hrfFit.mat','data2=s0401mc345hrfFit.mat');
 

function [rois cleanRois] = shufflecor(varargin)       

%% Load Data %%

getArgs(varargin);

load(data1); cleanRoisA = rois;
load(data2); cleanRoisB = rois;


%% Subtract mean model and plot noise %%
for roi = 1:length(cleanRoisA)
    for voxel = 1:length(cleanRoisA(roi).vox.r2)
        if (cleanRoisA(roi).vox.r2(voxel) > .5 & cleanRoisB(roi).vox.r2(voxel) > .5)
            figure(roi), hold on, scatter(zscore(cleanRoisA(roi).vox.tSeries(:,voxel)),zscore(cleanRoisA(roi).vox.tSeries(:,voxel))-zscore(cleanRoisB(roi).vox.tSeries(:,voxel)),1,'filled','k');
            xlabel('model time series % signal'),ylabel('residual (mean model - true)');
        end
    end
end





