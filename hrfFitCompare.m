%% hrfFitCompare.m %%
%       By: Josh Wilson
%       Created: April 2022
%
% Takes the rois/cleanRois of data1 (unfit hrf) and data2 (fit), saved from getprftseries and
% compares parameters. If fitting a 2 parameter hrf (lag and tau), should
% by 5 parameters total you can look at. Just have r2 coded up right now
% but can very easily add others.
%
% Example usage:
%   hrfFitCompare('s0401pRF','s0401hrfFitContrainedLM')
%

function hrfFitCompare(data1,data2)

% load data %
load(data1,'rois','cleanRois')
urois = rois; ucleanRois = cleanRois; %save as u(nfit hrf)rois
clear('rois','cleanRois')
load(data2,'rois','cleanRois')

%hardcoded for v1-v3, can change later if we to be more flexible...
goodVoxels{1} = urois(1).vox.r2 > .6; goodVoxels{2} = urois(2).vox.r2 > .6; goodVoxels{3} = urois(3).vox.r2 > .6; 
roiNames = {'v1','v2','v3'}


for i = 1:3;figure;histogram(cleanRois(i).vox.params(4,:),'BinWidth',.25);xlabel('hrfLag');title(roiNames{i});xlim([-5 5]);end;
for i = 1:3;figure;histogram(cleanRois(i).vox.params(5,:),'BinWidth',.1);xlabel('tau');title(roiNames{i});end;


figure;for i = 1:3; hold on; scatter(urois(i).vox.r2,rois(i).vox.r2);end; plot([0,1],[0,1]);xlabel('Unfit r2 (voxel)');ylabel('Fit r2 (voxel)');legend(roiNames);title('Voxel r2, fit & unfit hrf');