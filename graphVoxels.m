function graphVoxels(area)

%% Set a noise correlation level for search for%%
corLevel = input("Set a minimum correlation to search for (-1 to 1):");


%% Load variables from the environment. Need to pre-compute correlations and stuff (speed issue) %%
areaNoiseCorStr = strcat(area,'NoiseCor'); areaNoiseCor = evalin('base',areaNoiseCorStr);
areaTseriesCorStr = strcat(area,'tSeriesCor'); areaTseriesCor = evalin('base',areaTseriesCorStr);
cleanRois = evalin('base','cleanRois'); rois = evalin('base','rois');
% note - load cleanRois and run the correlation code in base environment from getpRFTSeries.m


%% Show number of noise cors > corLevel and pick the first voxel %%
sprintf('Number of noise correlations greater than %1.1d, by voxel:',corLevel);
bigCors = find(sum(areaNoiseCor(:,:) > corLevel)>1);
if isempty(bigCors);
    sprintf('!!!!! No 2 voxels have a noise correlation that high !!!!!')
    return; end
bigCors
vox = input("Enter a voxel you want to look at: ")


%% Show the voxels vox1 has high noise correlations with and pick one to look at %%
sprintf('Voxels that have %1.1d> noise correlation with voxel %i: ', corLevel, vox)
find(areaNoiseCor(vox,:) > corLevel)
vox1 = vox; vox2 = input('Enter a second voxel to look at: ')


%plot time series%
figure(100);close;figure(100);subplot(2,2,1)
s1 = scatter(1:240,cleanRois(1).vox.tSeries(:,vox1),'red');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,vox1),'red')
s2 = scatter(1:240,cleanRois(1).vox.tSeries(:,vox2),'black');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,vox2),'black')
legend([s1 s2],{sprintf('Voxel %1.0f',vox1),sprintf('Voxel %1.0f',vox2)})
title(sprintf('Voxels %1.0f,%1.0f time series (correlation = %.2f)',vox1,vox2,areaTseriesCor(vox1,vox2))); ylabel('BOLD signal (change from average)'); xlabel('Timecourse')
subplot(2,2,2)
plot(1:240,cleanRois(1).vox.baselineNoise(:,vox1),'red'); hold on
plot(1:240,cleanRois(1).vox.baselineNoise(:,vox2),'black');
title(sprintf('Voxels %1.0f,%1.0f noise (noise correlation = %.2f)',vox1,vox2,areaNoiseCor(vox1,vox2))); xlabel('Timecourse'); ylabel('Residual BOLD signal (change from average)')


%show correlations%
subplot(2,2,3)
s3 = scatter(cleanRois(1).vox.tSeries(:,vox1),cleanRois(1).vox.tSeries(:,vox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox2));
title('Correlation')

subplot(2,2,4)
s3 = scatter(cleanRois(1).vox.baselineNoise(:,vox1),cleanRois(1).vox.baselineNoise(:,vox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox2));
title('Correlation')