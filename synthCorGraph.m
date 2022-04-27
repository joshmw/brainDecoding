function synthCorGraph()

vox = evalin('caller','vox');
param = evalin('caller','param');
rfOverlapRec = evalin('caller','rfOverlapRec');
noiseCor = evalin('caller','noiseCor');


%% Set a noise correlation level for search for%%
corLevel = input(sprintf('Set a minimum correlation to search for (-1 to 1, maximum: %0.2i):',max(max(noiseCor(noiseCor<1)))));

%% Show number of noise cors > corLevel and pick the first voxel %%
sprintf('Number of noise correlations greater than %1.1d, by voxel:',corLevel);

bigCors = find(sum(noiseCor(:,:) > corLevel)>1);

if isempty(bigCors);
    sprintf('!!!!! No 2 voxels have a noise correlation that high !!!!!')
    return; end
bigCors
voxel = input("Enter a voxel you want to look at: ")


%% Show the voxels vox1 has high noise correlations with and pick one to look at %%
sprintf('Voxels that have %1.1d> noise correlation with voxel %i: ', corLevel, voxel)
find(noiseCor(:,voxel) > corLevel)
vox2 = voxel; vox1 = input('Enter a second voxel to look at: ')



%plot time series%
figure;subplot(2,3,1)
s1 = scatter(1:length(vox{1}.trueSeries),vox{vox1}.noisyTrueSeries,'red');hold on;plot(1:length(vox{1}.trueSeries),vox{vox1}.recoveredSeries,'red')
s2 = scatter(1:length(vox{1}.trueSeries),vox{vox2}.noisyTrueSeries,'black');hold on;plot(1:length(vox{1}.trueSeries),vox{vox2}.recoveredSeries,'black')
legend([s1 s2],{sprintf('Voxel %1.0f',vox1),sprintf('Voxel %1.0f',vox2)})
title(sprintf('Voxels %1.0f,%1.0f time series',vox1,vox2)); ylabel('BOLD signal (au)'); xlabel('Timecourse')
subplot(2,3,2)
plot(1:length(vox{1}.trueSeries),vox{vox1}.noiseSeries,'red'); hold on
plot(1:length(vox{1}.trueSeries),vox{vox2}.noiseSeries,'black');
title(sprintf('Voxels %1.0f,%1.0f noise (noise correlation = %.2f)',vox1,vox2,noiseCor(vox1,vox2))); xlabel('Timecourse'); ylabel('Residual BOLD signal (au)')


%show correlations%
subplot(2,3,4)
s3 = scatter(vox{vox1}.noisyTrueSeries,vox{vox2}.noisyTrueSeries,'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox2));
title('Correlation')

subplot(2,3,5)
s3 = scatter(vox{vox1}.noiseSeries,vox{vox2}.noiseSeries,'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',vox2));
title('Correlation')


subplot(4,6,5);imshow(rescale(vox{vox1}.OGrf)); title(sprintf('Voxel %1.0f True RF',vox1))
subplot(4,6,6);imshow(rescale(vox{vox1}.Rrf)); title(sprintf('Voxel %1.0f Recovered RF',vox1))

subplot(4,6,11);imshow(rescale(vox{vox2}.OGrf)); title(sprintf('Voxel %1.0f True RF',vox2))
subplot(4,6,12);imshow(rescale(vox{vox2}.Rrf)); title(sprintf('Voxel %1.0f Recovered RF',vox2)) 
