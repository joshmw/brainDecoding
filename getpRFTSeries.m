function [rois cleanRois] = getpRFtSeries(overlayNum,scanNum)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load scan info and set some parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = newView
v = viewSet(v,'curGroup','MotionComp');
nScans = viewGet(v,'nScans');
v = viewSet(v,'curScan',1); %remember to set the scan number and analysis you want. I will automate this later
v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

% manually define the ROIS you want to look at %
v1 = loadROITSeries(v,'v1',1,1)
v2 = loadROITSeries(v,'v2',1,1)
v3 = loadROITSeries(v,'v3',1,1)

rois = [v1 v2 v3]

graphStuff = 0; correlate = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the pRF-predicted time series of all voxels in the ROI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for roi = 1:length(rois)
for voxel = 1:rois(roi).n
    
    % grab computed analyses %
    x = rois(roi).scanCoords(1,voxel); y = rois(roi).scanCoords(2,voxel); z = rois(roi).scanCoords(3,voxel);
    a = viewGet(v,'Analysis');
    d = viewGet(v,'d',scanNum);
    r2 = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','r2'));
    thisR2 = r2(x,y,z);
    polarAngle = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','polarAngle'));
    thisPolarAngle = polarAngle(x,y,z);
    eccentricity = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','eccentricity'));
    thisEccentricity = eccentricity(x,y,z);
    rfHalfWidth = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','rfHalfWidth'));
    thisRfHalfWidth = rfHalfWidth(x,y,z);
    
    % get params  %
    scanDims = viewGet(v,'scanDims',scanNum);
    whichVoxel = find(d.linearCoords == sub2ind(scanDims,x,y,z));
    r = d.r(whichVoxel,:);
    params = d.params(:,whichVoxel);
    
    if isfield(d,'paramsInfo')
        paramsInfo = d.paramsInfo;
    else
        paramsInfo = [];
    end 
    
   % get the model time series for the voxel and parameters %
   m = pRFFit(v,scanNum,x,y,z,'stim',d.stim,'getModelResponse=1','params',params,'concatInfo',d.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo);
   rois(roi).vox.linearCoords(voxel) = whichVoxel;
   rois(roi).vox.params(:,voxel) = params;
   rois(roi).vox.r2(voxel) = thisR2;
   rois(roi).vox.polarAngle(voxel) = thisPolarAngle;
   rois(roi).vox.eccentricity(voxel) = thisEccentricity;
   rois(roi).vox.rfHalfWidth(voxel) = thisRfHalfWidth;
   rois(roi).vox.tSeries(:,voxel) = m.tSeries;
   rois(roi).vox.pRFtSeries(:,voxel) = m.modelResponse;
   rois(roi).vox.baselineNoise = rois(roi).vox.tSeries-rois(roi).vox.pRFtSeries;
   rois(roi).vox.measurementVar(voxel) = var(m.tSeries-m.modelResponse);
   %rois(roi).vox.rfStimOverlapTSeries(:,voxel) = m.rfStimCor(:); 
  
   if mod(voxel,50) == 0
       disp(sprintf('(getpRFTSeries) Roi %i of %i; Voxel %i of %i',roi,length(rois),voxel,rois(roi).n));
   end
   
end
end



%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Good Voxels %%
%%%%%%%%%%%%%%%%%%%%%%%%

% set the cutoffs for voxels you want to look at %
r2min = .6
eccMin = .1; eccMax = 25;
rfWmin = 1; rfWmax = 20;
cleanRois = rois

% filter for voxels that meet cutoffs %
for roi = 1:length(rois)
cleanRois(roi).vox.linearCoords = rois(roi).vox.linearCoords((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.rfHalfWidth = rois(roi).vox.rfHalfWidth((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.r2 = rois(roi).vox.r2((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.polarAngle = rois(roi).vox.polarAngle((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.eccentricity = rois(roi).vox.eccentricity((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.tSeries = rois(roi).vox.tSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.measurementVar = rois(roi).vox.measurementVar((eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.baselineNoise = rois(roi).vox.baselineNoise(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.x = rois(roi).vox.params(1,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth));
cleanRois(roi).vox.y = rois(roi).vox.params(2,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth));
cleanRois(roi).vox.rfstd = rois(roi).vox.params(3,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth));
%cleanRois(roi).vox.rfStimOverlapTSeries = rois(roi).vox.rfStimOverlapTSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.pRFtSeries = rois(roi).vox.pRFtSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.scanCoords = rois(roi).scanCoords(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
end



%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlate voxel RFs %%
%%%%%%%%%%%%%%%%%%%%%%%%%
if correlate
% trying KL divergence of multivariate gaussians %

%v1 rf self correlation
roi = 1; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        kld = mvgkl([cleanRois(roi).vox.x(column);cleanRois(roi).vox.y(column)], [cleanRois(roi).vox.rfstd(column)^2 1; 1 cleanRois(roi).vox.rfstd(column)^2], [cleanRois(roi).vox.x(row);cleanRois(roi).vox.y(row)], [cleanRois(roi).vox.rfstd(row)^2 1; 1 cleanRois(roi).vox.rfstd(row)^2]);
        v1kldCor(row,column) = kld;
    end
end

%v2 rf self correlation
roi = 2; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        kld = mvgkl([cleanRois(roi).vox.x(column);cleanRois(roi).vox.y(column)], [cleanRois(roi).vox.rfstd(column)^2 1; 1 cleanRois(roi).vox.rfstd(column)^2], [cleanRois(roi).vox.x(row);cleanRois(roi).vox.y(row)], [cleanRois(roi).vox.rfstd(row)^2 1; 1 cleanRois(roi).vox.rfstd(row)^2]);
        v2kldCor(row,column) = kld;
    end
end

%v3 rf self correlation
roi = 3; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        kld = mvgkl([cleanRois(roi).vox.x(column);cleanRois(roi).vox.y(column)], [cleanRois(roi).vox.rfstd(column)^2 1; 1 cleanRois(roi).vox.rfstd(column)^2], [cleanRois(roi).vox.x(row);cleanRois(roi).vox.y(row)], [cleanRois(roi).vox.rfstd(row)^2 1; 1 cleanRois(roi).vox.rfstd(row)^2]);
        v3kldCor(row,column) = kld;
    end
end


% v1/v2 
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    %disp(sprintf('(correlation) Voxel %i of %i',roi1vox,length(cleanRois(1).vox.linearCoords)));
    for roi2vox = 1:length(cleanRois(2).vox.linearCoords)
        kld = mvgkl([cleanRois(1).vox.x(roi1vox);cleanRois(1).vox.y(roi1vox)], [cleanRois(1).vox.rfstd(roi1vox)^2 1; 1 cleanRois(1).vox.rfstd(roi1vox)^2], [cleanRois(2).vox.x(roi2vox);cleanRois(2).vox.y(roi2vox)], [cleanRois(2).vox.rfstd(roi2vox)^2 1; 1 cleanRois(2).vox.rfstd(roi2vox)^2]);
        v1v2kldCor(roi2vox,roi1vox) = kld;
    end
end


% v1/v3
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    %disp(sprintf('(correlation) Voxel %i of %i',roi1vox,length(cleanRois(1).vox.linearCoords)));
    for roi2vox = 1:length(cleanRois(3).vox.linearCoords)
        kld = mvgkl([cleanRois(1).vox.x(roi1vox);cleanRois(1).vox.y(roi1vox)], [cleanRois(1).vox.rfstd(roi1vox)^2 1; 1 cleanRois(1).vox.rfstd(roi1vox)^2], [cleanRois(3).vox.x(roi2vox);cleanRois(3).vox.y(roi2vox)], [cleanRois(3).vox.rfstd(roi2vox)^2 1; 1 cleanRois(3).vox.rfstd(roi2vox)^2]);
        v1v3kldCor(roi2vox,roi1vox) = kld;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%% voxel distances %%
%%%%%%%%%%%%%%%%%%%%%

% v1 %
roi = 1; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        dist = sqrt((cleanRois(roi).vox.scanCoords(1,column)-cleanRois(roi).vox.scanCoords(1,row))^2 + (cleanRois(roi).vox.scanCoords(2,column)-cleanRois(roi).vox.scanCoords(2,row))^2 + (cleanRois(roi).vox.scanCoords(3,column)-cleanRois(roi).vox.scanCoords(3,row))^2);
        v1dist(row,column) = dist;
    end
end

% v2 %
roi = 2; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        dist = sqrt((cleanRois(roi).vox.scanCoords(1,column)-cleanRois(roi).vox.scanCoords(1,row))^2 + (cleanRois(roi).vox.scanCoords(2,column)-cleanRois(roi).vox.scanCoords(2,row))^2 + (cleanRois(roi).vox.scanCoords(3,column)-cleanRois(roi).vox.scanCoords(3,row))^2);
        v2dist(row,column) = dist;
    end
end

% v3 %
roi = 3; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        dist = sqrt((cleanRois(roi).vox.scanCoords(1,column)-cleanRois(roi).vox.scanCoords(1,row))^2 + (cleanRois(roi).vox.scanCoords(2,column)-cleanRois(roi).vox.scanCoords(2,row))^2 + (cleanRois(roi).vox.scanCoords(3,column)-cleanRois(roi).vox.scanCoords(3,row))^2);
        v3dist(row,column) = dist;
    end
end

% v1/v2 
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(2).vox.linearCoords)
        dist = sqrt((cleanRois(1).vox.scanCoords(1,roi1vox)-cleanRois(2).vox.scanCoords(1,roi2vox))^2 + (cleanRois(1).vox.scanCoords(2,roi1vox)-cleanRois(2).vox.scanCoords(2,roi2vox))^2 + (cleanRois(1).vox.scanCoords(3,roi1vox)-cleanRois(2).vox.scanCoords(3,roi2vox))^2);
        v1v2dist(roi2vox,roi1vox) = dist;
    end
end


% v1/v3
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(3).vox.linearCoords)
        dist = sqrt((cleanRois(1).vox.scanCoords(1,roi1vox)-cleanRois(3).vox.scanCoords(1,roi2vox))^2 + (cleanRois(1).vox.scanCoords(2,roi1vox)-cleanRois(3).vox.scanCoords(2,roi2vox))^2 + (cleanRois(1).vox.scanCoords(3,roi1vox)-cleanRois(3).vox.scanCoords(3,roi2vox))^2);
        v1v3dist(roi2vox,roi1vox) = dist;
    end
end

end

%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%
if graphStuff
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise std in signal absent and present time frames %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bin the signal change %
figure(9);
for roi = 1:length(cleanRois)
    x = []; y = [];
for bin = 0:.025:.95
    z = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    z = [z std(cleanRois(roi).vox.baselineNoise(quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin+.025)>cleanRois(roi).vox.pRFtSeries(:,voxel)&cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin),voxel))]; 
end
x = [x bin];
y = [y mean(z)];
end
    
subplot(1,length(cleanRois),roi); hold on;
scatter(x,y*100,'black');
title(sprintf('%s Noise by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Noise (% signal change)');
%ylim([0,1]);
%p = polyfit(x,y,1);
%plot(linspace(0,4),polyval(p,linspace(0,4)),'black');
%plot(linspace(0,4),linspace(0,4),'red');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);subplot(1,2,1);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1kldCor(i,:),v1NoiseCor(i,:)); end;

v1kldCorArr = reshape(v1kldCor,[1 length(v1kldCor)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1kldCorArr,sortOrder] = sort(v1kldCorArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .1
for bin = 0:step:40
    if sum( (bin < v1kldCorArr) & (v1kldCorArr < bin+step) ) > 20
        noiseCorAvg = median(v1NoiseCorArr((bin < v1kldCorArr) & (v1kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V1 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxel i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');

% second, do v2 %
figure(12);hold on;subplot(1,2,1);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2kldCor(i,:),v2NoiseCor(i,:)); end;

v2kldCorArr = reshape(v2kldCor,[1 length(v2kldCor)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2kldCorArr,sortOrder] = sort(v2kldCorArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldCorArr) & (v2kldCorArr < bin+step) ) > 20
        noiseCorAvg = median(v2NoiseCorArr((bin < v2kldCorArr) & (v2kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V2 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


% third, do v3 %
figure(13);hold on;subplot(1,2,1);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3kldCor(i,:),v3NoiseCor(i,:)); end;

v3kldCorArr = reshape(v3kldCor,[1 length(v3kldCor)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3kldCorArr,sortOrder] = sort(v3kldCorArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldCorArr) & (v3kldCorArr < bin+step) ) > 20
        noiseCorAvg = median(v3NoiseCorArr((bin < v3kldCorArr) & (v3kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V2 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


%between v1 and v2
figure(14);hold on;subplot(1,2,1);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2kldCor(i,:),v1v2NoiseCor(i,:)); end;

v1v2kldCorArr = reshape(v1v2kldCor,[1 min(size(v1v2kldCor))*max(size(v1v2kldCor))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2kldCorArr,sortOrder] = sort(v1v2kldCorArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2kldCorArr) & (v1v2kldCorArr < bin+step) ) > 20
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2kldCorArr) & (v1v2kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,20]);set(gca,'XDir','reverse');
title('Receptive Field and Noise Correlations between v1 and v2'); xlabel('KL Divergence between voxel V1i, V2j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V2j (mm)');


%between v1 and v3
figure(15);hold on;subplot(1,2,1);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3kldCor(i,:),v1v3NoiseCor(i,:)); end;

v1v3kldCorArr = reshape(v1v3kldCor,[1 min(size(v1v3kldCor))*max(size(v1v3kldCor))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3kldCorArr,sortOrder] = sort(v1v3kldCorArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3kldCorArr) & (v1v3kldCorArr < bin+step) ) > 20
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3kldCorArr) & (v1v3kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,20]);set(gca,'XDir','reverse');
title('Receptive Field and Noise Correlations between v1 and v3'); xlabel('KL Divergence between voxel V1i, V3j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V3j (mm)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distance and noise correlation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);hold on;subplot(1,2,2);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1dist(i,:),v1NoiseCor(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .25
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 20
        noiseCorAvg = median(v1NoiseCorArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);set(gca,'XDir','reverse');
title('V1 Distance and Noise Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');

% second, do v2 %
figure(12);hold on;subplot(1,2,2);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2dist(i,:),v2NoiseCor(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 20
        noiseCorAvg = median(v2NoiseCorArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);set(gca,'XDir','reverse');
title('V2 Distance and Noise Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');


% third, do v3 %
figure(13);hold on;subplot(1,2,2);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3dist(i,:),v3NoiseCor(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 20
        noiseCorAvg = median(v3NoiseCorArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('v3 Distance and Noise correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');


%between v1 and v2
figure(14);hold on;subplot(1,2,2);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2dist(i,:),v1v2NoiseCor(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 20
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('Distance and Noise Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('Noise correlation between voxels V1i, V2j');


%between v1 and v3
figure(15);hold on;subplot(1,2,2);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3dist(i,:),v1v3NoiseCor(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 20
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('Distance and Noise Correlations between v1 and v3'); xlabel('Distance between voxels V1i, V3j'); ylabel('Noise correlation between voxels V1i, V2j');




%%%%%%%%%%%%%%%%%%%%%%%%
%% distance vs rf cor %%
%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(25);hold on;subplot(2,3,1);hold on;

for i = 1:length(v1kldCor); scatter(v1dist(i,:),v1kldCor(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1kldCorArr = reshape(v1kldCor,[1 length(v1kldCor)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1kldCorArr = v1kldCorArr(sortOrder);
bins = [];kldCorAvgs = []; step = .25
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 20
        kldCorAvg = median(v1kldCorArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V1 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v2 %
figure(25);hold on;subplot(2,3,2);hold on;

for i = 1:length(v2kldCor); scatter(v2dist(i,:),v2kldCor(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2kldCorArr = reshape(v2kldCor,[1 length(v2kldCor)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2kldCorArr = v2kldCorArr(sortOrder);
bins = [];kldCorAvgs = []; step = .25
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 20
        kldCorAvg = median(v2kldCorArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V2 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v3 %
figure(25);hold on;subplot(2,3,3);hold on;

for i = 1:length(v3kldCor); scatter(v3dist(i,:),v3kldCor(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3kldCorArr = reshape(v3kldCor,[1 length(v3kldCor)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3kldCorArr = v3kldCorArr(sortOrder);
bins = [];kldCorAvgs = []; step = .25
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 20
        kldCorAvg = median(v3kldCorArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V3 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


%between v1 and v2
figure(25);hold on;subplot(2,3,4);hold on;

for i = 1:min(size(v1v2kldCor)); scatter(v1v2dist(i,:),v1v2kldCor(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2kldCorArr = reshape(v1v2kldCor,[1 min(size(v1v2kldCor))*max(size(v1v2kldCor))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2kldCorArr = v1v2kldCorArr(sortOrder);

bins = [];kldCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 20
        kldCorAvg = median(v1v2kldCorArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end

plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('Distance and RF Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('KL divergence between voxel RFs V1i, V2j (distance)');


%between v1 and v3
figure(25);hold on;subplot(2,3,5);hold on;

for i = 1:min(size(v1v3kldCor)); scatter(v1v3dist(i,:),v1v3kldCor(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3kldCorArr = reshape(v1v3kldCor,[1 min(size(v1v3kldCor))*max(size(v1v3kldCor))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3kldCorArr = v1v3kldCorArr(sortOrder);

bins = [];kldCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 20
        kldCorAvg = median(v1v3kldCorArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end

plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('Distance and RF Correlations between v1 and v3'); xlabel('Distance between voxels V1i, V3j'); ylabel('KL divergence between voxel RFs RFs V1i, V2j (distance)');




end  %%%%% end of graphing stuff














dontdo = 0;if dontdo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% old stuff i don't think is useful %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signal to noise ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
snrtSeries = graph.baselineNoise./graph.pRFtSeries;
snrtSeriesCor = corrcoef(snrtSeries);
figure(154); hold on;
for i = 1:length(graph.linearCoords);
    scatter(rfStimOverlapCor(i,:),snrtSeriesCor(i,:));
end



   
%%%%%%%%%%%%%%%%    
% Eccentricity %
%%%%%%%%%%%%%%%%
figure(5); hold on;   
scatter(graph.eccentricity,sqrt(graph.measurementVar)*100);
%ylim([0,.0001]);
xlabel('Eccentricity');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by Eccentricity // numvoxels = %i',length(graph.linearCoords)));
p = polyfit(graph.eccentricity,sqrt(graph.measurementVar),1);
plot(linspace(eccMin,eccMax),polyval(p,linspace(eccMin,eccMax)));


%%%%%%%%%%%%%%%%%%%%%%%%
% Receptive field size %
%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); hold on;
scatter(graph.rfHalfWidth,sqrt(graph.measurementVar)*100);
%ylim([0,.001]);
xlabel('rf HalfWidth');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by rf Size // numvoxels = %i',length(graph.linearCoords)));
p = polyfit(graph.rfHalfWidth,sqrt(graph.measurementVar),1);
plot(linspace(rfWmin,rfWmax),polyval(p,linspace(rfWmin,rfWmax)));


%%%%%%%%%%%%
% Model r2 %
%%%%%%%%%%%%
figure(7); hold on;   
scatter(graph.r2,sqrt(graph.measurementVar)*100);
%ylim([0,.001]);
xlabel('r2');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by r2 // numvoxels = %i',sum((eccMax>v1.vox.eccentricity>eccMin)&(v1.vox.r2>r2min))));
p = polyfit(graph.r2,sqrt(graph.measurementVar),1);
plot(linspace(r2min,1),polyval(p,linspace(r2min,1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise as a function of temporal dynamics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
derivativeQuantile = .75
figure(8); hold on;
x = []; y = [];
for voxel = 1:roi.n
    x = [x std(v1.vox.baselineNoise(abs(diff(v1.vox.pRFtSeries(:,voxel)))>quantile(abs(diff(v1.vox.pRFtSeries(:,voxel))),derivativeQuantile),voxel))];
    y = [y std(v1.vox.baselineNoise(abs(diff(v1.vox.pRFtSeries(:,voxel)))<quantile(abs(diff(v1.vox.pRFtSeries(:,voxel))),1-derivativeQuantile),voxel))];
end
scatter(x,y,'black');
p = polyfit(x,y,1);
plot(linspace(0,.04),polyval(p,linspace(0,.04)),'black');
plot(linspace(0,.04),linspace(0,.04),'red');
xlim([0,.04]);ylim([0,.04]);
title('Individual voxel baseline noise in high- vs low-derivative time series areas')
xlabel('Std of baseline noise: rapid response (% signal change)');
ylabel('Std of baseline noise: slow response shift (% signal change)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise reduction by averaging %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(9)
x = linspace(0,.0025); y = x*sqrt(2);
plot(x,y,'red'); hold on
scatter(v1.voxAv.measurementVar,v1.voxM.measurementVar,'black')
xlim([0,.0025]); ylim([0,.005]);
title('Average and non-averaged baseline Variance per voxel')
xlabel('Baseline Variance of voxel N in time series averaged across 2 scans (% signal change)')
ylabel('Baseline Variance of voxel N in single scan time series (% signal change)')

figure(10)
x = linspace(0,.0025); y = x*sqrt(5)/sqrt(2);
plot(x,y,'red'); hold on
scatter(v1.vox.measurementVar,v1.voxAv.measurementVar,'black')
xlim([0,.00075]); ylim([0,.0025]);
title('Average and non-averaged baseline Variance per voxel')
xlabel('Baseline Variance of voxel N in time series averaged across 5 scans (% signal change)')
ylabel('Baseline Variance of voxel N in time series averaged across 2 scans (% signal change)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus overlap correlation and noise correlation %     % think i commented out the calculations of this from the getmodelresponse part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rfStimOverlapCor = corrcoef(graph.rfStimOverlapTSeries);
rfStimOverlapCorArr = reshape(rfStimOverlapCor,[1 length(graph.linearCoords)^2]);
rfStimOverlapCorArr = rfStimOverlapCor(sortOrder);

figure(14); hold on;
for i = 1:length(graph.linearCoords);
    scatter(rfStimOverlapCor(i,:),tSeriesBaselineCor(i,:));
end

bins = [];overlapCorAvgs = [];
for bin = -1:.01:.99
    if sum( (bin < rfStimOverlapCorArr) & (rfStimOverlapCorArr < bin+.01) ) > 100
        overlapCorAvg = median(tSeriesBaselineCorArr((bin < rfStimOverlapCorArr) & (rfStimOverlapCorArr < bin+.01)));
        bins = [bins bin]; overlapCorAvgs = [overlapCorAvgs overlapCorAvg];
    end
end
plot(bins,overlapCorAvgs,'black','LineWidth',7);

title('Stimulus overlap time series and Baseline TimeSeries correlations')
xlabel('Stimulus overlap correlation (voxel i,j)')
ylabel('Baseline time series correlation (voxel i,j)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rf and stimulus overlap correlation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is a sanity check

figure(153); hold on;
for i = 1:length(graph.linearCoords);
    scatter(rfCor(i,:),rfStimOverlapCor(i,:));
end
title('Receptive Field and Stimulus Overlap Time Series Correlations')
xlabel('Receptive Field Correlation (voxel i,j)')
ylabel('Stimulus Overlap Correlation (voxel i,j)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% old correlation graphs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I am drawing compressed gaussians in mgl then calculating the correlation between the images. i think overlap/kl divergence is better and muuuch faster
    
% using correlate %    
%correlate receptive fields between areas
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    disp(sprintf('(correlation) Voxel %i of %i',roi1vox,length(cleanRois(1).vox.linearCoords)));
    for roi2vox = 1:length(cleanRois(2).vox.linearCoords)
        rfCorVal = corrcoef(mglMakeGaussian(5,5,cleanRois(1).vox.rfstd(roi1vox)/5,cleanRois(1).vox.rfstd(roi1vox)/5,cleanRois(1).vox.x(roi1vox)/5,cleanRois(1).vox.y(roi1vox)/5),mglMakeGaussian(5,5,cleanRois(2).vox.rfstd(roi2vox)/5,cleanRois(2).vox.rfstd(roi2vox)/5,cleanRois(2).vox.x(roi2vox)/5,cleanRois(2).vox.y(roi2vox)/5));
        areaCor(roi2vox,roi1vox) = rfCorVal(1,2);
    end
end

%v1 rf self correlation
roi = 1; for row = 1:length(cleanRois(roi).vox.linearCoords)
    disp(sprintf('(correlation) Row %i of %i',row,length(cleanRois(roi).vox.linearCoords)));
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        rfCorVal = corrcoef(mglMakeGaussian(5,5,cleanRois(roi).vox.rfstd(column)/5,cleanRois(roi).vox.rfstd(column)/5,cleanRois(roi).vox.x(column)/5,cleanRois(roi).vox.y(column)/5),mglMakeGaussian(5,5,cleanRois(roi).vox.rfstd(row)/5,cleanRois(roi).vox.rfstd(row)/5,cleanRois(roi).vox.x(row)/5,cleanRois(roi).vox.y(row)/5));
        v1SelfCor(row,column) = rfCorVal(1,2);
    end
end

%v2 rf self correlation
roi = 2; for row = 1:length(cleanRois(roi).vox.linearCoords)
    disp(sprintf('(correlation) Row %i of %i',row,length(cleanRois(roi).vox.linearCoords)));
    for column = 1:length(cleanRois(roi).vox.linearCoords)
        rfCorVal = corrcoef(mglMakeGaussian(5,5,cleanRois(roi).vox.rfstd(column)/5,cleanRois(roi).vox.rfstd(column)/5,cleanRois(roi).vox.x(column)/5,cleanRois(roi).vox.y(column)/5),mglMakeGaussian(5,5,cleanRois(roi).vox.rfstd(row)/5,cleanRois(roi).vox.rfstd(row)/5,cleanRois(roi).vox.x(row)/5,cleanRois(roi).vox.y(row)/5));
        v2SelfCor(row,column) = rfCorVal(1,2);
    end
end

% first, do v1 %
figure(11);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1SelfCor(i,:),v1NoiseCor(i,:)); end;

v1SelfCorArr = reshape(v1SelfCor,[1 length(v1SelfCor)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1SelfCorArr,sortOrder] = sort(v1SelfCorArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -1:.01:.99
    if sum( (bin < v1SelfCorArr) & (v1SelfCorArr < bin+.01) ) > 100
        noiseCorAvg = median(v1NoiseCorArr((bin < v1SelfCorArr) & (v1SelfCorArr < bin+.01)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);
title('V1 Receptive Field and Baseline TimeSeries correlations'); xlabel('Receptive field correlation (voxel i,j)'); ylabel('Baseline time series correlation (voxel i,j)');

% second, do v2 %
figure(12);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2SelfCor(i,:),v2NoiseCor(i,:)); end;

v2SelfCorArr = reshape(v2SelfCor,[1 length(v2SelfCor)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2SelfCorArr,sortOrder] = sort(v2SelfCorArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -1:.01:.99
    if sum( (bin < v2SelfCorArr) & (v2SelfCorArr < bin+.01) ) > 100
        noiseCorAvg = median(v2NoiseCorArr((bin < v2SelfCorArr) & (v2SelfCorArr < bin+.01)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);
title('V2 Receptive Field and Baseline TimeSeries correlations'); xlabel('Receptive field correlation (voxel i,j)'); ylabel('Baseline time series correlation (voxel i,j)');


% noise and rf correlation between areas %
figure(13);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(areaCor(i,:),v1v2NoiseCor(i,:)); end;

areaCorArr = reshape(areaCor,[1 min(size(areaCor))*max(size(areaCor))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[areaCorArr,sortOrder] = sort(areaCorArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -1:.01:.99
    if sum( (bin < areaCorArr) & (areaCorArr < bin+.01) ) > 100
        noiseCorAvg = median(v1v2NoiseCorArr((bin < areaCorArr) & (areaCorArr < bin+.01)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);
title('Receptive Field and Noise Correlations between v1 and v2'); xlabel('Receptive field correlation (voxel V1i, V2j)'); ylabel('Baseline time series correlation (voxel V1i, V2j)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise std in signal absent and present time frames %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
responseQ = .85
figure(10);

for roi = 1:length(cleanRois)

    x = []; y = []
    
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    x = [x std(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),responseQ),voxel))];
    y = [y std(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)<quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),1-responseQ),voxel))];
end

subplot(1,2,roi); hold on;
scatter(x*100,y*100,'black');
title(sprintf('%s Baseline Std: top and bottom %1.1i percent voxel activity',cleanRois(roi).name,(1-responseQ)*100))
xlabel('Baseline Std: signal PRESENT time series (% signal change)');
ylabel('Baseline Std: signal ABSENT time series (% signal change)');
xlim([0,2.5]);ylim([0,2.5]);
p = polyfit(x,y,1);
plot(linspace(0,4),polyval(p,linspace(0,4)),'black');
plot(linspace(0,4),linspace(0,4),'red');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean data for leave1out noise modeling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%run average/concatenation and save data as avgRois, then run motioncomp or raw

cleanRois(1).vox.tSeries = rois(1).vox.tSeries(:,voxels1)/100+1;
cleanRois(1).vox.pRFtSeries = avgRois(1).vox.tSeries;
cleanRois(1).vox.baselineNoise = cleanRois(1).vox.tSeries-cleanRois(1).vox.pRFtSeries;

cleanRois(2).vox.tSeries = rois(2).vox.tSeries(:,voxels2)/100+1;
cleanRois(2).vox.pRFtSeries = avgRois(2).vox.tSeries;
cleanRois(2).vox.baselineNoise = cleanRois(2).vox.tSeries-cleanRois(2).vox.pRFtSeries;

cleanRois(3).vox.tSeries = rois(3).vox.tSeries(:,voxels3)/100+1;
cleanRois(3).vox.pRFtSeries = avgRois(3).vox.tSeries;
cleanRois(3).vox.baselineNoise = cleanRois(3).vox.tSeries-cleanRois(3).vox.pRFtSeries;





end






