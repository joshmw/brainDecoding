function [rois cleanRois] = getpRFtSeries(overlayNum,scanNum)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load scan info and set some parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = newView
v = viewSet(v,'curGroup','Concatenation');
nScans = viewGet(v,'nScans');
v = viewSet(v,'curScan',1); %remember to set the scan number and analysis you want. I will automate this later
v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

% manually define the ROIS you want to look at %
v1 = loadROITSeries(v,'v1')
v2 = loadROITSeries(v,'v2')
v3 = loadROITSeries(v,'v3')

rois = [v1 v2 v3]

graphStuff = 1; correlate = 1;
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
  
   if mod(voxel,100) == 0
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



%%%%%%%%%%%%%%%%
%% RF Overlap %%
%%%%%%%%%%%%%%%%

% receptive field overlaps %
sprintf('Correlation RFs and calculating distances (takes about a minute)...')

%v1
roi = 1; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)   
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v1rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end

% v2 %
roi = 2; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)      
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v2rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v2rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end


roi = 3; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)       
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v3rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end




% v1/v2 
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(2).vox.linearCoords)       
        mu1 = 0; s1 = cleanRois(1).vox.rfstd(roi1vox); s2 = cleanRois(2).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRois(1).vox.x(roi1vox)-cleanRois(2).vox.x(roi2vox))^2 + (cleanRois(1).vox.y(roi1vox)-cleanRois(2).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1v2rfOverlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v1v2rfOverlap = transpose(v1v2rfOverlap);



% v1/v3
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(3).vox.linearCoords)
        mu1 = 0; s1 = cleanRois(1).vox.rfstd(roi1vox); s2 = cleanRois(3).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRois(1).vox.x(roi1vox)-cleanRois(3).vox.x(roi2vox))^2 + (cleanRois(1).vox.y(roi1vox)-cleanRois(3).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1v3rfOverlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v1v3rfOverlap = transpose(v1v3rfOverlap);


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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%
if graphStuff
    
    
sprintf('Graphing things (also takes about a minute)...')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise std in signal absent and present time frames %%

figure(9);
for roi = 1:length(cleanRois)
    x = []; y = []; avg = []; errhigh = []; errlow = []; errmeanhigh = []; errmeanlow = []; yAll = []; avgAll = [];
seg = .1; for bin = 0:seg:(1-seg)
    z = []; avgs = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    if mean(cleanRois(roi).vox.pRFtSeries(:,voxel)) > .5 %this fixes an issue where some model responses are centered around 0
    lowbound = min(cleanRois(roi).vox.pRFtSeries(:,voxel))+bin*(max(cleanRois(roi).vox.pRFtSeries(:,voxel))-min(cleanRois(roi).vox.pRFtSeries(:,voxel)));    
    highbound = min(cleanRois(roi).vox.pRFtSeries(:,voxel))+(bin+seg)*(max(cleanRois(roi).vox.pRFtSeries(:,voxel))-min(cleanRois(roi).vox.pRFtSeries(:,voxel)));
    z = [z std(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & cleanRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    avgs = [avgs mean(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & cleanRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    end
end

x = [x bin];
y = [y mean(z)];
yAll{length(x)} = z;
avg = [avg mean(avgs)];
avgAll{length(x)} = avgs;

% errorbars %
for i = 1:1000; err(i) = mean(datasample(z,length(cleanRois(roi).vox.linearCoords)));end; err = sort(err);
errhigh = [errhigh quantile(err,.95)-mean(z)];
errlow = [errlow mean(z)-quantile(err,.05)];

for i = 1:1000; err(i) = mean(datasample(avgs,length(cleanRois(roi).vox.linearCoords)));end; err = sort(err);
errmeanhigh = [errmeanhigh quantile(err,.95)-mean(avgs)];
errmeanlow = [errmeanlow mean(avgs)-quantile(err,.05)];
end

% plot %
subplot(2,length(cleanRois),roi); hold on;
%for i = 1:length(x); scatter(repmat(x(i),1,length(cleanRois(roi).vox.linearCoords)),yAll{i}*100); hold on; end;
scatter(x,y*100,60,'black','filled'); 
eb = errorbar(x,y*100,errlow*100,errhigh*100,'o'); eb.Color = 'black';

title(sprintf('%s Residual Std by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Noise (% signal change)');

subplot(2,length(cleanRois),roi+3); hold on;
%for i = 1:length(x); scatter(repmat(x(i),1,length(cleanRois(roi).vox.linearCoords)),avgAll{i}*100); hold on; end;
scatter(x,avg*100,60,'black','filled');plot([0,1],[0,0],'black--');
eb = errorbar(x,avg*100,errmeanlow*100,errmeanhigh*100,'o'); eb.Color = 'black';


title(sprintf('%s Average Residual by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Average Residual (% signal change)');

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field overlaps cutoff by distance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distanceHigh = 10; distanceLow = 3;

% v1 %
figure(30); v1NoiseArrClose = v1NoiseCor(v1dist<distanceLow & v1dist>0); v1NoiseArrFar = v1NoiseCor(v1dist>distanceHigh);
v1kldArrClose = v1rfOverlap(v1dist<distanceLow & v1dist>0); v1kldArrFar = v1rfOverlap(v1dist>distanceHigh);

subplot(1,2,1);scatter(v1kldArrClose,v1NoiseArrClose);hold on;
bins = [];noiseCorAvgs = []; step = .01;
for bin = 0:step:40
    if sum( (bin < v1kldArrClose) & (v1kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v1NoiseArrClose((bin < v1kldArrClose) & (v1kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('V1 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v1kldArrFar,v1NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v1kldArrFar) & (v1kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v1NoiseArrFar((bin < v1kldArrFar) & (v1kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('V1 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


% v2 %
figure(31); v2NoiseArrClose = v2NoiseCor(v2dist<distanceLow & v2dist>0); v2NoiseArrFar = v2NoiseCor(v2dist>distanceHigh);
v2kldArrClose = v2rfOverlap(v2dist<distanceLow & v2dist>0); v2kldArrFar = v2rfOverlap(v2dist>distanceHigh);

subplot(1,2,1);scatter(v2kldArrClose,v2NoiseArrClose);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldArrClose) & (v2kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v2NoiseArrClose((bin < v2kldArrClose) & (v2kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v2 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v2kldArrFar,v2NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldArrFar) & (v2kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v2NoiseArrFar((bin < v2kldArrFar) & (v2kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v2 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


% v3 %
figure(32); v3NoiseArrClose = v3NoiseCor(v3dist<distanceLow & v3dist>0); v3NoiseArrFar = v3NoiseCor(v3dist>distanceHigh);
v3kldArrClose = v3rfOverlap(v3dist<distanceLow & v3dist>0); v3kldArrFar = v3rfOverlap(v3dist>distanceHigh);

subplot(1,2,1);scatter(v3kldArrClose,v3NoiseArrClose);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldArrClose) & (v3kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v3NoiseArrClose((bin < v3kldArrClose) & (v3kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v3 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v3kldArrFar,v3NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldArrFar) & (v3kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v3NoiseArrFar((bin < v3kldArrFar) & (v3kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v3 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);subplot(1,2,1);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1rfOverlap(i,:),v1NoiseCor(i,:)); end;

v1rfOverlapArr = reshape(v1rfOverlap,[1 length(v1rfOverlap)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1rfOverlapArr,sortOrder] = sort(v1rfOverlapArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v1rfOverlapArr) & (v1rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1rfOverlapArr) & (v1rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V1 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j (mm)');

% second, do v2 %
figure(12);hold on;subplot(1,2,1);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2rfOverlap(i,:),v2NoiseCor(i,:)); end;

v2rfOverlapArr = reshape(v2rfOverlap,[1 length(v2rfOverlap)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2rfOverlapArr,sortOrder] = sort(v2rfOverlapArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2rfOverlapArr) & (v2rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2rfOverlapArr) & (v2rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V2 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


% third, do v3 %
figure(13);hold on;subplot(1,2,1);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3rfOverlap(i,:),v3NoiseCor(i,:)); end;

v3rfOverlapArr = reshape(v3rfOverlap,[1 length(v3rfOverlap)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3rfOverlapArr,sortOrder] = sort(v3rfOverlapArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3rfOverlapArr) & (v3rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3rfOverlapArr) & (v3rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V3 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


%between v1 and v2
figure(14);hold on;subplot(1,2,1);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2rfOverlap(i,:),v1v2NoiseCor(i,:)); end;

v1v2rfOverlapArr = reshape(v1v2rfOverlap,[1 min(size(v1v2rfOverlap))*max(size(v1v2rfOverlap))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2rfOverlapArr,sortOrder] = sort(v1v2rfOverlapArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2rfOverlapArr) & (v1v2rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2rfOverlapArr) & (v1v2rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,1]);
title('Receptive Field and Noise Correlations between v1 and v2'); xlabel('KL Divergence between voxel V1i, V2j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V2j (mm)');


%between v1 and v3
figure(15);hold on;subplot(1,2,1);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3rfOverlap(i,:),v1v3NoiseCor(i,:)); end;

v1v3rfOverlapArr = reshape(v1v3rfOverlap,[1 min(size(v1v3rfOverlap))*max(size(v1v3rfOverlap))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3rfOverlapArr,sortOrder] = sort(v1v3rfOverlapArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3rfOverlapArr) & (v1v3rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3rfOverlapArr) & (v1v3rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,1]);
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

bins = [];noiseCorAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 15
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
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 15
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
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 15
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
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 15
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
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 15
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

for i = 1:length(v1rfOverlap); scatter(v1dist(i,:),v1rfOverlap(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1rfOverlapArr = reshape(v1rfOverlap,[1 length(v1rfOverlap)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1rfOverlapArr = v1rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1rfOverlapArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V1 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v2 %
figure(25);hold on;subplot(2,3,2);hold on;

for i = 1:length(v2rfOverlap); scatter(v2dist(i,:),v2rfOverlap(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2rfOverlapArr = reshape(v2rfOverlap,[1 length(v2rfOverlap)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2rfOverlapArr = v2rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 15
        rfOverlapAvg = median(v2rfOverlapArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V2 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v3 %
figure(25);hold on;subplot(2,3,3);hold on;

for i = 1:length(v3rfOverlap); scatter(v3dist(i,:),v3rfOverlap(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3rfOverlapArr = reshape(v3rfOverlap,[1 length(v3rfOverlap)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3rfOverlapArr = v3rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = []; 
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 15
        rfOverlapAvg = median(v3rfOverlapArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V3 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


%between v1 and v2
figure(25);hold on;subplot(2,3,4);hold on;

for i = 1:min(size(v1v2rfOverlap)); scatter(v1v2dist(i,:),v1v2rfOverlap(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2rfOverlapArr = reshape(v1v2rfOverlap,[1 min(size(v1v2rfOverlap))*max(size(v1v2rfOverlap))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2rfOverlapArr = v1v2rfOverlapArr(sortOrder);

bins = [];rfOverlapAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1v2rfOverlapArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end

plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('Distance and RF Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('KL divergence between voxel RFs V1i, V2j (distance)');


%between v1 and v3
figure(25);hold on;subplot(2,3,5);hold on;

for i = 1:min(size(v1v3rfOverlap)); scatter(v1v3dist(i,:),v1v3rfOverlap(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3rfOverlapArr = reshape(v1v3rfOverlap,[1 min(size(v1v3rfOverlap))*max(size(v1v3rfOverlap))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3rfOverlapArr = v1v3rfOverlapArr(sortOrder);

bins = [];rfOverlapAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1v3rfOverlapArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end

plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('Distance and RF Correlations between v1 and v3'); xlabel('Distance between voxels V1i, V3j'); ylabel('KL divergence between voxel RFs RFs V1i, V2j (distance)');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time series and noise correlation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% V1
figure(26);hold on;subplot(1,3,1);hold on;

v1tSeriesCor = corrcoef(cleanRois(1).vox.pRFtSeries);

for i = 1:length(v1NoiseCor); scatter(v1tSeriesCor(i,:),v1NoiseCor(i,:)); end;

v1tSeriesCorArr = reshape(v1tSeriesCor,[1 length(v1tSeriesCor)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1tSeriesCorArr,sortOrder] = sort(v1tSeriesCorArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .05;
for bin = -.5:step:1
    if sum( (bin < v1tSeriesCorArr) & (v1tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1tSeriesCorArr) & (v1tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V1 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');


%% V2
figure(26);hold on;subplot(1,3,2);hold on;

v2tSeriesCor = corrcoef(cleanRois(2).vox.pRFtSeries);

for i = 1:length(v2NoiseCor); scatter(v2tSeriesCor(i,:),v2NoiseCor(i,:)); end;

v2tSeriesCorArr = reshape(v2tSeriesCor,[1 length(v2tSeriesCor)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2tSeriesCorArr,sortOrder] = sort(v2tSeriesCorArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; 
for bin = -.5:step:1
    if sum( (bin < v2tSeriesCorArr) & (v2tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2tSeriesCorArr) & (v2tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V2 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');


%% V3
figure(26);hold on;subplot(1,3,3);hold on;

v3tSeriesCor = corrcoef(cleanRois(3).vox.pRFtSeries);

for i = 1:length(v3NoiseCor); scatter(v3tSeriesCor(i,:),v3NoiseCor(i,:)); end;

v3tSeriesCorArr = reshape(v3tSeriesCor,[1 length(v3tSeriesCor)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3tSeriesCorArr,sortOrder] = sort(v3tSeriesCorArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -.5:step:1
    if sum( (bin < v3tSeriesCorArr) & (v3tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3tSeriesCorArr) & (v3tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V3 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');








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

roi = 1; voxels1 = (eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth);
roi = 2; voxels2 = (eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth);
roi = 3; voxels3 = (eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth);
for roi = 1:3; avgRois(roi) = cleanRois(roi);end;


cleanRois(1).vox.tSeries = rois(1).vox.tSeries(:,voxels1)/100+1;
cleanRois(1).vox.pRFtSeries = avgRois(1).vox.tSeries;
cleanRois(1).vox.baselineNoise = cleanRois(1).vox.tSeries-cleanRois(1).vox.pRFtSeries;

cleanRois(2).vox.tSeries = rois(2).vox.tSeries(:,voxels2)/100+1;
cleanRois(2).vox.pRFtSeries = avgRois(2).vox.tSeries;
cleanRois(2).vox.baselineNoise = cleanRois(2).vox.tSeries-cleanRois(2).vox.pRFtSeries;

cleanRois(3).vox.tSeries = rois(3).vox.tSeries(:,voxels3)/100+1;
cleanRois(3).vox.pRFtSeries = avgRois(3).vox.tSeries;
cleanRois(3).vox.baselineNoise = cleanRois(3).vox.tSeries-cleanRois(3).vox.pRFtSeries;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% what is anticorrelation %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%these voxels have anticorrelated noise

s 
sum((v1NoiseCor(:,:) < -.2) & (v1tSeriesCor(:,:) > .9))
find((v1NoiseCor(vox,:) < -.2) & (v1tSeriesCor(vox,:) > .9))





badvox1 = 24;badvox2 = 8;

%plot time series%
figure(100);close;figure(100);subplot(2,2,1)
s1 = scatter(1:240,cleanRois(1).vox.tSeries(:,badvox1),'red');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,badvox1),'red')
s2 = scatter(1:240,cleanRois(1).vox.tSeries(:,badvox2),'black');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,badvox2),'black')
legend([s1 s2],{sprintf('Voxel %1.0f',badvox1),sprintf('Voxel %1.0f',badvox2)})
title(sprintf('Voxels %1.0f,%1.0f time series (correlation = %.2f)',badvox1,badvox2,v1tSeriesCor(badvox1,badvox2))); ylabel('BOLD signal (change from average)'); xlabel('Timecourse')
subplot(2,2,2)
plot(1:240,cleanRois(1).vox.baselineNoise(:,badvox1),'red'); hold on
plot(1:240,cleanRois(1).vox.baselineNoise(:,badvox2),'black');
title(sprintf('Voxels %1.0f,%1.0f noise (noise correlation = %.2f)',badvox1,badvox2,v1NoiseCor(badvox1,badvox2))); xlabel('Timecourse'); ylabel('Residual BOLD signal (change from average)')

%show correlations%
subplot(2,2,3)
s3 = scatter(cleanRois(1).vox.tSeries(:,badvox1),cleanRois(1).vox.tSeries(:,badvox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',badvox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',badvox2));
title('Correlation')

subplot(2,2,4)
s3 = scatter(cleanRois(1).vox.baselineNoise(:,badvox1),cleanRois(1).vox.baselineNoise(:,badvox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',badvox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',badvox2));
title('Correlation')




%plot time series of good voxels%
goodvox1 = 50;goodvox2 = 56;

figure(101);close;figure(101)
subplot(2,2,1)
s1 = scatter(1:240,cleanRois(1).vox.tSeries(:,goodvox1),'red');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,goodvox1),'red')
s2 = scatter(1:240,cleanRois(1).vox.tSeries(:,goodvox2),'black');hold on;plot(1:240,cleanRois(1).vox.pRFtSeries(:,goodvox2),'black')
legend([s1 s2],{sprintf('Voxel %1.0f',goodvox1),sprintf('Voxel %1.0f',goodvox2)})
title(sprintf('Voxels %1.0f,%1.0f time series (correlation = %.2f)',goodvox1,goodvox2,v1tSeriesCor(goodvox1,goodvox2))); ylabel('BOLD signal (change from average)'); xlabel('Timecourse')
subplot(2,2,2)
plot(1:240,cleanRois(1).vox.baselineNoise(:,goodvox1),'red'); hold on
plot(1:240,cleanRois(1).vox.baselineNoise(:,goodvox2),'black');
title(sprintf('Voxels %1.0f,%1.0f noise (noise correlation = %.2f)',goodvox1,goodvox2,v1NoiseCor(goodvox1,goodvox2))); xlabel('Timecourse'); ylabel('Residual BOLD signal (change from average)')


subplot(2,2,3)
s3 = scatter(cleanRois(1).vox.tSeries(:,goodvox1),cleanRois(1).vox.tSeries(:,goodvox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',goodvox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',goodvox2));
title('Correlation')

subplot(2,2,4)
s3 = scatter(cleanRois(1).vox.baselineNoise(:,goodvox1),cleanRois(1).vox.baselineNoise(:,goodvox2),'black');
xlabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',goodvox1));ylabel(sprintf('Voxel %1.0f Activity (BOLD signal (%)',goodvox2));
title('Correlation')





%%%%%%%%%% figure 9 with equal bins
figure(9);
for roi = 1:length(cleanRois)
    x = []; y = []; avg = []; errhigh = []; errlow = [];
for bin = 0:.1:.9
    z = []; avgs = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    if mean(cleanRois(roi).vox.pRFtSeries(:,voxel)) > .5 %this fixes an issue where some model responses are centered around 0
    z = [z std(cleanRois(roi).vox.baselineNoise(quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin+.025)>cleanRois(roi).vox.pRFtSeries(:,voxel)&cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin),voxel))]; 
    avgs = [avgs mean(cleanRois(roi).vox.baselineNoise(quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin+.025)>cleanRois(roi).vox.pRFtSeries(:,voxel)&cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin),voxel))];
    end
end
x = [x bin];
y = [y mean(z)];
errhigh = [errhigh quantile(sort(z),.95)];
errlow = [errlow quantile(sort(z),.05)];
avg = [avg mean(avgs)];
end

subplot(2,length(cleanRois),roi); hold on;
scatter(x,y*100,'black');
%errorbar(x,y*100,errlow*100,errhigh*100);
title(sprintf('%s Residual Std by Activity Quantile',cleanRois(roi).name))

xlabel('Quantile of Voxel Activity');
ylabel('Noise (% signal change)');

subplot(2,length(cleanRois),roi+3); hold on;
scatter(x,avg*100,'black');plot([0,1],[0,0],'black--');
title(sprintf('%s Average Residual by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Average Residual (% signal change)');

end





%end dontDo
end






