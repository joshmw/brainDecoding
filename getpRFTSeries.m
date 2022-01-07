function [rois cleanRois] = getpRFtSeries(overlayNum,scanNum)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load scan info and set some parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = newView
v = viewSet(v,'curGroup','Concatenation');
nScans = viewGet(v,'nScans');
v = viewSet(v,'curScan',nScans);
v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

% manually define the ROIS you want to look at %
v1 = loadROITSeries(v,'v1',1,1)
v2 = loadROITSeries(v,'v2',1,1)
rois = [v1 v2]



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
   rois(roi).vox.rfStimOverlapTSeries(:,voxel) = m.rfStimCor(:); 
  
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
cleanRois(roi).vox.rfStimOverlapTSeries = rois(roi).vox.rfStimOverlapTSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.pRFtSeries = rois(roi).vox.pRFtSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
end



%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%
graphStuff = 0
if graphStuff

% Eccentricity %
figure(5); hold on;   
scatter(graph.eccentricity,sqrt(graph.measurementVar)*100);
%ylim([0,.0001]);
xlabel('Eccentricity');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by Eccentricity // numvoxels = %i',length(graph.linearCoords)));
p = polyfit(graph.eccentricity,sqrt(graph.measurementVar),1);
plot(linspace(eccMin,eccMax),polyval(p,linspace(eccMin,eccMax)));

% Receptive field size %
figure(6); hold on;
scatter(graph.rfHalfWidth,sqrt(graph.measurementVar)*100);
%ylim([0,.001]);
xlabel('rf HalfWidth');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by rf Size // numvoxels = %i',length(graph.linearCoords)));
p = polyfit(graph.rfHalfWidth,sqrt(graph.measurementVar),1);
plot(linspace(rfWmin,rfWmax),polyval(p,linspace(rfWmin,rfWmax)));

% Model r2 %
figure(7); hold on;   
scatter(graph.r2,sqrt(graph.measurementVar)*100);
%ylim([0,.001]);
xlabel('r2');ylabel('Baseline std (% signal change)');title(sprintf('Baseline std by r2 // numvoxels = %i',sum((eccMax>v1.vox.eccentricity>eccMin)&(v1.vox.r2>r2min))));
p = polyfit(graph.r2,sqrt(graph.measurementVar),1);
plot(linspace(r2min,1),polyval(p,linspace(r2min,1)));
   



% Noise as a function of temporal dynamics %
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





%% noise std in signal absent and present time frames %%
responseQ = .7
figure(20); hold on;
x = []; y = [];

for voxel = 1:length(graph.linearCoords)
    x = [x std(graph.baselineNoise(graph.pRFtSeries(:,voxel)>quantile(graph.pRFtSeries(:,voxel),responseQ),voxel))];
    y = [y std(graph.baselineNoise(graph.pRFtSeries(:,voxel)<quantile(graph.pRFtSeries(:,voxel),1-responseQ),voxel))];
end

scatter(x*100,y*100,'black');
title(sprintf('Signal Present vs Signal Absent Baseline Std: top and bottom %0f percent voxel activity',(1-responseQ)*100))
xlabel('Baseline Std: signal PRESENT time series (% signal change)');
ylabel('Baseline Std: signal ABSENT time series (% signal change)');
xlim([0,2.5]);ylim([0,2.5]);
p = polyfit(x,y,1);
%plot(linspace(0,4),polyval(p,linspace(0,4)),'black');
plot(linspace(0,4),linspace(0,4),'red');




% noise and receptive field correlation %
%%% i need to change this to a parfor loop, currently takes like 3-4 hours for 300 voxels %%%
for row = 1:length(graph.linearCoords)
    disp(sprintf('(correlation) Row %i of %i',row,length(graph.linearCoords)));
    for column = 1:length(graph.linearCoords)
        rfCorVal = corrcoef(mglMakeGaussian(25,25,graph.rfstd(column),graph.rfstd(column),graph.x(column),graph.y(column)),mglMakeGaussian(25,25,graph.rfstd(row),graph.rfstd(row),graph.x(row),graph.y(row)));
        rfCor(row,column) = rfCorVal(1,2);
    end
end
tSeriesBaselineCor = corrcoef(graph.baselineNoise);
figure(151);hold on;
for i = 1:length(graph.linearCoords);
    scatter(rfCor(i,:),tSeriesBaselineCor(i,:));
end
title('Receptive Field and Baseline TimeSeries correlations')
xlabel('Receptive field correlation (voxel i,j)')
ylabel('Baseline time series correlation (voxel i,j)')

rfCorArr = reshape(rfCor,[1 length(graph.linearCoords)^2]);
tSeriesBaselineCorArr = reshape(tSeriesBaselineCor,[1 length(graph.linearCoords)^2]);
[rfCorArr,sortOrder] = sort(rfCorArr);
tSeriesBaselineCorArr = tSeriesBaselineCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -1:.01:.99
    if sum( (bin < rfCorArr) & (rfCorArr < bin+.01) ) > 100
        noiseCorAvg = median(tSeriesBaselineCorArr((bin < rfCorArr) & (rfCorArr < bin+.01)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end
plot(bins,noiseCorAvgs,'black','LineWidth',8);





% stimulus overlap correlation and noise correlation %
rfStimOverlapCor = corrcoef(graph.rfStimOverlapTSeries);
rfStimOverlapCorArr = reshape(rfStimOverlapCor,[1 length(graph.linearCoords)^2]);
rfStimOverlapCorArr = rfStimOverlapCor(sortOrder);

figure(152); hold on;
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signal to noise ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
snrtSeries = graph.baselineNoise./graph.pRFtSeries;
snrtSeriesCor = corrcoef(snrtSeries);
figure(154); hold on;
for i = 1:length(graph.linearCoords);
    scatter(rfStimOverlapCor(i,:),snrtSeriesCor(i,:));
end







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



end

