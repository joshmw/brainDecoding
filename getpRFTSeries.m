function [roi pRFtSeries] = getpRFtSeries(v,overlayNum,scanNum,roi)

graphStuff = 0
%roi = loadROITSeries(v,'roi',1,1);

v = viewSet(v,'curGroup','Concatenation');
nScans = viewGet(v,'nScans');
v = viewSet(v,'curScan',nScans);

v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

v1 = loadROITSeries(v,'v1',1,1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the pRF-predicted time series of all voxels in the ROI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for voxel = 1:roi.n

    %% Get info on the voxel %%
    
    % voxel coordinates
    x = roi.scanCoords(1,voxel); y = roi.scanCoords(2,voxel); z = roi.scanCoords(3,voxel);
    
    % grab computed analyses %
    scanNum = 1;  %% need to set this properly
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
    
    %% get params  %%
    scanDims = viewGet(v,'scanDims',scanNum);
    whichVoxel = find(d.linearCoords == sub2ind(scanDims,x,y,z));
    r = d.r(whichVoxel,:);
    params = d.params(:,whichVoxel);
    
    if isfield(d,'paramsInfo')
        paramsInfo = d.paramsInfo;
    else
        paramsInfo = [];
    end 
    
   %% get the model time series for the voxel and parameters **
   m = pRFFit(v,scanNum,x,y,z,'stim',d.stim,'getModelResponse=1','params',params,'concatInfo',d.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo);
   pRFtSeries(x,y,z,:) = m.modelResponse;
   roi.vox.linearCoords(voxel) = whichVoxel;
   roi.vox.params(:,voxel) = params;
   roi.vox.r2(voxel) = thisR2;
   roi.vox.polarAngle(voxel) = thisPolarAngle;
   roi.vox.eccentricity(voxel) = thisEccentricity;
   roi.vox.rfHalfWidth(voxel) = thisRfHalfWidth;
   roi.vox.tSeries(:,voxel) = m.tSeries;
   roi.vox.pRFtSeries(:,voxel) = m.modelResponse;
   roi.vox.baselineNoise = roi.vox.tSeries-roi.vox.pRFtSeries;
   roi.vox.measurementVar(voxel) = var(m.tSeries-m.modelResponse);
  
   
   
   if mod(voxel,50) == 0
       disp(sprintf('(getpRFTSeries) Voxel %i of %i',voxel,roi.n));
   end
   
end


if graphStuff

%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%

% set the cutoffs for voxels you want to look at %
r2min = .6
eccMin = .1; eccMax = 25;
rfWmin = 1; rfWmax = 20;

graph.linearCoords = roi.vox.linearCoords((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.rfHalfWidth = roi.vox.rfHalfWidth((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.r2 = roi.vox.r2((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.polarAngle = roi.vox.polarAngle((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.eccentricity = roi.vox.eccentricity((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.tSeries = roi.vox.tSeries(:,(eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.measurementVar = roi.vox.measurementVar((eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.baselineNoise = roi.vox.baselineNoise(:,(eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth))
graph.x = roi.vox.params(1,(eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth));
graph.y = roi.vox.params(2,(eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth));
graph.rfstd = roi.vox.params(3,(eccMax>roi.vox.eccentricity>eccMin)&(roi.vox.r2>r2min)&(roi.vox.rfHalfWidth>rfWmin)&(rfWmax>roi.vox.rfHalfWidth));

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





% noise std in signal absent and present time frames %
responseQ = .9
figure(20); hold on;
x = []; y = [];
for voxel = 1:length(graph.linearCoords)
    x = [x std(graph.baselineNoise(graph.tSeries(:,voxel)>quantile(graph.tSeries(:,voxel),responseQ),voxel))];
    y = [y std(graph.baselineNoise(graph.tSeries(:,voxel)<quantile(graph.tSeries(:,voxel),1-responseQ),voxel))];
end
scatter(x*100,y*100,'blue');
title('Signal Present vs Signal Absent Baseline Std')
xlabel('Baseline Std: signal PRESENT time series (% signal change)');
ylabel('Baseline Std: signal ABSENT time series (% signal change)');
xlim([0,2.5]);ylim([0,2.5]);
p = polyfit(x,y,1);
plot(linspace(0,4),polyval(p,linspace(0,4)),'blue');
plot(linspace(0,4),linspace(0,4),'red');





% noise and receptive field correlation %
tSeriesCor = corrcoef(graph.tSeries);
%%% i need to change this to a parfor loop, currently takes like 10 hours %%%
for row = 1:length(graph.linearCoords)
    for column = 1:length(graph.linearCoords)
        rfCorVal = corrcoef(mglMakeGaussian(25,25,graph.rfstd(column),graph.rfstd(column),graph.x(column),graph.y(column)),mglMakeGaussian(25,25,graph.rfstd(row),graph.rfstd(row),graph.x(row),graph.y(row)));
        rfCor(row,column) = rfCorVal(1,2);
    end
end
figure(151);hold on;
for i = 1:50;
    scatter(rfCor(i,:),tSeriesBaselineCor(i,:));
end
title('Receptive Field and Baseline TimeSeries correlations')
xlabel('Receptive field correlation (voxel i,j)')
ylabel('Baseline time series correlation (voxel i,j)')




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

