function [rois cleanRois] = getpRFtSeries()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load scan info and set some parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = newView
v = viewSet(v,'curGroup','Concatenation');
v = viewSet(v,'curScan',groupScanNum); %remember to set the scan number and analysis you want. I will automate this later
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
parfor voxel = 1:rois(roi).n
    
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
cleanRois(roi).vox.pRFtSeries = rois(roi).vox.pRFtSeries(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
cleanRois(roi).vox.scanCoords = rois(roi).scanCoords(:,(eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth))
end



groupRois = cleanRois;
whichVoxels = (eccMax>rois(roi).vox.eccentricity>eccMin)&(rois(roi).vox.r2>r2min)&(rois(roi).vox.rfHalfWidth>rfWmin)&(rfWmax>rois(roi).vox.rfHalfWidth);


cleanRois = []; rois = [];
deleteView(v)





%%%%%%%%%%%%%%%%% 
% get single params %
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load scan info and set some parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = newView
v = viewSet(v,'curGroup','Concatenation');
v = viewSet(v,'curScan',1); %remember to set the scan number and analysis you want. I will automate this later
v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

% manually define the ROIS you want to look at %
v1 = loadROITSeries(v,'v1')
v2 = loadROITSeries(v,'v2')
v3 = loadROITSeries(v,'v3')

rois = [v1 v2 v3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the pRF-predicted time series of all voxels in the ROI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for roi = 1:length(rois)
parfor voxel = 1:rois(roi).n
    
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
cleanRois(roi).vox.linearCoords = rois(roi).vox.linearCoords(whichVoxels)
cleanRois(roi).vox.rfHalfWidth = rois(roi).vox.rfHalfWidth(whichVoxels)
cleanRois(roi).vox.r2 = rois(roi).vox.r2(whichVoxels)
cleanRois(roi).vox.polarAngle = rois(roi).vox.polarAngle(whichVoxels)
cleanRois(roi).vox.eccentricity = rois(roi).vox.eccentricity(whichVoxels)
cleanRois(roi).vox.tSeries = rois(roi).vox.tSeries(:,whichVoxels)
cleanRois(roi).vox.measurementVar = rois(roi).vox.measurementVar(whichVoxels)
cleanRois(roi).vox.x = rois(roi).vox.params(1,whichVoxels);
cleanRois(roi).vox.y = rois(roi).vox.params(2,whichVoxels);
cleanRois(roi).vox.rfstd = rois(roi).vox.params(3,whichVoxels);
cleanRois(roi).vox.scanCoords = rois(roi).scanCoords(:,whichVoxels)


cleanRois(roi).vox.pRFtSeries = groupRois(roi).vox.tSeries
cleanRois(roi).vox.baselineNoise = cleanRois(roi).vox.tSeries - cleanRois(roi).vox.pRFtSeries
end



%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%

cleanvoxs = 0;
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
    cleanvoxs = cleanvoxs + sum(cleanRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & cleanRois(roi).vox.pRFtSeries(:,voxel)<=highbound);
    end
end

%if roi == 1; vox2show = 80; 
%showvoxclean(cleanRois, groupRois, bin, roi, vox2show, lowbound, highbound); plot(cleanRois(roi).vox.pRFtSeries(:,vox2show)); end;

z = z(~isnan(z)); avgs = avgs(~isnan(avgs));
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

figure(9);
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bin by group prf t series %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


groupvoxs = 0;
for roi = 1:length(cleanRois)
    x = []; y = []; avg = []; errhigh = []; errlow = []; errmeanhigh = []; errmeanlow = []; yAll = []; avgAll = [];
seg = .1; for bin = 0:seg:(1-seg)
    z = []; avgs = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    if mean(groupRois(roi).vox.pRFtSeries(:,voxel)) > .5 %this fixes an issue where some model responses are centered around 0
    lowbound = min(groupRois(roi).vox.pRFtSeries(:,voxel))+bin*(max(groupRois(roi).vox.pRFtSeries(:,voxel))-min(groupRois(roi).vox.pRFtSeries(:,voxel)));    
    highbound = min(groupRois(roi).vox.pRFtSeries(:,voxel))+(bin+seg)*(max(groupRois(roi).vox.pRFtSeries(:,voxel))-min(groupRois(roi).vox.pRFtSeries(:,voxel)));
    z = [z std(cleanRois(roi).vox.baselineNoise(groupRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & groupRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    avgs = [avgs mean(cleanRois(roi).vox.baselineNoise(groupRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & groupRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    groupvoxs = groupvoxs + sum(groupRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & groupRois(roi).vox.pRFtSeries(:,voxel)<=highbound);
    end
end

%if roi == 1; vox2show = 80;
%showvoxclean(cleanRois, groupRois, bin, roi, vox2show, lowbound, highbound); plot(groupRois(roi).vox.pRFtSeries(:,vox2show)); end;

z = z(~isnan(z)); avgs = avgs(~isnan(avgs));
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


figure(10)
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

