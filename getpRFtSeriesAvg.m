function [roi pRFtSeries] = getpRFtSeries(v,overlayNum,scanNum,roi)

%% make sure you are on the right scan numbers!! %%

%roi = loadROITSeries(v,'roi',1,1);

v = viewSet(v,'curGroup','Averages');
nScans = viewGet(v,'nScans');
v = viewSet(v,'curScan',nScans);

v = loadAnalysis(v,'pRFAnal/pRF');
tSeries = loadTSeries(v);

v1 = loadROITSeries(v,'v1',1,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the pRF-predicted time series of all voxels in the ROI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for voxel = 1:roi.n

    %% Get info on the voxel %%
    
    % voxel coordinates
    x = roi.scanCoords(1,voxel); y = roi.scanCoords(2,voxel); z = roi.scanCoords(3,voxel);
    
    % grab computed analyses %
    scanNum = scanNum; 
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
   roi.voxAv.linearCoords(voxel) = whichVoxel;
   roi.voxAv.params(:,voxel) = params;
   roi.voxAv.r2(voxel) = thisR2;
   roi.voxAv.polarAngle(voxel) = thisPolarAngle;
   roi.voxAv.eccentricity(voxel) = thisEccentricity;
   roi.voxAv.rfHalfWidth(voxel) = thisRfHalfWidth;
   roi.voxAv.tSeries(:,voxel) = m.tSeries;
   roi.voxAv.pRFtSeries(:,voxel) = m.modelResponse;
   roi.voxAv.baselineNoise = roi.voxAv.tSeries-roi.voxAv.pRFtSeries;
   roi.voxAv.measurementVar(voxel) = var(m.tSeries-m.modelResponse);
  
   
   
   if mod(voxel,50) == 0
       disp(sprintf('(getpRFTSeries) Voxel %i of %i',voxel,roi.n));
   end
   
end


v1.voxAv.tSeries = v1.voxAv.tSeries/100+1;
v1.voxAv.pRFtSeries = v1.voxAv.pRFtSeries/100+1;
v1.voxAv.baselineNoise = v1.voxAv.tSeries-v1.voxAv.pRFtSeries;
v1.voxAv.measurementVar = var(v1.voxAv.baselineNoise);
