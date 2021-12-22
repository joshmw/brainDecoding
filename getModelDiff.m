function pRFTSeries = getpRFTSeries(v,overlayNum,scanNum,roi)

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
    
    
   %% get the model time series for the voxel **
   m = pRFFit(v,scanNum,x,y,z,'stim',d.stim,'getModelResponse=1','params',params,'concatInfo',d.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo);
   pRFTSeries(x,y,z,:) = m.modelResponse;
   
   if mod(voxel,10)
       disp(sprintf('(getpRFTSeries) Voxel %s of %s, %s% complete',voxel,roi.n,voxel/roi.n);
       disp(sprintf('(pRF) Running on scans %s:%s (restrict %s)',params.groupName,num2str(params.scanNum,'%i '),params.restrict ));
   end
end
k=2

   
   
   
   