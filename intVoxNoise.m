function intVoxNoise(v,overlayNum,scanNum,x,y,z,roi)

v = getMLRView; scanNum = v.curScan;

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
m.baselineNoise = m.tSeries - m.modelResponse;


%% graph time series %%
figure

subplot(2,2,1)
scatter(1:length(m.tSeries),m.tSeries,'black'); hold on;
plot(1:length(m.modelResponse),m.modelResponse,'black')

title(sprintf('[%i %i %i] r^2=%0.2f eccentricity=%0.2f rfHalfWidth=%0.2f %s',x,y,z,thisR2,thisEccentricity,thisRfHalfWidth,a.params.pRFFit.rfType));


subplot(2,2,3)
plot(1:length(m.tSeries),m.baselineNoise,'black'); hold on;


%% calc noise dynamics %%


    yAll = []; avgsAll = []; x = []; y= []; avgs = []; meds = [];
    
    seg = .1; for bin = 0:seg:(1-seg)
    
    lowbound = min(m.modelResponse)+bin*(max(m.modelResponse)-min(m.modelResponse));    
    highbound = min(m.modelResponse)+(bin+seg)*(max(m.modelResponse)-min(m.modelResponse));
    sd = std(m.baselineNoise(m.modelResponse>=lowbound & m.modelResponse<=highbound)); 
    avg = mean(m.baselineNoise(m.modelResponse>=lowbound & m.modelResponse<=highbound)); 
    med = median(m.baselineNoise(m.modelResponse>=lowbound & m.modelResponse<=highbound)); 

    sd = sd(~isnan(sd)); avg = avg(~isnan(avg)); med = med(~isnan(med));
    
    x = [x lowbound];
    y = [y sd];
    avgs = [avgs avg];
    meds = [meds med];
        
    yAll{length(x)} = m.baselineNoise(m.modelResponse>=lowbound & m.modelResponse<=highbound);
    avgsAll{length(x)} = m.baselineNoise(m.modelResponse>=lowbound & m.modelResponse<=highbound);
    
end

%% plot %%
subplot(2,2,2); hold on;
scatter(x*100,y*100,60,'black','filled'); 

xlim([min(x)*100,max(x)*100]);
title('Residual Std by Activity')
xlabel('BOLD signal (percent average)');
ylabel('Noise (% signal change)');

subplot(2,2,4); hold on;
scatter(m.modelResponse*100,m.baselineNoise*100,30,'red'); hold on;
scatter(x*100,avgs*100,60,'black','filled');
scatter(x*100,meds*100,60,'blue','filled');
plot([80,120],[0,0],'black--'); 

xlim([min(x)*100,max(x)*100]);
title('Residual by model-predicted activity')
xlabel('BOLD signal (percent average)');
ylabel('Residual (% signal change)');
