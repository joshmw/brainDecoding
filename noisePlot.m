function noisePlot(varargin)

% grab arguments by on if you called via command line or interrogator
if length(varargin)==7;
x = varargin{4}; y = varargin{5}; z = varargin{6};
elseif length(varargin) == 3;
length(varargin)==7;
x = varargin{1}; y = varargin{2}; z = varargin{3};
end

% grab other scan info
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

subplot(2,3,1)
scatter(1:length(m.tSeries),m.tSeries*100,'black'); hold on;
plot(1:length(m.modelResponse),m.modelResponse*100,'black')

title(sprintf('[%i %i %i] r^2=%0.2f eccentricity=%0.2f rfHalfWidth=%0.2f %s',x,y,z,thisR2,thisEccentricity,thisRfHalfWidth,a.params.pRFFit.rfType));
xlabel('Time');
ylabel('BOLD activity');


subplot(2,3,4)
plot(1:length(m.tSeries),m.baselineNoise*100,'black'); hold on;
title('Residual timecourse');
xlabel('Time');
ylabel('Residual (% signal change)');

%% graph slope %%
subplot(2,3,3)

slope = m.modelResponse;
    for elem = 2:239;
        slope(elem) = (m.modelResponse(elem+1) - m.modelResponse(elem-1))/2;  
    end

scatter(slope(2:239)*100,m.baselineNoise(2:239)*100);
P = polyfit(slope(2:239)*100,m.baselineNoise(2:239)*100,1); yfit = P(1)*slope(2:239)*100+P(2);
    hold on; plot(slope(2:239)*100,yfit,'b-.');
title('Residuals by slope of pRF model'); xlabel('Slope of pRF Model ([t+1] - [t-1])'); ylabel('Residual (% activity)');

subplot(2,3,6); scatter(abs(slope(2:239))*100,m.baselineNoise(2:239)*100)
title('Residuals by Absolute slope of pRF model'); xlabel('Absolute value slope of pRF Model ([t+1] - [t-1])'); ylabel('Residual (% activity)');    

subplot(2,3,4);plot((2:239),slope(2:239)*100, 'red');

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
subplot(2,3,5); hold on;
scatter(x*100,y*100,60,'black','filled'); 

xlim([min(x)*100,max(x)*100]);
title('Residual Std by Activity')
xlabel('BOLD signal (percent average)');
ylabel('Noise (% signal change)');

subplot(2,3,2); hold on;
scatter(m.modelResponse*100,m.baselineNoise*100,30,'red'); hold on;
scatter(x*100,avgs*100,60,'black','filled');
scatter(x*100,meds*100,60,'blue','filled');
plot([80,120],[0,0],'black--'); 

xlim([min(x)*100,max(x)*100]);
title('Residual by model-predicted activity')
xlabel('BOLD signal (percent average)');
ylabel('Residual (% signal change)');
