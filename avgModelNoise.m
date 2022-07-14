%% avgModelNoise.m %%
%
%       By: Josh Wilson
%       Created: June 2022
%
%   Subtracts the mean model (1 set of avg scans) from another scan and plots the residuals.
%   At some point, I will do a model comparison of the types of noise.
%
%   The data2 input is the mean model; data1 is the true time series. You can set r2 cutoffs (and other stuff) in the script.
%
%
%   Usage: avgModelNoise('data1=s0401mc12hrfFit.mat','data2=s0401mc345hrfFit.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scanA, scanB, additiveFit, fullFit] = avgModelNoise(varargin);       

%% Load Data %%
getArgs(varargin);
load(data1); cleanRoisA = rois;
load(data2); cleanRoisB = rois;



%% Clean Data %%
for roi = 1:length(cleanRoisA)

% set cutoffs
r2min = .6;
analysisCutoff = (cleanRoisA(roi).vox.r2 > r2min) & (cleanRoisB(roi).vox.r2 > r2min);

% create new data structures with cutoff stats
scanA(roi).r2 = cleanRoisA(roi).vox.r2(analysisCutoff); scanB(roi).r2 = cleanRoisB(roi).vox.r2(analysisCutoff);
scanA(roi).tSeries = cleanRoisA(roi).vox.tSeries(:,analysisCutoff)'; scanB(roi).tSeries = cleanRoisB(roi).vox.tSeries(:,analysisCutoff)';
scanA(roi).pRFtSeries = cleanRoisA(roi).vox.pRFtSeries(:,analysisCutoff)'; scanB(roi).pRFtSeries = cleanRoisB(roi).vox.pRFtSeries(:,analysisCutoff)';

% get number of voxels and number of frames
scanA(roi).nFrames = cleanRoisA(roi).nFrames; scanB(roi).nFrames = cleanRoisB(roi).nFrames;
scanA(roi).nVoxels = sum(analysisCutoff); scanB(roi).nVoxels = sum(analysisCutoff);
end



%% Fit additive noise to fixed multiplicative scaling values and plot log likelihoods %%

% parse arguments 
startparams(1) = 1;
opts = optimset('display','off','maxIter',1000000,'MaxFunEvals',1000000,'DiffMinChange',0);
minsearch = [0]; maxsearch = [inf];
multiplicativeScales = -2:1:2;

% go through each roi and voxel individually
for roi = 1:length(scanA)
for voxel = 1:scanA(roi).nVoxels;
    iteration = 1;

    % go through the multiplicative scales and fit the additive noise component to find log likelihood
    for multiplicativeScale = multiplicativeScales;
     [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
         lsqnonlin(@fitNoiseParameters,startparams(1),minsearch(1),maxsearch(1),opts,scanA,scanB,voxel,roi,0,1,multiplicativeScale);

    % save the parameters, likelihoods, exit flag, jacobian
    additiveFit(roi).parameters{voxel}{iteration} = params; additiveFit(roi).resnorms{voxel}{iteration} = resnorm; additiveFit(roi).residual{voxel}{iteration} = residual(1)-1000; additiveFit(roi).exitflag{voxel}{iteration} = exitflag; additiveFit(roi).output{voxel}{iteration} = output; additiveFit(roi).lambda{voxel}{iteration} = lambda; additiveFit(roi).jacobian{voxel}{iteration} = jacobian;
    iteration = iteration+1;
    end %scales
end %voxels
end %rois



%% graph the loglikelihoods for each voxel at each multi scale %%
figure

% go through every roi
for roi = 1:length(scanA)
subplot(1,length(scanA),roi); averageResiduals = 0; hold on;

% start count for total graphed voxels and set the r2 cutoff you want to graph
graphCutoff = 0; totalGraphedVoxels = 0

% plot voxels with r2 greater than graphCutoff and tally number of plotted voxels
for voxel = 1:scanA(roi).nVoxels
    if (scanA(roi).r2(voxel) > graphCutoff & scanB(roi).r2(voxel > graphCutoff))
        plot(multiplicativeScales,cell2mat(additiveFit(roi).residual{voxel}),'lineWidth',1,'color',[.8 .8 .8])
        averageResiduals = averageResiduals + cell2mat(additiveFit(roi).residual{voxel});
        totalGraphedVoxels = totalGraphedVoxels+1;
    end %cutoff check
end %voxel iteration

% plot the average 
plot(multiplicativeScales,averageResiduals/totalGraphedVoxels,'lineWidth',3,'color',[0 0 0])
xlabel('Multiplicative Noise Scale'),ylabel('Log Likelihood'),title('Log likelihoods with fit additive noise');
end



%% find best parameters for every voxel, fitting additive noise and multiplicative scaling %%

% parse arguments 
startparams(1) = .1; startparams(2) = 0;
opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000);
minsearch = [0 -inf]; maxsearch = [inf inf];

% go through each voxel individually
for roi = 1:length(scanA)
for voxel = 1:scanA(roi).nVoxels;

    % fit the additive noise component and multiplicative scaling to find log likelihood
     [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
         lsqnonlin(@fitNoiseParameters,startparams,minsearch,maxsearch,opts,scanA,scanB,voxel,roi,1,0,0);

    % save the parameters, likelihoods, exit flag, jacobian
    fullFit(roi).parameters{voxel} = params; fullFit(roi).resnorms{voxel} = resnorm; fullFit(roi).residual{voxel} = residual(1)-1000; fullFit(roi).exitflag{voxel} = exitflag; fullFit(roi).output{voxel} = output; fullFit(roi).lambda{voxel} = lambda; fullFit(roi).jacobian{voxel} = jacobian;
end %voxel iteration
end %roi iteration  




%% END OF SCRIPT %%











%% helper functions %%


%%%%%%%%%%%%%%%%%%%%%%%%
%% fitNoiseParameters %%
%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = fitNoiseParameters(params,scanA,scanB,voxel,roi,fitAll,additiveFit,multiplicativeScale)

% Parsing parameters and get tseries length
if fitAll;additiveNoiseStd = params(1); multiplicativeNoiseScale = params(2); end
if additiveFit; additiveNoiseStd = params(1); multiplicativeNoiseScale = multiplicativeScale; end
timeSeriesLength = scanA(roi).nFrames;

% get residual T series as mean model Tseries - measured single Tseries
residualTSeries = scanA(roi).tSeries(voxel,:) - scanB(roi).tSeries(voxel,:);

% calculated the multiplicative noise series scaled by the pRF-predicted response (zerod by min of prf prediction)
multiplicativeNoiseSeries = multiplicativeNoiseScale*(scanA(roi).pRFtSeries(voxel,:)-min(scanA(roi).pRFtSeries(voxel,:)));

% calculate Std at every time point as a function of the model BOLD activity. Second term is a variance term, so take the sqrt (poisson mean-variance relationship) (took out for now)
noiseStdTimeSeries = additiveNoiseStd*ones(1,timeSeriesLength) + multiplicativeNoiseSeries;

% compute the likelihood of each residual time point %
logLikelihood = sum(-log(normpdf(residualTSeries,0,noiseStdTimeSeries)))+1000;


