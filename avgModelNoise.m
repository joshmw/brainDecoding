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
 

function [rois cleanRois] = shufflecor(varargin)       

%% Load Data %%

getArgs(varargin);

load(data1); cleanRoisA = rois;
load(data2); cleanRoisB = rois;






% parse arguments 
startparams(1) = 1; startparams(2) = 0;
opts = optimset('display','off','maxIter',1000000,'MaxFunEvals',1000000);
minsearch = [0 -inf]; maxsearch = [inf inf];

%% iterate through fixed multiplicative noise scales, just varying additive noise, and find the additive noise magnitudes and likelihoods
% multiplicative scale to iterate through
multiplicativeScales = -.5:.5:.5;
keyboard
% go through voxels 
for voxel = 1:cleanRoisA(1).n
    iteration = 1;
    for multiplicativeScale = multiplicativeScales;
     [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
         lsqnonlin(@fitNoiseParameters,startparams(1),minsearch(1),maxsearch(1),opts,cleanRoisA,cleanRoisB,voxel,0,1,multiplicativeScale);

    additiveFit.parameters{voxel}{iteration} = params; additiveFit.resnorms{voxel}{iteration} = resnorm; additiveFit.residual{voxel}{iteration} = residual(1)-1000000; additiveFit.exitflag{voxel}{iteration} = exitflag; additiveFit.output{voxel}{iteration} = output; additiveFit.lambda{voxel}{iteration} = lambda; additiveFit.jacobian{voxel}{iteration} = jacobian;
    iteration = iteration+1;
    end
end
keyboard
figure, hold on, averageResiduals = 0;
for voxel = 1:cleanRoisA(1).n
    plot(multiplicativeScales,cell2mat(additiveFit.residual{voxel}{1}),'lineWidth',1,'color',[.8 .8 .8])
    averageResiduals = averageResiduals + cell2mat(additiveFit.residual{voxel});
end
plot(multiplicativeScales,averageResiduals/cleanRoisA(1).n,'lineWidth',3,'color',[0 0 0])
xlabel('Multiplicative Noise Scale'),ylabel('Log Likelihood'),title('Log Likelihoods: Fit additive noise with varying multiplicative scaling');




















%%%%%%%%%%%%%%%%%%%%%%%%
%% fitNoiseParameters %%
%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = fitNoiseParameters(params,cleanRoisA,cleanRoisB,voxel,fitAll,additiveFit,multiplicativeScale)

% Parsing parameters and get tseries length
if fitAll; additiveNoiseStd = params(1); multiplicativeNoiseScale = params(2); end
if additiveFit; additiveNoiseStd = params(1); multiplicativeNoiseScale = multiplicativeScale; end
timeSeriesLength = cleanRoisA(1).nFrames;

% get residual T series as mean model Tseries - measured single Tseries
residualTSeries = cleanRoisA(1).vox.tSeries(:,voxel) - cleanRoisB(1).vox.tSeries(:,voxel);

% calculate Std at every time point as a function of the model BOLD activity. Second term is a variance term, so take the sqrt (poisson mean-variance relationship)
noiseStdTimeSeries = additiveNoiseStd*ones(1,timeSeriesLength) + multiplicativeNoiseScale*cleanRoisA(1).vox.pRFtSeries(:,voxel);

% compute the likelihood of each residual time point %
logLikelihood = sum(-log(normpdf(residualTSeries,0,noiseStdTimeSeries)))+1000000;


