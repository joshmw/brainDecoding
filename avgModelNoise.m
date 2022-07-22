%% avgModelNoise.m %%
%
%       By: Josh Wilson
%       Created: June 2022
%
%   Look at noise characteristics by using a mean model from another scan. Requires 2 scans (averages of scans) as input.
%
%   The data2 input is the mean model; data1 is the true time series. You can set r2 cutoffs (and other stuff) in the script.
%
%
%
%
%
%
%   Usage: [scanA, scanB, additiveFit, fullFit] = avgModelNoise('data1=s0401mc12hrfFit.mat','data2=s0401mc345hrfFit.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scanA, scanB, additiveFit, fullFit] = avgModelNoise(varargin);       



%%%%%%%%%%%%%%%%%%
%% Prepare Data %%
%%%%%%%%%%%%%%%%%%

%% Load Data %%
getArgs(varargin);
load(data1); cleanRoisA = rois;
load(data2); cleanRoisB = rois;
roiNames = {'V1';'V2';'V3'};
correctSignalChangeScale = 1;


%% Clean Data %%
for roi = 1:length(cleanRoisA)

% set cutoffs
r2min = .2;
analysisCutoff = (cleanRoisA(roi).vox.r2 > r2min) & (cleanRoisB(roi).vox.r2 > r2min);

% create new data structures with cutoff stats
scanA(roi).r2 = cleanRoisA(roi).vox.r2(analysisCutoff); scanB(roi).r2 = cleanRoisB(roi).vox.r2(analysisCutoff);
scanA(roi).tSeries = cleanRoisA(roi).vox.tSeries(:,analysisCutoff)'; scanB(roi).tSeries = cleanRoisB(roi).vox.tSeries(:,analysisCutoff)';
scanA(roi).pRFtSeries = cleanRoisA(roi).vox.pRFtSeries(:,analysisCutoff)'; scanB(roi).pRFtSeries = cleanRoisB(roi).vox.pRFtSeries(:,analysisCutoff)';

% correct the scale if needed
if correctSignalChangeScale
    scanA(roi).tSeries = scanA(roi).tSeries/100+1; scanB(roi).tSeries = scanB(roi).tSeries/100+1;
    scanA(roi).pRFtSeries = scanA(roi).pRFtSeries/100+1; scanB(roi).pRFtSeries = scanB(roi).pRFtSeries/100+1;
end

% get number of voxels and number of frames
scanA(roi).nFrames = cleanRoisA(roi).nFrames; scanB(roi).nFrames = cleanRoisB(roi).nFrames;
scanA(roi).nVoxels = sum(analysisCutoff); scanB(roi).nVoxels = sum(analysisCutoff);
scanA(roi).r2min = r2min; scanB(roi).r2min = r2min;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiplicative noise fitting %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fit additive noise to fixed multiplicative scaling values and plot log likelihoods %%

% parse arguments 
startparams(1) = 1;
opts = optimset('display','off','maxIter',1000000,'MaxFunEvals',1000000,'DiffMinChange',0);
minsearch = [0]; maxsearch = [inf];
multiplicativeScales = -1:.25:1;

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
additiveFit(roi).multiplicativeScales = multiplicativeScales;
end %rois


%% graph the log likelihoods for each voxel at each multi scale %%
figure

% go through every roi
for roi = 1:length(scanA)
subplot(1,length(scanA),roi); averageResiduals = 0; hold on;

% start count for total graphed voxels and set the r2 cutoff you want to graph
graphMinCutoff = .6; graphMaxCutoff=1; totalGraphedVoxels = 0;

% plot voxels with r2 greater than graphCutoff and tally number of plotted voxels
for voxel = 1:scanA(roi).nVoxels

    %check if over min and under max r2 cutoffs
    if ((scanA(roi).r2(voxel) > graphMinCutoff & scanB(roi).r2(voxel) > graphMinCutoff) & ...
       (scanA(roi).r2(voxel) < graphMaxCutoff & scanB(roi).r2(voxel) < graphMaxCutoff))

        % if so, plot 
        plot(additiveFit(roi).multiplicativeScales,cell2mat(additiveFit(roi).residual{voxel}),'lineWidth',1,'color',[1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4 1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4  1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4])
        
        % add to average residual
        averageResiduals = averageResiduals + cell2mat(additiveFit(roi).residual{voxel});
        totalGraphedVoxels = totalGraphedVoxels+1;

    end %cutoff check
end %voxel iteration

% plot the average 
plot(additiveFit(roi).multiplicativeScales,averageResiduals/totalGraphedVoxels,'lineWidth',4,'color',[0 0 0])
xlabel('Multiplicative Noise Scale'),ylabel('Log Likelihood'),title(sprintf('%s Log likelihoods with fit additive noise',roiNames{roi})),
end


%% graph the std of the noise at binned activity levels %%
bins = 10; figure;

for roi = 1:length(scanA)

    % start count for total graphed voxels and set the r2 cutoff you want to graph
    graphMinCutoff = .6; graphMaxCutoff=1; totalGraphedVoxels = 0; avgStds = zeros(1,bins); binPlot = []; binStdPlot = [];
    subplot(1,3,roi); hold on;

    for voxel = 1:length(scanA(roi).r2)


        % check for r2 cutoffs
        if ((scanA(roi).r2(voxel) > graphMinCutoff & scanB(roi).r2(voxel) > graphMinCutoff) & ...
        (scanA(roi).r2(voxel) < graphMaxCutoff & scanB(roi).r2(voxel) < graphMaxCutoff))

            % count this voxel
            totalGraphedVoxels = totalGraphedVoxels +1;
            
            % iterate through bins
            for bin = 0:1/bins:1-1/bins
                    
                % find the percentage (by bin) cutoffs in the pRFtSeries
                lowPercentile = prctile(scanA(roi).pRFtSeries(voxel,:),bin*100);
                highPercentile = prctile(scanA(roi).pRFtSeries(voxel,:),(bin+(1/bins))*100);
                residualTSeries = scanA(roi).tSeries(voxel,:) - scanB(roi).tSeries(voxel,:);
                
                % calculate std of residuals in the bin boundaries
                binStd = std(residualTSeries((scanA(roi).pRFtSeries(voxel,:) >= lowPercentile) & (scanA(roi).pRFtSeries(voxel,:) <= highPercentile)));
                avgStds(round(bin*bins+1)) = avgStds(round(bin*bins+1)) + binStd;

                % plot
                binPlot = [binPlot bin]; 
                binStdPlot = [binStdPlot binStd];
                end % r2 check

        end % bin iteration

    end % voxel iteration
    
    % plot average
    scatter(binPlot,binStdPlot,1,'filled','k')
    scatter(0:1/bins:1-1/bins, avgStds/totalGraphedVoxels,75,'filled','k');

    % label stuff
    xlabel(['pRF activity bin']),ylabel('Std of residuals'); title(sprintf('Std of residuals by activity level: %S',roiNames{roi}));

end % roi iteration



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
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    % save the parameters, likelihoods, exit flag, jacobian
    fullFit(roi).parameters{voxel} = params; fullFit(roi).resnorms{voxel} = resnorm; fullFit(roi).residual{voxel} = residual(1)-1000; fullFit(roi).exitflag{voxel} = exitflag; fullFit(roi).output{voxel} = output; fullFit(roi).lambda{voxel} = lambda; fullFit(roi).jacobian{voxel} = jacobian;
end %voxel iteration
end %roi iteration  



%%%%%%%%%%%%%%%%%%%%%%
%% Autocorrelations %%
%%%%%%%%%%%%%%%%%%%%%%
zScoreTseries = 0; graphMinCutoff = .6; graphMaxCutoff=1;
figure

% iterate through rois
for roi = 1:length(scanA)

    % reset counts for each roi and change subplots
    subplot(2,length(scanA),roi); hold on, allCorrelations = 0; totalGraphedVoxels = 0;
    
    % go through every voxel
    for voxel = 1:scanA(roi).nVoxels
    
        % check it it meets r2 cutoff
        if ((scanA(roi).r2(voxel) > graphMinCutoff & scanB(roi).r2(voxel) > graphMinCutoff) & ...
           (scanA(roi).r2(voxel) < graphMaxCutoff & scanB(roi).r2(voxel) < graphMaxCutoff))
        
            % if so, initiate params
            lags = []; correlations = [];
        
            % calculate residual time series
            residualTSeries = scanA(roi).tSeries(voxel,:) - scanB(roi).tSeries(voxel,:);
            if zScoreTseries; residualTSeries = zscore(scanA(roi).tSeries(voxel,:)) - zscore(scanB(roi).tSeries(voxel,:)); end
        
            % iterate through lags and calculate correlation
            for lag = 0:scanA(roi).nFrames-1
                shiftedResidualTSeries = [residualTSeries(lag+1:end) residualTSeries(1:lag)];
                lags = [lags lag]; correlations = [correlations corr2(residualTSeries,shiftedResidualTSeries)];
            end
        
            % plot it
            plot(lags,correlations,'LineWidth',1,'color',[1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4 1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4  1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4])
        
            % add to group data
            allCorrelations = allCorrelations + correlations;
            totalGraphedVoxels = totalGraphedVoxels + 1;
            
        end %r2 cutoff check
    
    end %voxel iteration
    
    % plot group data
    plot(lags,allCorrelations/totalGraphedVoxels,'LineWidth',4,'color','k')
    title(sprintf('s01 Autocorrelations: %s',roiNames{roi}));xlabel('Lag (TRs)');ylabel('Correlation');
    
    % fourier transform
    subplot(2,length(scanA),roi+length(scanA));
    fftTSeries = fft(allCorrelations);
    % set mean to zero
    fftTSeries(1) = 0;
    % plot it 
    plot(1:(length(fftTSeries)/2)-1,abs(fftTSeries(2:length(fftTSeries)/2)),'k.-');
    xlabel('FFT components'); ylabel('FFT of fMRI Signal'); title('Fourier Components')
    axis tight; zoom on

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


