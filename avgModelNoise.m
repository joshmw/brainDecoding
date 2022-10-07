%% avgModelNoise.m %%
%
%       By: Josh Wilson
%       Created: June 2022
%
%   This script is supposed to show that there isn't stimulus-driven multiplicative noise in fmri - that the data
%   aren't noisier when there is a signal present or the voxel is reponding. You can get the log likelihood of the data by fitting some
%   baseline variance and a scaling parameter that lets the variance vary as a function of the signal-RF overlap. The model best fits the data
%   when that extra parameter is 0 - meaning that the voxel isn't more or less noisy when there is stimulus overlap.
%
%   The scaling parameter values you test are held in the "multiplicativeScales" variable. If you want a publication graph, set
%   that really high (500 values or so) otherwise just set like 5 or else it takes forever and the file is like 10 GB.
%
%   The data2 input is the mean model; data1 is the true time series; base is the original pRF you are using for cutoffs. This script finds which
%   voxels the original pRF puts in the original cleanRois for the base prf and uses those same voxels in the independent scans.
%
%   Usage: [scanA, scanB, additiveFit, fullFit] = avgModelNoise('data1=s0403mc345GaussianHdrNM.mat','data2=s0403mc12GaussianHdrNM.mat','base=s0403prf.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scanA, scanB, additiveFit, fullFit] = avgModelNoise(varargin);       



%%%%%%%%%%%%%%%%%%
%% Prepare Data %%
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%% Load data %%
%%%%%%%%%%%%%%%

getArgs(varargin);

load(data1);
roisA = rois;
load(data2);
roisB = rois;
load(base,'rois','cleanRois');


% get the other scan clean rois
for roi = 1:length(rois)
cleanRoisA(roi).vox.linearCoords = roisA(roi).vox.linearCoords(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.rfHalfWidth = roisA(roi).vox.rfHalfWidth(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.r2 = roisA(roi).vox.r2(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.polarAngle = roisA(roi).vox.polarAngle(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.eccentricity = roisA(roi).vox.eccentricity(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.tSeries = roisA(roi).vox.tSeries(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.measurementVar = roisA(roi).vox.measurementVar(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.baselineNoise = roisA(roi).vox.baselineNoise(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.x = roisA(roi).vox.params(1,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.y = roisA(roi).vox.params(2,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.rfstd = roisA(roi).vox.params(3,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.pRFtSeries = roisA(roi).vox.pRFtSeries(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).vox.scanCoords = roisA(roi).scanCoords(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisA(roi).nVoxels = length(cleanRoisA(roi).vox.linearCoords);
cleanRoisA(roi).nFrames = length(cleanRoisA(roi).vox.pRFtSeries(:,1))
end

for roi = 1:length(rois)
cleanRoisB(roi).vox.linearCoords = roisB(roi).vox.linearCoords(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.rfHalfWidth = roisB(roi).vox.rfHalfWidth(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.r2 = roisB(roi).vox.r2(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.polarAngle = roisB(roi).vox.polarAngle(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.eccentricity = roisB(roi).vox.eccentricity(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.tSeries = roisB(roi).vox.tSeries(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.measurementVar = roisB(roi).vox.measurementVar(ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.baselineNoise = roisB(roi).vox.baselineNoise(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.x = roisB(roi).vox.params(1,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.y = roisB(roi).vox.params(2,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.rfstd = roisB(roi).vox.params(3,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.pRFtSeries = roisB(roi).vox.pRFtSeries(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).vox.scanCoords = roisB(roi).scanCoords(:,ismember(rois(roi).vox.linearCoords,cleanRois(roi).vox.linearCoords));
cleanRoisB(roi).nVoxels = length(cleanRoisB(roi).vox.linearCoords);
cleanRoisB(roi).nFrames = length(cleanRoisB(roi).vox.pRFtSeries(:,1))
end


roiNames = {'V1';'V2';'V3'};

for roi = 1:length(cleanRoisA)
scanA(roi).r2 = cleanRoisA(roi).vox.r2; scanB(roi).r2 = cleanRoisB(roi).vox.r2;
scanA(roi).tSeries = cleanRoisA(roi).vox.tSeries'; scanB(roi).tSeries = cleanRoisB(roi).vox.tSeries';
scanA(roi).pRFtSeries = cleanRoisA(roi).vox.pRFtSeries'; scanB(roi).pRFtSeries = cleanRoisB(roi).vox.pRFtSeries';
scanA(roi).x = cleanRoisA(roi).vox.x; scanB(roi).x = cleanRoisB(roi).vox.x;
scanA(roi).y = cleanRoisA(roi).vox.y; scanB(roi).y = cleanRoisB(roi).vox.y;
scanA(roi).std = cleanRoisA(roi).vox.rfstd; scanB(roi).std = cleanRoisB(roi).vox.rfstd;
scanA(roi).nVoxels = cleanRoisA(roi).nVoxels;scanB(roi).nVoxels = cleanRoisB(roi).nVoxels;
scanA(roi).nFrames = cleanRoisA(roi).nFrames;scanB(roi).nFrames = cleanRoisB(roi).nFrames;
end

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiplicative noise fitting %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fit additive noise to fixed multiplicative scaling values and plot log likelihoods %%

% parse arguments 
startparams(1) = 1;
opts = optimset('display','off','maxIter',1000000,'MaxFunEvals',1000000,'DiffMinChange',0);
minsearch = [0]; maxsearch = [inf];
multiplicativeScales = -1:1:1;

% go through each roi and voxel individually
for roi = 1:length(scanA)
for voxel = 1:scanA(roi).nVoxels;
    iteration = 1;

    % go through the multiplicative scales and fit the additive noise component to find log likelihood
    for multiplicativeScale = multiplicativeScales;
     [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
         lsqnonlin(@fitNoiseParameters,startparams(1),minsearch(1),maxsearch(1),opts,scanA,scanB,voxel,roi,0,1,multiplicativeScale);

    % save the parameters, likelihoods, exit flag, jacobian
    additiveFit(roi).parameters{voxel}{iteration} = params; additiveFit(roi).resnorms{voxel}{iteration} = resnorm; additiveFit(roi).residual{voxel}{iteration} = residual(1)-10000; additiveFit(roi).exitflag{voxel}{iteration} = exitflag; additiveFit(roi).output{voxel}{iteration} = output; additiveFit(roi).lambda{voxel}{iteration} = lambda; additiveFit(roi).jacobian{voxel}{iteration} = jacobian;
    iteration = iteration+1;
    end %scales
end %voxels
additiveFit(roi).multiplicativeScales = multiplicativeScales;
end %rois


%% graph the log likelihoods for each voxel at each multi scale %%
gcf = imgcf;

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
        zerovalue = cell2mat(additiveFit(roi).residual{voxel});
        zerovalue = zerovalue(round(length(zerovalue)/2));
        plot(additiveFit(roi).multiplicativeScales,cell2mat(additiveFit(roi).residual{voxel})-zerovalue,'lineWidth',1,'color',[1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4 1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4  1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4])
        
        % add to average residual
        averageResiduals = averageResiduals + cell2mat(additiveFit(roi).residual{voxel})-zerovalue;
        totalGraphedVoxels = totalGraphedVoxels+1;

    end %cutoff check
end %voxel iteration

% plot the average 
plot(additiveFit(roi).multiplicativeScales,averageResiduals/totalGraphedVoxels,'lineWidth',4,'color',[0 0 0])
xlabel('Multiplicative Noise Scale'),ylabel('Log likelihood difference from 0 scale'),title(sprintf('%s Log likelihoods',roiNames{roi})),
plot([min(additiveFit(1).multiplicativeScales) max(additiveFit(1).multiplicativeScales)],[0 0],'r');
ylim([-10 60]);
drawPublishAxis('labelFontSize=14');legend('off');
end
set(gcf,'renderer','Painters')
set(gcf,'Position', [58.7375 5.22111111111111 30.8680555555556 20.6022222222222]);


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
    fullFit(roi).parameters{voxel} = params; fullFit(roi).resnorms{voxel} = resnorm; fullFit(roi).residual{voxel} = residual(1)-100000; fullFit(roi).exitflag{voxel} = exitflag; fullFit(roi).output{voxel} = output; fullFit(roi).lambda{voxel} = lambda; fullFit(roi).jacobian{voxel} = jacobian;
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




%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal dynamics %%
%%%%%%%%%%%%%%%%%%%%%%%

previousMeasurementWeights = -1:1:1;

% go through each roi and voxel individually
for roi = 1:length(scanA)
for voxel = 1:scanA(roi).nVoxels;
    iteration = 1;

    % go through the multiplicative scales and fit the additive noise component to find log likelihood
    for previousMeasurementWeight = previousMeasurementWeights;

        % get residual T series as mean model Tseries - measured single Tseries
        residualTSeries = scanA(roi).tSeries(voxel,:) - scanB(roi).tSeries(voxel,:);

        %cut off first value
        residualTSeriesX = residualTSeries(2:240);

        %Mu value is the same as X value but 1 earler
        residualTSeriesMu = residualTSeries(1:239);

        % compute the likelihood of each residual time point %
        logLikelihood = sum(-log(normpdf(residualTSeriesX,previousMeasurementWeight*residualTSeriesMu,std(residualTSeries))));


    % save the likelihood
    temporalShift(roi).likelihood{voxel}{iteration} = logLikelihood; 
    iteration = iteration+1;
    end %scales

end %voxels
temporalShift(roi).previousMeasurementWeights = previousMeasurementWeights;
end %rois


%% graph the log likelihoods for each voxel at each weight%%
gcf = imgcf;

% go through every roi
for roi = 1:length(scanA)
subplot(1,length(scanA),roi); averageResiduals = 0; hold on;

% start count for total graphed voxels and set the r2 cutoff you want to graph
graphMinCutoff = .2; graphMaxCutoff=1; totalGraphedVoxels = 0;

% plot voxels with r2 greater than graphCutoff and tally number of plotted voxels
for voxel = 1:scanA(roi).nVoxels

    %check if over min and under max r2 cutoffs
    if ((scanA(roi).r2(voxel) > graphMinCutoff & scanB(roi).r2(voxel) > graphMinCutoff) & ...
       (scanA(roi).r2(voxel) < graphMaxCutoff & scanB(roi).r2(voxel) < graphMaxCutoff))

        % if so, plot 
        zerovalue = cell2mat(temporalShift(roi).likelihood{voxel});
        zerovalue = zerovalue(round(length(zerovalue)/2));
        plot(temporalShift(roi).previousMeasurementWeights,cell2mat(temporalShift(roi).likelihood{voxel})-zerovalue,'lineWidth',1,'color',[1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4 1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4  1-(scanA(roi).r2(voxel)+scanB(roi).r2(voxel))/4])
        
        % add to average residual
        averageResiduals = averageResiduals + cell2mat(temporalShift(roi).likelihood{voxel})-zerovalue;
        totalGraphedVoxels = totalGraphedVoxels+1;

    end %cutoff check
end %voxel iteration

% plot the average 
plot(temporalShift(roi).previousMeasurementWeights,averageResiduals/totalGraphedVoxels,'lineWidth',4,'color',[0 0 0])
xlabel('Weight of t-1 noise value'),ylabel('Log likelihood difference from 0 weight'),title(sprintf('%s Log likelihoods',roiNames{roi})),
plot([min(temporalShift(1).previousMeasurementWeights) max(temporalShift(1).previousMeasurementWeights)],[0 0],'r');
drawPublishAxis('labelFontSize=14');legend('off');
end
set(gcf,'renderer','Painters')
set(gcf,'Position', [58.7375 5.22111111111111 30.8680555555556 20.6022222222222]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean model correlations in different areas %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% rf overlap %%%
roi = 1; for row = 1:length(scanA(roi).r2)
    for column = 1:length(scanA(roi).r2)   
        mu1 = 0; s1 = scanA(roi).std(column); s2 = scanA(roi).std(row);
        mu2 = sqrt((scanA(roi).x(column)-scanA(roi).x(row))^2 + (scanA(roi).y(column)-scanA(roi).y(row))^2);
        if mu1 == mu2 & s1 == s2; v1rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v1rfOverlap = v1rfOverlap.*100;

% v2
roi = 2; for row = 1:length(scanA(roi).r2)
    for column = 1:length(scanA(roi).r2)   
        mu1 = 0; s1 = scanA(roi).std(column); s2 = scanA(roi).std(row);
        mu2 = sqrt((scanA(roi).x(column)-scanA(roi).x(row))^2 + (scanA(roi).y(column)-scanA(roi).y(row))^2);
        if mu1 == mu2 & s1 == s2; v2rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v2rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v2rfOverlap = v2rfOverlap.*100;

% v3
roi = 3; for row = 1:length(scanA(roi).r2)
    for column = 1:length(scanA(roi).r2)   
        mu1 = 0; s1 = scanA(roi).std(column); s2 = scanA(roi).std(row);
        mu2 = sqrt((scanA(roi).x(column)-scanA(roi).x(row))^2 + (scanA(roi).y(column)-scanA(roi).y(row))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v3rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v3rfOverlap = v3rfOverlap.*100;

%v1 v3
 for row = 1:length(scanA(1).r2)
    for column = 1:length(scanA(3).r2)   
        mu1 = 0; s1 = scanA(3).std(column); s2 = scanA(1).std(row);
        mu2 = sqrt((scanA(3).x(column)-scanA(1).x(row))^2 + (scanA(3).y(column)-scanA(1).y(row))^2);
        if mu1 == mu2 & s1 == s2; v1v3rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1v3rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
 end
 v1v3rfOverlap = v1v3rfOverlap.*100;

%% noise correlations %%
% calculate noise correlations in scan A by subtracting the mean model of scan B
v1NoiseCor = corr(scanA(1).tSeries'-scanB(1).tSeries');
v2NoiseCor = corr(scanA(2).tSeries'-scanB(2).tSeries');
v3NoiseCor = corr(scanA(3).tSeries'-scanB(3).tSeries');
v1v3NoiseCor = corr(scanA(1).tSeries'-scanB(1).tSeries',scanA(3).tSeries'-scanB(3).tSeries');



%% graph %%
% first, do v1 %
v1rfOverlapArr = reshape(v1rfOverlap,[1 length(v1rfOverlap)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1rfOverlapArr,sortOrder] = sort(v1rfOverlapArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

v1NoiseCorArr(v1rfOverlapArr==100) = []; v1rfOverlapArr(v1rfOverlapArr==100) = [];

figure(11); hold on; scatter(v1rfOverlapArr,v1NoiseCorArr,1,'filled','k');

%bootstrap
getFitConfidenceInterval(v1NoiseCorArr',v1rfOverlapArr');

expFit = fit(v1rfOverlapArr',v1NoiseCorArr','exp1');
legendColor = plot(expFit); legendColor.Color = [.3, .5, .3]; legendColor.LineWidth = 2;
v1expFit = plot(expFit); v1expFit.Color = [0, 0.4470, 0.7410]; v1expFit.LineWidth = 2;
%confidence internal
confInt = confint(expFit);
lowerFitBound = @(x) confInt(1,1)*exp(confInt(1,2)*x); upperFitBound = @(x) confInt(2,1)*exp(confInt(2,2)*x);
fplot(lowerFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]); fplot(upperFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]);

legend('Exponential Fit');
title('V1 residual correlations'); xlabel('Receptive field overlap between voxels (percent)'); ylabel('Residual correlation between voxels (Pearson r)'); ylim([-.4 1]);

drawPublishAxis('labelFontSize=14');set(gcf,'renderer','Painters')

leg = legend([v1expFit legendColor],'Exponential Fit', 'Bootstrapped Fits'); leg.Position = [0.17 .77 0.2685 0.1003];

% second, do v2 %
v2rfOverlapArr = reshape(v2rfOverlap,[1 length(v2rfOverlap)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2rfOverlapArr,sortOrder] = sort(v2rfOverlapArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

v2NoiseCorArr(v2rfOverlapArr==100) = []; v2rfOverlapArr(v2rfOverlapArr==100) = [];

figure(12); hold on; scatter(v2rfOverlapArr,v2NoiseCorArr,1,'filled','k');

%bootstrap
getFitConfidenceInterval(v2NoiseCorArr',v2rfOverlapArr');

expFit = fit(v2rfOverlapArr',v2NoiseCorArr','exp1');
legendColor = plot(expFit); legendColor.Color = [.3, .5, .3]; legendColor.LineWidth = 2;
v2expFit = plot(expFit); v2expFit.Color = [0, 0.4470, 0.7410]; v2expFit.LineWidth = 2;
%confidence internal
confInt = confint(expFit);
lowerFitBound = @(x) confInt(1,1)*exp(confInt(1,2)*x); upperFitBound = @(x) confInt(2,1)*exp(confInt(2,2)*x);
fplot(lowerFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]); fplot(upperFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]);

legend('Exponential Fit');
title('V2 residual correlations'); xlabel('Receptive field overlap between voxels (percent)'); ylabel('Residual correlation between voxels (Pearson r)');ylim([-.4 1]);

drawPublishAxis('labelFontSize=14');set(gcf,'renderer','Painters')

leg = legend([v2expFit legendColor],'Exponential Fit', 'Bootstrapped Fits'); leg.Position = [0.17 .77 0.2685 0.1003];

% third, do v3 %
v3rfOverlapArr = reshape(v3rfOverlap,[1 length(v3rfOverlap)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3rfOverlapArr,sortOrder] = sort(v3rfOverlapArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

v3NoiseCorArr(v3rfOverlapArr==100) = []; v3rfOverlapArr(v3rfOverlapArr==100) = [];

figure(13); hold on; scatter(v3rfOverlapArr,v3NoiseCorArr,1,'filled','k');

%bootstrap
getFitConfidenceInterval(v3NoiseCorArr',v3rfOverlapArr');

expFit = fit(v3rfOverlapArr',v3NoiseCorArr','exp1');
legendColor = plot(expFit); legendColor.Color = [.3, .5, .3]; legendColor.LineWidth = 2;
v3expFit = plot(expFit); v3expFit.Color = [0, 0.4470, 0.7410]; v3expFit.LineWidth = 2;
%confidence internal
confInt = confint(expFit);
lowerFitBound = @(x) confInt(1,1)*exp(confInt(1,2)*x); upperFitBound = @(x) confInt(2,1)*exp(confInt(2,2)*x);
fplot(lowerFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]); fplot(upperFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]);

legend('Exponential Fit');
title('V3 residual correlations'); xlabel('Receptive field overlap between voxels (percent)'); ylabel('Residual correlation between voxels (Pearson r)');ylim([-.4 1]);

drawPublishAxis('labelFontSize=14');set(gcf,'renderer','Painters')

leg = legend([v3expFit legendColor],'Exponential Fit', 'Bootstrapped Fits'); leg.Position = [0.17 .77 0.2685 0.1003];

%between v1 and v3
v1v3rfOverlapArr = reshape(v1v3rfOverlap,[1 min(size(v1v3rfOverlap))*max(size(v1v3rfOverlap))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3rfOverlapArr,sortOrder] = sort(v1v3rfOverlapArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

v1v3NoiseCorArr(v1v3rfOverlapArr==100) = []; v1v3rfOverlapArr(v1v3rfOverlapArr==100) = [];

v1v3NoiseCorArr(isnan(v1v3rfOverlapArr))=[];
v1v3rfOverlapArr(isnan(v1v3rfOverlapArr))=[];

figure(15); hold on; scatter(v1v3rfOverlapArr,v1v3NoiseCorArr,1,'filled','k');

getFitConfidenceInterval(v1v3NoiseCorArr',v1v3rfOverlapArr');

expFit = fit(v1v3rfOverlapArr',v1v3NoiseCorArr','exp1');
legendColor = plot(expFit); legendColor.Color = [.3, .5, .3]; legendColor.LineWidth = 2;
v1v3expFit = plot(expFit); v1v3expFit.Color = [0, 0.4470, 0.7410]; v1v3expFit.LineWidth = 2;
%confidence internal
confInt = confint(expFit);
lowerFitBound = @(x) confInt(1,1)*exp(confInt(1,2)*x); upperFitBound = @(x) confInt(2,1)*exp(confInt(2,2)*x);
fplot(lowerFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]); fplot(upperFitBound,[0,100],'--','color',[0, 0.4470, 0.7410]);

legend('Exponential Fit');
title('V1/V3 residual correlations'); xlabel('Receptive field overlap between voxels (percent)'); ylabel('Residual correlation between voxels (Pearson r)');ylim([-.4 1]);

drawPublishAxis('labelFontSize=14');set(gcf,'renderer','Painters')

leg = legend([v1v3expFit legendColor],'Exponential Fit', 'Bootstrapped Fits'); leg.Position = [0.17 .77 0.2685 0.1003];





keyboard



%% helper functions %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getFitConfidenceInterval %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getFitConfidenceInterval(correlationArray,overlapArray);

numBootstraps = 50;

expFitA = zeros(1,numBootstraps); expFitB = zeros(1,numBootstraps); expFitMax = zeros(1,numBootstraps);

for bootstrap = 1:numBootstraps;
    shuffledOverlap = overlapArray(randperm(length(overlapArray)));
    expFit = fit(shuffledOverlap,correlationArray,'exp1');
    expFitA(bootstrap) = expFit.a;
    expFitB(bootstrap) = expFit.b;
    expFitMax(bootstrap) = expFit.a*exp(expFit.b);
    
    legend('off');
    shuffleFit = plot(expFit);
    shuffleFit.Color = [.3, .5, .3 .05]; shuffleFit.LineWidth = 3;

end

[expFitA, sortOrder] = sort(expFitA);
expFitB = expFitB(sortOrder);
expFitMax = expFitMax(sortOrder);



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
shiftedpRFtSeries = (scanA(roi).pRFtSeries(voxel,:)-min(scanA(roi).pRFtSeries(voxel,:)));
multiplicativeNoiseSeries = multiplicativeNoiseScale*shiftedpRFtSeries*std(residualTSeries)/(max(shiftedpRFtSeries));

% calculate Std at every time point as a function of the model BOLD activity. Second term is a variance term, so take the sqrt (poisson mean-variance relationship) (took out for now)
noiseStdTimeSeries = additiveNoiseStd*ones(1,timeSeriesLength) + multiplicativeNoiseSeries;
% compute the likelihood of each residual time point %
logLikelihood = sum(-log(normpdf(residualTSeries,0,noiseStdTimeSeries)))+10000;


