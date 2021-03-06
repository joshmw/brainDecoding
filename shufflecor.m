%% getpRFTSeries.m %%
%
%       By: Josh Wilson
%       Created: March 2022
%
% Takes a scan and analysis in mrTools, grabs info (time series, parameters, etc) from voxels, and graphs some stuff. 
% cleanRois is a data frame that indexes into different rois. You define what rois you want and the order in the script. Defaults to v1,v2,v3.
%
% Part 1: Saves voxel data in structures 'rois' (all  voxels) and 'cleanRois' (filtered version of 'rois' that meets 
% certain cutoffs you can define like r2, RF width, etc). Need to do unless you have structures saved already.
% 
% Part 2: Takes the cleanROIs data and graphs things like RF overlap, distance, noise correlations between voxels. If you 
% don't want to do this (ex. you just want to save the data), you can set "graphStuff = 0". But might as well graph if extracting for first time.
%
%   Arguments you probably want to pass in:
%
%   *loadData*: pass in 1 if you to load a saved structure, else 0 to extract new from mrTools analysis.
%   
%       If 1, also pass in:
%       *data*: name of file with rois + cleanRois (ex. 's0401pRF.mat'). Need to have rois + cleanRois saved in that file.
%
%       If 0, also pass in:
%       *scanNum*: Scan number. Defaults to the concat group. I always concat everything (even single scans) to filter them.
%       *analysis*: Name of the anaylsis you want (ex: 'pRFDoG', 'pRF').
%       You should also save the 'rois' and 'cleanRois' structures in a .mat file for easier loading/analysis later.
%
%
%       Example usage:
%           Pre-extracted data:              getpRFTSeries('loadData=1','data=s0401pRF.mat')
%           Unextracted (mrTools open):      [rois cleanRois] = getpRFTSeries('scanNum=2','analysis=pRFDoG'
%
%       
 

function [rois cleanRois] = shufflecor(varargin)       

%%%%%%%%%%%%%%%%%%%
%% Get Arguments %%
%%%%%%%%%%%%%%%%%%%


graphStuff = 1;

%%%%%%%%%%%%%%%
%% Load data %%
%%%%%%%%%%%%%%%

getArgs(varargin);

load(data1);
cleanRoisA = cleanRois;
load(data2);
cleanRoisB = cleanRois;


% Start of part 2 (calculating correlations, distances, then graphing).
%%%%%%%%%%%%%%%%
%% RF Overlap %%
%%%%%%%%%%%%%%%%

% receptive field overlaps %
sprintf('Correlation RFs and calculating distances (takes about a minute)...')

%v1
for roi1vox = 1:length(cleanRoisA(1).vox.linearCoords);
    for roi2vox = 1:length(cleanRoisB(1).vox.linearCoords);       
        mu1 = 0; s1 = cleanRoisA(1).vox.rfstd(roi1vox); s2 = cleanRoisB(1).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRoisA(1).vox.x(roi1vox)-cleanRoisB(1).vox.x(roi2vox))^2 + (cleanRoisA(1).vox.y(roi1vox)-cleanRoisB(1).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v1rfOverlap = transpose(v1Overlap);

%v2 
for roi1vox = 1:length(cleanRoisA(2).vox.linearCoords);
    for roi2vox = 1:length(cleanRoisB(2).vox.linearCoords);       
        mu1 = 0; s1 = cleanRoisA(2).vox.rfstd(roi1vox); s2 = cleanRoisB(2).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRoisA(2).vox.x(roi1vox)-cleanRoisB(2).vox.x(roi2vox))^2 + (cleanRoisA(2).vox.y(roi1vox)-cleanRoisB(2).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v2Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v2rfOverlap = transpose(v2Overlap);

%v3
for roi1vox = 1:length(cleanRoisA(3).vox.linearCoords);
    for roi2vox = 1:length(cleanRoisB(3).vox.linearCoords);       
        mu1 = 0; s1 = cleanRoisA(3).vox.rfstd(roi1vox); s2 = cleanRoisB(3).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRoisA(3).vox.x(roi1vox)-cleanRoisB(3).vox.x(roi2vox))^2 + (cleanRoisA(3).vox.y(roi1vox)-cleanRoisB(3).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v3Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v3rfOverlap = transpose(v3Overlap);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get noise correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1NoiseCor = transpose(corr(cleanRoisA(1).vox.baselineNoise,cleanRoisB(1).vox.baselineNoise));
v2NoiseCor = transpose(corr(cleanRoisA(2).vox.baselineNoise,cleanRoisB(2).vox.baselineNoise)); 
v3NoiseCor = transpose(corr(cleanRoisA(3).vox.baselineNoise,cleanRoisB(3).vox.baselineNoise)); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %

v1rfOverlapArr = reshape(v1rfOverlap,[1 min(size(v1rfOverlap))*max(size(v1rfOverlap))]);
v1NoiseCorArr = reshape(v1NoiseCor,[1 min(size(v1NoiseCor))*max(size(v1NoiseCor))]);
[v1rfOverlapArr,sortOrder] = sort(v1rfOverlapArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

v1NoiseCorArr(v1rfOverlapArr==1) = []; v1rfOverlapArr(v1rfOverlapArr==1) = [];

figure(11); hold on; scatter(v1rfOverlapArr,v1NoiseCorArr,1,'filled','k');

expFit = fit(v1rfOverlapArr',v1NoiseCorArr','exp1');
v1expFit = plot(expFit,'predobs'); for i = 1:3, v1expFit(i).Color = [0, 0.4470, 0.7410]; v1expFit(i).LineWidth = 2; end
for i = 2:3, v1expFit(i).LineStyle = '--'; v1expFit(i).LineWidth = .75; end

legend('','Linear Fit', 'Exponential Fit');
title('V1 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j');

drawPublishAxis('labelFontSize=14');
leg = legend('','Exponential Fit', '95% Prediction bounds'); leg.Position = [0.6 0.2 0.2685 0.1003];

% second, do v2 %
v2rfOverlapArr = reshape(v2rfOverlap,[1 min(size(v2rfOverlap))*max(size(v2rfOverlap))]);
v2NoiseCorArr = reshape(v2NoiseCor,[1 min(size(v2NoiseCor))*max(size(v2NoiseCor))]);
[v2rfOverlapArr,sortOrder] = sort(v2rfOverlapArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

v2NoiseCorArr(v2rfOverlapArr==1) = []; v2rfOverlapArr(v2rfOverlapArr==1) = [];

figure(12); hold on; scatter(v2rfOverlapArr,v2NoiseCorArr,1,'filled','k');

expFit = fit(v2rfOverlapArr',v2NoiseCorArr','exp1');
v2expFit = plot(expFit,'predobs'); for i = 1:3, v2expFit(i).Color = [0, 0.4470, 0.7410]; v2expFit(i).LineWidth = 2; end
for i = 2:3, v2expFit(i).LineStyle = '--'; v2expFit(i).LineWidth = .75; end

legend('','Linear Fit', 'Exponential Fit');
title('V2 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j');

drawPublishAxis('labelFontSize=14');
leg = legend('','Exponential Fit', '95% Prediction bounds'); leg.Position = [0.6 0.2 0.2685 0.1003];

% third, do v3 %
v3rfOverlapArr = reshape(v3rfOverlap,[1 min(size(v3rfOverlap))*max(size(v3rfOverlap))]);
v3NoiseCorArr = reshape(v3NoiseCor,[1 min(size(v3NoiseCor))*max(size(v3NoiseCor))]);
[v3rfOverlapArr,sortOrder] = sort(v3rfOverlapArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

v3NoiseCorArr(v3rfOverlapArr==1) = []; v3rfOverlapArr(v3rfOverlapArr==1) = [];

figure(13); hold on; scatter(v3rfOverlapArr,v3NoiseCorArr,1,'filled','k');

expFit = fit(v3rfOverlapArr',v3NoiseCorArr','exp1');
v3expFit = plot(expFit,'predobs'); for i = 1:3, v3expFit(i).Color = [0, 0.4470, 0.7410]; v3expFit(i).LineWidth = 2; end
for i = 2:3, v3expFit(i).LineStyle = '--'; v3expFit(i).LineWidth = .75; end

legend('','Linear Fit', 'Exponential Fit');
title('V3 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j');

drawPublishAxis('labelFontSize=14');
leg = legend('','Exponential Fit', '95% Prediction bounds'); leg.Position = [0.6 0.2 0.2685 0.1003];


keyboard


