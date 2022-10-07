%% shuffleCor.m %%
%
%       By: Josh Wilson
%       Created: March 2022
%
%       Takes 2 scans and calculates the residual correlations between them. Used to show that voxels with overlapping receptive fields
%       have correlated noise between different scans, where noise is completely independent. 
%       Thus, any correlations arise from model failure (both receptive field and hemodynamic response).
%
%       Needs 3 scan inputs: base is the original pRF on the average of all the scans, and scan 1 and scan 2 are however you divide of the
%       scans (see example usage).
%
%       Example usage:
%           shufflecor('data1=s0401mc12GaussianHdrNM','data2=s0401mc345GaussianHdrNM','base=s0401pRF.mat')       
%           note: each scan should be the cleanRois and rois data from getprftseries. Set your voxel cutoffs there.

 

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
roisA = rois;
load(data2);
roisB = rois;
load(base,'rois','cleanRois');
keyboard



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
end



keyboard




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
        if mu1 == mu2 & s1 == s2; v1rfOverlap(roi1vox,roi2vox) = 1;
        elseif s1 < 0 || s2 < 0; v1rfOverlap(roi1vox,roi2vox) = 1; else %takes out voxels with negative std fits because we filter out the full overlaps later
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v1rfOverlap = transpose(v1Overlap).*100;

%v2 
for roi1vox = 1:length(cleanRoisA(2).vox.linearCoords);
    for roi2vox = 1:length(cleanRoisB(2).vox.linearCoords);       
        mu1 = 0; s1 = cleanRoisA(2).vox.rfstd(roi1vox); s2 = cleanRoisB(2).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRoisA(2).vox.x(roi1vox)-cleanRoisB(2).vox.x(roi2vox))^2 + (cleanRoisA(2).vox.y(roi1vox)-cleanRoisB(2).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v2rfOverlap(roi1vox,roi2vox) = 1;
        elseif s1 < 0 || s2 < 0; v2rfOverlap(roi1vox,roi2vox) = 1; else
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v2Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v2rfOverlap = transpose(v2Overlap).*100;

%v3
for roi1vox = 1:length(cleanRoisA(3).vox.linearCoords);
    for roi2vox = 1:length(cleanRoisB(3).vox.linearCoords);       
        mu1 = 0; s1 = cleanRoisA(3).vox.rfstd(roi1vox); s2 = cleanRoisB(3).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRoisA(3).vox.x(roi1vox)-cleanRoisB(3).vox.x(roi2vox))^2 + (cleanRoisA(3).vox.y(roi1vox)-cleanRoisB(3).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;
        elseif s1 < 0 || s2 < 0; v3rfOverlap(roi1vox,roi2vox) = 1; else
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v3Overlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
v3rfOverlap = transpose(v3Overlap).*100;



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
v2rfOverlapArr = reshape(v2rfOverlap,[1 min(size(v2rfOverlap))*max(size(v2rfOverlap))]);
v2NoiseCorArr = reshape(v2NoiseCor,[1 min(size(v2NoiseCor))*max(size(v2NoiseCor))]);
[v2rfOverlapArr,sortOrder] = sort(v2rfOverlapArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

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
v3rfOverlapArr = reshape(v3rfOverlap,[1 min(size(v3rfOverlap))*max(size(v3rfOverlap))]);
v3NoiseCorArr = reshape(v3NoiseCor,[1 min(size(v3NoiseCor))*max(size(v3NoiseCor))]);
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



%%%%%%%%%%%%%%%%%%%
%% END OF SCRIPT %%
%%%%%%%%%%%%%%%%%%%














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

