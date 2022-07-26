%% noiseStruct.m 
%
%       By: Josh Wilson
%       Created: March 2022
%
%       I don't know what this does yet.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function noiseStruct(varargin)

%% Set up data, matrices: residual correlations, rf overlap, distance, variance 
% get arguments
getArgs(varargin);
load(data);

% structure data - here, put ALL the voxels from every ROI you want into
[xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords] = structureData(cleanRois);
numVoxels = length(xCenters);

% get noise Correlation, overlap, and distance matrices (n x n)
pRFresidualCorrMatrix = corr(pRFresidualtSeries);
rfOverlapMatrix = findOverlapMatrix(xCenters,yCenters,stdCenters);
distanceMatrix = findDistanceMatrix(scanCoords);
varianceMatrix = std(pRFresidualtSeries)' * std(pRFresidualtSeries);

% plot it, just to make sure it's not fucked up
figure,subplot(1,3,1);imshow(pRFresidualCorrMatrix);title('pRF residual correlations');...
    subplot(1,3,2);imshow(rfOverlapMatrix);title('Receptive field Overlap');...
    subplot(1,3,3);imshow(rescale(distanceMatrix)); title('3d Voxel distance');



%% do stuff%%

manualVarCovarMatrix = pRFresidualCorrMatrix.*varianceMatrix;
matlabVarCovarMatrix = cov(pRFresidualtSeries);

keyboard


logmvnpdf(pRFresidualtSeries(1,:),zeros(1,numVoxels),makePosSemiDef(manualVarCovarMatrix))

logmvnpdf(pRFresidualtSeries(1,:),zeros(1,numVoxels),makePosSemiDef(matlabVarCovarMatrix))


keyboard













%% END OF SCRIPT %%








%%%%%%%%%%%%%%%%%%%
%% structureData %%
%%%%%%%%%%%%%%%%%%%

function [xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords] = structureData(cleanRois);

% make empty arrays/matrices of what you want
xCenters = []; yCenters = []; stdCenters = []; tSeries = [];  pRFtSeries  = []; pRFresidualtSeries = []; scanCoords = [];

%loop through the ROIs you want to include
for roi = 1:1

    % add things to the centralized matrices
    xCenters = [xCenters cleanRois(roi).vox.x];
    yCenters = [yCenters cleanRois(roi).vox.y];
    stdCenters = [stdCenters cleanRois(roi).vox.rfstd];
    tSeries = [tSeries cleanRois(roi).vox.tSeries];
    pRFtSeries = [pRFtSeries cleanRois(roi).vox.pRFtSeries];
    pRFresidualtSeries = [pRFresidualtSeries cleanRois(roi).vox.tSeries - cleanRois(roi).vox.pRFtSeries];
    scanCoords = [scanCoords cleanRois(roi).vox.scanCoords];

end



%%%%%%%%%%%%%%%%%%%%%%%
%% findOverlapMatrix %%
%%%%%%%%%%%%%%%%%%%%%%%

function rfOverlapMatrix = findOverlapMatrix(xCenters, yCenters, stdCenters)

% pre-allocate
rfOverlapMatrix = zeros(length(xCenters));

for row = 1:length(xCenters)
    for column = 1:length(xCenters)   
        mu1 = 0; s1 = stdCenters(column); s2 = stdCenters(row);
        mu2 = sqrt((xCenters(column)-xCenters(row))^2 + (yCenters(column)-yCenters(row))^2);
        if mu1 == mu2 & s1 == s2; rfOverlapMatrix(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        rfOverlapMatrix(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%
%% findDistanceMatrix %%
%%%%%%%%%%%%%%%%%%%%%%%%

function distanceMatrix = findDistanceMatrix(scanCoords)

% pre-allocate
distanceMatrix = zeros(length(scanCoords));

% calculate pythagorean distance between all voxels
for row = 1:length(scanCoords)
    for column = 1:length(scanCoords)
        dist = sqrt((scanCoords(1,column)-scanCoords(1,row))^2 + (scanCoords(2,column)-scanCoords(2,row))^2 + (scanCoords(3,column)-scanCoords(3,row))^2);
        distanceMatrix(row,column) = dist;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makePositiveSemiDefinite %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputMatrix = makePosSemiDef(inputMatrix);

% check if minimum eigenvalue is below zero, else end
if min(eig(inputMatrix)) > 0;
    
    outputMatrix = inputMatrix;

% if not, do adjustment
else
    
    % check the adjustment and issue warning if it is large
    eigToAdd = min(eig(inputMatrix));
    if median(abs(inputMatrix(:)))/eigToAdd > -10^10;
        disp(sprintf(['Identity matrix you are adding is less than 10 orders of magnitude less than the median covariance value. \n' ...
            'Consider re-checking whether your matrix is non positive semi-definite due to rounding error, or actual contents.']))
    end

    % add an identify matrix multiplied by the minimum eigenvalue of the input matrix
    outputMatrix = inputMatrix+eye(length(inputMatrix))*abs(min(eig(inputMatrix)));

end







