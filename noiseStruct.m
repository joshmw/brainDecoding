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
covarianceMatrix = std(pRFresidualtSeries)' * std(pRFresidualtSeries);

% plot it, just to make sure it's not fucked up
figure,subplot(1,3,1);imshow(pRFresidualCorrMatrix);title('pRF residual correlations');...
    subplot(1,3,2);imshow(rfOverlapMatrix);title('Receptive field Overlap');...
    subplot(1,3,3);imshow(rescale(distanceMatrix)); title('3d Voxel distance');



%% calculate log likelihoods of different models %%
% set mean for mvnPdf as zero
pdfMean = zeros(1,numVoxels);

% Known covariance model %
    
    % calculate varCovar matrix using known variance and correlations
    trueVarCovarMatrix = makePosSemiDef( pRFresidualCorrMatrix.*covarianceMatrix );
    
    % calculate likelihood with known variance-covariance matrix
    knownCovarianceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarCovarMatrix))
    
    
% Variance only model %

    % calculate only variance matrix from measured variance
    trueVarianceMatrix = makePosSemiDef( covarianceMatrix .* eye(numVoxels) );

    % calculate likelihood with known variance matrix
    knownVarianceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarianceMatrix))


% Fit single-term covariance Model %

    % fit covariance parameter Rho
    startparams(1) = 0; minsearch = [0]; maxsearch = [1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10);
    
    [Rho, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeVarCovarLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix );

    % calc likelihood with Jehee-inspired variance covariance matrix
    jeheeVarCovarLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix))


% Variance, global covariance, and rf covariance %

    % fit parameters %
    startparams = [0 0]; minsearch = [0 0]; maxsearch = [1 1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10);
    
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeFullLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning
    
    Rho = params(1); Sigma = params(2);

    % calc varCovar matrix as weighted sum of variance, covariance, and rf overlap
    jeheeFullMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix );

    % calc likelihood with the covariance matrix
    jeheeFullLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix))


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

% loop through each receptive field (row x column) and caluculate overlap as 2d gaussian intersection
for row = 1:length(xCenters)
    for column = 1:length(xCenters)   
        
        % draw rfs as 2d gaussians
        mu1 = 0; s1 = stdCenters(column); s2 = stdCenters(row);
        mu2 = sqrt((xCenters(column)-xCenters(row))^2 + (yCenters(column)-yCenters(row))^2);
        
        % calculate overlap
        if mu1 == mu2 & s1 == s2; rfOverlapMatrix(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        rfOverlapMatrix(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    
    end %column of overlap matrix

end %row of overlap matrix



%%%%%%%%%%%%%%%%%%%%%%%%
%% findDistanceMatrix %%
%%%%%%%%%%%%%%%%%%%%%%%%

function distanceMatrix = findDistanceMatrix(scanCoords)

% pre-allocate
distanceMatrix = zeros(length(scanCoords));

% calculate pythagorean distance between all voxels
for row = 1:length(scanCoords)
    for column = 1:length(scanCoords)

        %calculate dist as the sqrt of squared sum of xyz distances
        dist = sqrt((scanCoords(1,column)-scanCoords(1,row))^2 + (scanCoords(2,column)-scanCoords(2,row))^2 + (scanCoords(3,column)-scanCoords(3,row))^2);
        distanceMatrix(row,column) = dist;

    end %column of distance matrix

end %row of distance matrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makePositiveSemiDefinite %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makes a matrix positive semi definite by adding the minimum eigenvalue. This occurs in matlab with very small values due to rounding errors.
% If the minmum eigenvalue is negaative enough to have potentially arisen from something other than rounding error, it will offer you a warning.
% Theoretically, all covariance matrices should be positive semi definite so this correction should only be fixing rounding errors if you pass in covar matrices.

function outputMatrix = makePosSemiDef(inputMatrix);

% check if minimum eigenvalue is below zero, else end
if min(eig(inputMatrix)) > 0;
    
    outputMatrix = inputMatrix;

% if so, do adjustment
else
    
    % offer a warning if the adjustment identity matrix is within 5 OOM of minimum covar value
    eigToAdd = min(eig(inputMatrix));
    if min(abs(inputMatrix(:)))/eigToAdd > -10^5;
        disp(sprintf(['*** Identity matrix you are adding is less than 5 orders of magnitude less than the minimum covariance value. \n' ...
            'Consider re-checking whether your matrix is non positive semi-definite due to rounding error, or actual contents. ***']))

    end

    % add an identify matrix multiplied by the absolute value of the minimum eigenvalue of the input matrix to make all eigenvalues positive
    outputMatrix = inputMatrix+eye(length(inputMatrix))*abs(min(eig(inputMatrix)));

end


% sometimes the adjustment is too small for matab to make and the output matrix will still be non positive. 
% If so, add a very small value (I/10^18) (I've found this usually works)
if min(eig(outputMatrix)) < 0;

    % offer a warning if the adjustment identity matrix is within 5 OOM of minimum covar value
    disp(sprintf('Adding Identity matrix of minimum eigenvalue didnt make matrix positive semi definite. \n Adding a very small Identity value (1/10^18) instead.'))

    if min(abs(outputMatrix(:)))/(1/10^18) < 10^5;
        disp(sprintf(['*** Identity matrix you are adding is less than 5 orders of magnitude less than the minimum covariance value. \n' ...
            'Consider re-checking whether your matrix is non positive semi-definite due to rounding error, or actual contents. ***']))
    end

    % add the small identity matrix
    outputMatrix = outputMatrix + eye(length(outputMatrix))/10^18;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeVarCovarLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeVarCovarLikelihood(Rho,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean)

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix );

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeFullLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeFullLikelihood(params,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix)
    
    Rho = params(1); Sigma = params(2);

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeFullMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix );

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix));






