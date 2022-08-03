%% noiseStruct.m 
%
%       By: Josh Wilson
%       Created: March 2022
%
%       Takes two retinotopy scans (training and testing) and evaluates
%       different noise models. Calculates individual voxel variance and
%       covariances from training data. I port over the variance values of
%       the voxels for the testing data from the training data, scaled by the number of scans (but you can
%       fit if you want).
%
%       Covariance models I test on the training and testing data here:
%           - exact residaul covariance from the training data (should do horribly on test data but best on training)
%           - variance only
%           - global covariance/global covariance + channel covariance (both fit on the training data, from Jehee)
%           - rf Overlap and/or voxel distance exponential fits (alpha/beta/Nu parameters fit to training data)
%
%       Example Usage:
%
%           noiseStruct('trainData=s0401mc345GaussianHdrNM.mat','testData=s0401mc12GaussianHdrNM.mat','numTrain=3','numTest=2')
%           note: each scan should be the cleanRois and rois data from getprftseries. Set your voxel cutoffs there.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [train test] = noiseStruct(varargin)



%% Set up data, matrices: residual correlations, rf overlap, distance, variance 
% get the data
getArgs(varargin);
trainData = load(trainData);
testData = load(testData);

% structure data - here, put ALL the voxels from every ROI you want into
[xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords testpRFresidualtSeries] = structureData(trainData,testData);
numVoxels = length(xCenters);

% get noise Correlation, overlap, and distance matrices (n x n)
pRFresidualCorrMatrix = corr(pRFresidualtSeries);
rfOverlapMatrix = findOverlapMatrix(xCenters,yCenters,stdCenters);
distanceMatrix = findDistanceMatrix(scanCoords);
covarianceMatrix = std(pRFresidualtSeries)' * std(pRFresidualtSeries);

% plot it, just to make sure it's not fucked up
%figure,subplot(1,3,1);imshow(pRFresidualCorrMatrix);title('pRF residual correlations');...
%    subplot(1,3,2);imshow(rfOverlapMatrix);title('Receptive field Overlap');...
%    subplot(1,3,3);imshow(rescale(distanceMatrix)); title('3d Voxel distance');



%% calculate log likelihoods of different models %%
% get value to scale test covariance by number of scans
scaleByNumScans = sqrt(numTrainScans)/sqrt(numTestScans);

% set zero mean values for logmvnpdf input
pdfMean = zeros(1,numVoxels);


% Known covariance model %
    
    % calculate varCovar matrix using known variance and correlations
    trueVarCovarMatrix = makePosSemiDef( pRFresidualCorrMatrix.*covarianceMatrix, 1);
    
    % calculate likelihood with known variance-covariance matrix
    train.KnownVarCovarLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarCovarMatrix));
    test.KnownVarCovarLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, trueVarCovarMatrix*scaleByNumScans));
    
    
% Variance only model %

    % calculate only variance matrix from measured variance
    trueVarianceMatrix = makePosSemiDef( covarianceMatrix .* eye(numVoxels), 1);

    % calculate likelihood with known variance matrix
    train.KnownVarianceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarianceMatrix));
    test.KnownVarianceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, trueVarianceMatrix*scaleByNumScans));


% Fit single-term covariance Model %

    % fit covariance parameter Rho
    startparams(1) = 0; minsearch = [0]; maxsearch = [1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10,'UseParallel',1);
    
    [Rho, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeVarCovarLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix, 1);

    % calc likelihood with Jehee-inspired variance covariance matrix
    train.JeheeVarCovarLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix));
    test.JeheeVarCovarLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, jeheeVarCovarMatrix*scaleByNumScans));


% Variance, global covariance, and rf covariance model %

    % fit parameters %
    startparams = [0 0]; minsearch = [0 0]; maxsearch = [.5 1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10,'UseParallel',1);
    
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeFullLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning
    
    Rho = params(1); Sigma = params(2);

    % calc varCovar matrix as weighted sum of variance, covariance, and rf overlap
    jeheeFullMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix , 1);

    % calc likelihood with the covariance matrix
    train.JeheeFullLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix));
    test.JeheeFullLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, jeheeFullMatrix*scaleByNumScans));


% Individual covariance terms from overlap-corr fit model %

    % filter out overlap = 1 values to fit the exponential overlap/corfunction
    rfOverlapArray = rfOverlapMatrix(:); pRFresidualCorrArray = pRFresidualCorrMatrix(:);
    pRFresidualCorrArray = pRFresidualCorrArray(rfOverlapArray < 1); rfOverlapArray = rfOverlapArray(rfOverlapArray < 1);

    % fit the exponential
    expFit = fit(rfOverlapArray,pRFresidualCorrArray,'exp1');
    
    % make the matrix: get covar values from exponential fit, then 0 out the diagonal and add the variance matrix to get full varCovar
    rfOverlapCovarMatrix = arrayfun(@(x1)getCorFromFit(x1,expFit), rfOverlapMatrix);

    % adjust diagonal to 1 to keep variance
    rfOverlapCovarMatrix = rfOverlapCovarMatrix - rfOverlapCovarMatrix.*eye(numVoxels) + eye(numVoxels);

    % multiply by the covar matrix
    rfOverlapVarCovarMatrix = makePosSemiDef(rfOverlapCovarMatrix .* covarianceMatrix, 1);

    % calc likelihood with covariance matrix
    train.RfOverlapLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, rfOverlapVarCovarMatrix));
    test.RfOverlapLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, rfOverlapVarCovarMatrix*scaleByNumScans));


% Individual covariance terms from distance-corr fit model %

    % filter out overlap = 1 values to fit the exponential overlap/corfunction
    distanceArray = distanceMatrix(:); pRFresidualCorrArray = pRFresidualCorrMatrix(:);
    pRFresidualCorrArray = pRFresidualCorrArray(distanceArray > 0); distanceArray = distanceArray(distanceArray > 0);

    % fit the exponential
    expFit = fit(distanceArray,pRFresidualCorrArray,'exp1');
    
    % make the matrix: get covar values from exponential fit, then 0 out the diagonal and add the variance matrix to get full varCovar
    distanceCovarMatrix = arrayfun(@(x1)getCorFromFit(x1,expFit), distanceMatrix);

    % adjust diagonal to 1 to keep variance
    distanceCovarMatrix = distanceCovarMatrix - distanceCovarMatrix.*eye(numVoxels) + eye(numVoxels);

    % multiply by the covar matrix
    distanceVarCovarMatrix = makePosSemiDef(distanceCovarMatrix .* covarianceMatrix, 1);

    % calc likelihood with covariance matrix
    train.DistanceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, distanceVarCovarMatrix));
    test.DistanceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, distanceVarCovarMatrix*scaleByNumScans));


% Individual terms based on overlap and distance %

    % fit parameters %
    startparams = [.5]; minsearch = [0]; maxsearch = [1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10,'UseParallel',1);
    
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getOverlapDistanceLikelihood,startparams,minsearch,maxsearch,opts,pRFresidualtSeries,pdfMean,distanceVarCovarMatrix,rfOverlapVarCovarMatrix); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    Nu = params(1);

    % get the matrix with Nu %
    overlapDistanceVarCovarMatrix = makePosSemiDef( Nu*rfOverlapVarCovarMatrix + (1-Nu)*distanceVarCovarMatrix, 1);

    % calc likelihood with covariance matrix
    train.OverlapDistanceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix))
    test.OverlapDistanceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix*scaleByNumScans))







keyboard










%% END OF SCRIPT %%








%%%%%%%%%%%%%%%%%%%
%% structureData %%
%%%%%%%%%%%%%%%%%%%

function [xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords testpRFresidualtSeries] = structureData(trainData, testData);

% make empty arrays/matrices of what you want
xCenters = []; yCenters = []; stdCenters = []; tSeries = [];  pRFtSeries  = []; pRFresidualtSeries = []; scanCoords = []; testpRFresidualtSeries = [];

%loop through the ROIs you want to include
for roi = 3

    % Training data: add things to the centralized matrices
    xCenters = [xCenters trainData.cleanRois(roi).vox.x];
    yCenters = [yCenters trainData.cleanRois(roi).vox.y];
    stdCenters = [stdCenters trainData.cleanRois(roi).vox.rfstd];
    tSeries = [tSeries trainData.cleanRois(roi).vox.tSeries] ;
    pRFtSeries = [pRFtSeries trainData.cleanRois(roi).vox.pRFtSeries];
    pRFresidualtSeries = [pRFresidualtSeries trainData.cleanRois(roi).vox.tSeries - trainData.cleanRois(roi).vox.pRFtSeries];
    scanCoords = [scanCoords trainData.cleanRois(roi).vox.scanCoords];

    % Testing data: get the time series
    testpRFresidualtSeries = [testpRFresidualtSeries testData.rois(roi).vox.baselineNoise(:,ismember(testData.rois(roi).vox.linearCoords, trainData.cleanRois(roi).vox.linearCoords)) ];

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

function outputMatrix = makePosSemiDef(inputMatrix, checkWarning);

% initial warning flag that will trigger if your ending matrix is adjusted by a large amount
warningFlag = 0; smallWarning = 0;

% check if minimum eigenvalue is below zero, else end
if min(eig(inputMatrix)) > 0;
    
    outputMatrix = inputMatrix;

% if so, do adjustment
else
    
    % offer a warning if the adjustment identity matrix is within 5 OOM of minimum covar value
    eigToAdd = min(eig(inputMatrix));
    
    if min(abs(inputMatrix(:)))/eigToAdd > -10^5;
        %disp(sprintf(['*** Identity matrix you are adding is less than 5 orders of magnitude less than the minimum covariance value. \n' ...
        %    'Consider re-checking whether your matrix is non positive semi-definite due to rounding error, or actual contents. ***']))
        %disp(sprintf('Adding %i. \nMinimum covariance value: %i. \nMedian covariance value: %i.',abs(eigToAdd),min(min(abs(inputMatrix(:)))),median(abs(inputMatrix), 'all')))
        
        % trigger warning flag
        if checkWarning; warningFlag = 1; end

    end

    % add an identify matrix multiplied by the absolute value of the minimum eigenvalue of the input matrix to make all eigenvalues positive
    outputMatrix = inputMatrix+eye(length(inputMatrix))*abs(min(eig(inputMatrix)));
    if checkWarning; smallWarning = 1; end

end


% sometimes the adjustment is too small for matab to make and the output matrix will still be non positive. 
% If so, add a very small value (I/10^18) (I've found this usually works)
if min(eig(outputMatrix)) < 0;

    % offer a warning if the adjustment identity matrix is within 5 OOM of minimum covar value
    %disp(sprintf('Adding Identity matrix of minimum eigenvalue didnt make matrix positive semi definite. \n Adding a very small Identity value (1/10^18) instead.'))

    if min(abs(outputMatrix(:)))/(1/10^18) < 10^5;

         % trigger warning flag
        if checkWarning; warningFlag = 1; end

    end

    % add the small identity matrix
    outputMatrix = outputMatrix + eye(length(outputMatrix))/10^18;
    if checkWarning; smallWarning = 1; end
end

% send a warning if final output was adjusted by a large amount
if warningFlag
sprintf('!!!!! WARNING !!!!! \nYour final covariance matrix had to be adjusted by a non-trivial amount in order to make it Positive Semi Definite. \nIt is likely that the reason it was not PSD was due to something other than a rounding error.')
end

if smallWarning
sprintf('!!!!! WARNING !!!!! You added a small value to make the var/covar matrix positive definite.')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeVarCovarLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeVarCovarLikelihood(Rho,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean)

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix , 0);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeFullLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeFullLikelihood(params,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix)
    
    Rho = params(1); Sigma = params(2);

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeFullMatrix = makePosSemiDef( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix , 0);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix));



%%%%%%%%%%%%%%%%%%%%%%%
%% getCorFromOverlap %%
%%%%%%%%%%%%%%%%%%%%%%%

function y = getCorFromFit(x,expFit)

% output of each matrix entry is the exponential fit to the data
y = expFit.a*exp(x*expFit.b);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getOverlapDistanceLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getOverlapDistanceLikelihood(params,pRFresidualtSeries,pdfMean,distanceVarCovarMatrix,rfOverlapVarCovarMatrix)
    
    Nu = params(1);

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    overlapDistanceVarCovarMatrix = makePosSemiDef( Nu*rfOverlapVarCovarMatrix + (1-Nu)*distanceVarCovarMatrix, 0);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix));
