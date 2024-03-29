%% noiseStruct.m 
%
%       By: Josh Wilson
%       Created: july 2022
%
%       Takes two retinotopy scans (training and testing) and evaluates different noise models. Calculates individual voxel variance and
%       covariances from training data. I port over the variance values of the voxels for the testing data from the training data, scaled 
%       by the number of scans (but you can fit if you want).
%
%       Covariance models I test on the training and testing data here:
%           - exact residaul covariance from the training data (should do horribly on test data but best on training)
%           - variance only
%           - global covariance/global covariance + channel covariance (both fit on the training data, from Jehee)
%           - rf Overlap and/or voxel distance exponential fits (alpha/beta/Nu parameters fit to training data)
%
%       Example Usage:
%
%           noiseStruct(training data, testing data, number of scans averaged in training data, number of scans averaged in testing data, rois you want to use)
%           noiseStruct('trainData=s0401mc345GaussianHdrNM.mat','testData=s0401mc12GaussianHdrNM.mat','numTrain=3','numTest=2','rois=[1 2 3]')
%           note: each scan should be the cleanRois and rois data from getprftseries. Set your voxel cutoffs and list your rois there.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [train test] = noiseStruct(varargin)



%% Set up data, matrices: residual correlations, rf overlap, distance, variance 
% get the data
getArgs(varargin); subject = trainData(1:5);
trainData = load(trainData);
testData = load(testData);

% structure data - here, put ALL the voxels from every ROI you want into
[xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords testpRFresidualtSeries] = structureData(trainData,testData,rois);
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
% get value to scale test covariance by number of scans
scaleByNumScans = sqrt(numTrainScans)/sqrt(numTestScans)

% set zero mean values for logmvnpdf input
pdfMean = zeros(1,numVoxels);


% Known covariance model %
    
    % calculate varCovar matrix using known variance and correlations
    trueVarCovarMatrix = nearestSPD( pRFresidualCorrMatrix.*covarianceMatrix);
    minEig.KnwonVarCovarLikelihood = min( eig( pRFresidualCorrMatrix.*covarianceMatrix ));

    % calculate likelihood with known variance-covariance matrix
    train.KnownVarCovarLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarCovarMatrix));
    test.KnownVarCovarLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, trueVarCovarMatrix*scaleByNumScans));
    
    
% Variance only model %

    % calculate only variance matrix from measured variance
    trueVarianceMatrix = nearestSPD( covarianceMatrix .* eye(numVoxels));
    minEig.trueVarianceMatrix = min( eig (covarianceMatrix .* eye(numVoxels)));

    % calculate likelihood with known variance matrix
    train.KnownVarianceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, trueVarianceMatrix));
    test.KnownVarianceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, trueVarianceMatrix*scaleByNumScans));


% Fit single-term covariance Model %

    % fit covariance parameter Rho
    startparams(1) = 0; minsearch = [0]; maxsearch = [1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10,'UseParallel',1);
    
    [singleTermRho, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeVarCovarLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = nearestSPD( (1-singleTermRho)*covarianceMatrix.*eye(numVoxels) + singleTermRho*covarianceMatrix);
    minEig.jeheeVarCovarMatrix = min( eig( (1-singleTermRho)*covarianceMatrix.*eye(numVoxels) + singleTermRho*covarianceMatrix));

    % calc likelihood with Jehee-inspired variance covariance matrix
    train.JeheeVarCovarLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix));
    test.JeheeVarCovarLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, jeheeVarCovarMatrix*scaleByNumScans));


% Variance, global covariance, and rf covariance model %

    % fit parameters %
    startparams = [.2 .2]; minsearch = [0 0]; maxsearch = [1 10];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-10,'TolFun',10^-10,'UseParallel',1);
    
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getJeheeFullLikelihood,startparams,minsearch,maxsearch,opts,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning
    
    fullRho = params(1); Sigma = params(2)^2;

    % calc varCovar matrix as weighted sum of variance, covariance, and rf overlap
    jeheeFullMatrix = nearestSPD( (1-fullRho-Sigma)*covarianceMatrix.*eye(numVoxels) + fullRho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix );
    minEig.jeheeFullMatrix = min( eig ((1-fullRho-Sigma)*covarianceMatrix.*eye(numVoxels) + fullRho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix));

    % calc likelihood with the covariance matrix
    train.JeheeFullLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix));
    test.JeheeFullLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, jeheeFullMatrix*scaleByNumScans));


% Individual covariance terms from overlap-corr fit model %

    % filter out overlap = 1 values to fit the exponential overlap/corfunction
    rfOverlapArray = rfOverlapMatrix(:); pRFresidualCorrArray = pRFresidualCorrMatrix(:);
    pRFresidualCorrArray = pRFresidualCorrArray(rfOverlapArray < 1); rfOverlapArray = rfOverlapArray(rfOverlapArray < 1);

    % fit the exponential
    curveFitOpts = fitoptions('exp2','lower',[-5 -5 -5 -5],'upper',[5 5 5 5]);
    overlapCurveFit = fit(rfOverlapArray,pRFresidualCorrArray,'exp2',curveFitOpts);
    showFit(rfOverlapArray,pRFresidualCorrArray,overlapCurveFit,1);
    
    % make the matrix: get covar values from exponential fit, then 0 out the diagonal and add the variance matrix to get full varCovar
    rfOverlapCovarMatrix = arrayfun(@(x1)getCorFromFit(x1,overlapCurveFit), rfOverlapMatrix);

    % adjust diagonal to 1 to keep variance
    rfOverlapCovarMatrix = nearestSPD(rfOverlapCovarMatrix - rfOverlapCovarMatrix.*eye(numVoxels) + eye(numVoxels));

    % multiply by the covar matrix
    rfOverlapVarCovarMatrix = nearestSPD(rfOverlapCovarMatrix .* covarianceMatrix);
    minEig.rfOverlapVarCovarMatrix = min( eig( rfOverlapCovarMatrix .* covarianceMatrix));

    % calc likelihood with covariance matrix
    train.RfOverlapLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, rfOverlapVarCovarMatrix));
    test.RfOverlapLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, rfOverlapVarCovarMatrix*scaleByNumScans));


% Individual covariance terms from distance-corr fit model %

    % filter out overlap = 1 values to fit the exponential overlap/corfunction
    distanceArray = distanceMatrix(:); pRFresidualCorrArray = pRFresidualCorrMatrix(:);
    pRFresidualCorrArray = pRFresidualCorrArray(distanceArray > 0); distanceArray = distanceArray(distanceArray > 0);

    % fit the exponential
    curveFitOpts = fitoptions('exp2','lower',[-5 -5 -5 -5],'upper',[5 5 5 5]);
    distanceCurveFit = fit(distanceArray,pRFresidualCorrArray,'exp2');
    showFit(distanceArray,pRFresidualCorrArray,distanceCurveFit,2);

    % make the matrix: get covar values from exponential fit, then 0 out the diagonal and add the variance matrix to get full varCovar
    distanceCovarMatrix = arrayfun(@(x1)getCorFromFit(x1,distanceCurveFit), distanceMatrix);

    % adjust diagonal to 1 to keep variance
    distanceCovarMatrix = nearestSPD(distanceCovarMatrix - distanceCovarMatrix.*eye(numVoxels) + eye(numVoxels));

    % multiply by the covar matrix
    distanceVarCovarMatrix = nearestSPD(distanceCovarMatrix .* covarianceMatrix);
    minEig.distanceVarCovarMatrix = min( eig( distanceCovarMatrix .* covarianceMatrix));

    % calc likelihood with covariance matrix
    train.DistanceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, distanceVarCovarMatrix));
    test.DistanceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, distanceVarCovarMatrix*scaleByNumScans));


% Individual terms based on overlap and distance %

    % fit parameters %
    startparams = [0]; minsearch = [0]; maxsearch = [1];
    opts = optimset('display','off','maxIter',100000,'MaxFunEvals',100000,'TolX',10^-30,'TolFun',10^-20,'UseParallel',1);
    
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = ... 
        lsqnonlin(@getOverlapDistanceLikelihood,startparams,minsearch,maxsearch,opts,pRFresidualtSeries,pdfMean,distanceVarCovarMatrix,rfOverlapVarCovarMatrix); 
         [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning

    Nu = params(1);

    % get the matrix with Nu %
    overlapDistanceVarCovarMatrix = nearestSPD( Nu*rfOverlapVarCovarMatrix + (1-Nu)*distanceVarCovarMatrix);
    minEig.overlapDistanceVarCovarMatrix = min( eig(  Nu*rfOverlapVarCovarMatrix + (1-Nu)*distanceVarCovarMatrix));

    % calc likelihood with covariance matrix
    train.OverlapDistanceLikelihood = sum( logmvnpdf( pRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix))
    test.OverlapDistanceLikelihood = sum( logmvnpdf( testpRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix*scaleByNumScans))


%% Show the matrices %%
figure;
subplot(3,3,5); imshow(imshowScale(corr(testpRFresidualtSeries))); title('True residual var/covar matrix');
subplot(3,3,1); imshow(imshowScale(trueVarCovarMatrix./covarianceMatrix)); title('First scan true var/covar matrix');
subplot(3,3,2); imshow(imshowScale(eye(numVoxels))); title('Variance Only');
subplot(3,3,3); imshow(imshowScale(jeheeVarCovarMatrix./covarianceMatrix)); title('Global Covariance');
subplot(3,3,4); imshow(imshowScale(jeheeFullMatrix./covarianceMatrix)); title('Jehee Matrix');
subplot(3,3,6), imshow(imshowScale(rfOverlapVarCovarMatrix./covarianceMatrix)); title('RF Overlap Fit Matrix');
subplot(3,3,7); imshow(imshowScale(distanceVarCovarMatrix./covarianceMatrix)); title('distance Matrix');
subplot(3,3,8); imshow(imshowScale(overlapDistanceVarCovarMatrix./covarianceMatrix)); title('Distance/Overlap fit matrix');
set(gcf,'Position',[1442 -161 1328 958]);
sgtitle(strcat(subject,' correlation matrices'));



%% Calculate AIC and BIC %%
    numObservations = size(pRFresidualtSeries,1);

    % calculate number of parameters in each model %
    numParams.KnownVarCovar = numVoxels*numVoxels;
    numParams.KnownVariance = numVoxels;
    numParams.JeheeVarCovar = numVoxels + 1; %voxels + rho
    numParams.JeheeFull = numVoxels + 2; %voxels + rho + sigma
    numParams.RfOverlap = numVoxels + length(coeffvalues(overlapCurveFit)); % variances + curve fit
    numParams.Distance = numVoxels + length(coeffvalues(distanceCurveFit)); % variances + curve fit
    numParams.OverlapDistance = numVoxels + length(coeffvalues(overlapCurveFit)) + length(coeffvalues(distanceCurveFit)) + 1; % variances, both curves, Nu

    % Calculate AIC %
        
        % training AIC %
        aicTrain.KnownVarCovar = 2*numParams.KnownVarCovar - 2*train.KnownVarCovarLikelihood;
        aicTrain.KnownVariance =  2*numParams.KnownVariance - 2*train.KnownVarianceLikelihood;
        aicTrain.JeheeVarCovar = 2*numParams.JeheeVarCovar - 2*train.JeheeVarCovarLikelihood;
        aicTrain.JeheeFull = 2*numParams.JeheeFull - 2*train.JeheeFullLikelihood;
        aicTrain.RfOverlap = 2*numParams.RfOverlap - 2*train.RfOverlapLikelihood;
        aicTrain.Distance = 2*numParams.Distance - 2*train.DistanceLikelihood;
        aicTrain.OverlapDistance = 2*numParams.OverlapDistance - 2*train.OverlapDistanceLikelihood;

        % testing AIC %
        aicTest.KnownVarCovar = 2*numParams.KnownVarCovar - 2*test.KnownVarCovarLikelihood;
        aicTest.KnownVariance =  2*numParams.KnownVariance - 2*test.KnownVarianceLikelihood;
        aicTest.JeheeVarCovar = 2*numParams.JeheeVarCovar - 2*test.JeheeVarCovarLikelihood;
        aicTest.JeheeFull = 2*numParams.JeheeFull - 2*test.JeheeFullLikelihood;
        aicTest.RfOverlap = 2*numParams.RfOverlap - 2*test.RfOverlapLikelihood;
        aicTest.Distance = 2*numParams.Distance - 2*test.DistanceLikelihood;
        aicTest.OverlapDistance = 2*numParams.OverlapDistance - 2*test.OverlapDistanceLikelihood;

    % Calculate BIC %
        
        % training BIC %
        bicTrain.KnownVarCovar = log(numObservations)*numParams.KnownVarCovar - 2*train.KnownVarCovarLikelihood;
        bicTrain.KnownVariance =  log(numObservations)*numParams.KnownVariance - 2*train.KnownVarianceLikelihood;
        bicTrain.JeheeVarCovar = log(numObservations)*numParams.JeheeVarCovar - 2*train.JeheeVarCovarLikelihood;
        bicTrain.JeheeFull = log(numObservations)*numParams.JeheeFull - 2*train.JeheeFullLikelihood;
        bicTrain.RfOverlap = log(numObservations)*numParams.RfOverlap - 2*train.RfOverlapLikelihood;
        bicTrain.Distance = log(numObservations)*numParams.Distance - 2*train.DistanceLikelihood;
        bicTrain.OverlapDistance = log(numObservations)*numParams.OverlapDistance - 2*train.OverlapDistanceLikelihood;

        % testing BIC %
        bicTest.KnownVarCovar = log(numObservations)*numParams.KnownVarCovar - 2*test.KnownVarCovarLikelihood;
        bicTest.KnownVariance =  log(numObservations)*numParams.KnownVariance - 2*test.KnownVarianceLikelihood;
        bicTest.JeheeVarCovar = log(numObservations)*numParams.JeheeVarCovar - 2*test.JeheeVarCovarLikelihood;
        bicTest.JeheeFull = log(numObservations)*numParams.JeheeFull - 2*test.JeheeFullLikelihood;
        bicTest.RfOverlap = log(numObservations)*numParams.RfOverlap - 2*test.RfOverlapLikelihood;
        bicTest.Distance = log(numObservations)*numParams.Distance - 2*test.DistanceLikelihood;
        bicTest.OverlapDistance = log(numObservations)*numParams.OverlapDistance - 2*test.OverlapDistanceLikelihood;




%% Plot the results %%

    % training data %
        figure; subplot(1,3,1);
        
        % rename stuff and graph
        names = categorical({'True var/covar','Variance Only','Variance + global covariance','Full Jehee Model','RF Overlap Model','Distance Model','Combined Distance RF Model'});
        names = reordercats(names,{'True var/covar','Variance Only','Variance + global covariance','Full Jehee Model','RF Overlap Model','Distance Model','Combined Distance RF Model'});
        bar(names, ...
            [aicTrain.KnownVarCovar bicTrain.KnownVarCovar; aicTrain.KnownVariance bicTrain.KnownVariance; aicTrain.JeheeVarCovar bicTrain.JeheeVarCovar; ...
             aicTrain.JeheeFull bicTrain.JeheeFull; aicTrain.RfOverlap bicTrain.RfOverlap; aicTrain.Distance bicTrain.Distance; aicTrain.OverlapDistance aicTrain.OverlapDistance])
        
        % label
        title(strcat(subject,' Training data')),xlabel('Model'),ylabel('Average loglikelihood of one observation'); legend('AIC','BIC');
        
    % testing data %
        subplot(1,3,2)
    
        % rename stuff and graph
        names = categorical({'Variance Only','Variance + global covariance','Full Jehee Model','RF Overlap Model','Distance Model','Combined Distance RF Model'});
        names = reordercats(names,{'Variance Only','Variance + global covariance','Full Jehee Model','RF Overlap Model','Distance Model','Combined Distance RF Model'});
        bar(names, ...
            [aicTest.KnownVariance bicTest.KnownVariance; aicTest.JeheeVarCovar bicTest.JeheeVarCovar; ...
             aicTest.JeheeFull bicTest.JeheeFull; aicTest.RfOverlap bicTest.RfOverlap; aicTest.Distance bicTest.Distance; aicTest.OverlapDistance aicTest.OverlapDistance])
        
        % label
        title(strcat(subject,' Testing data')),xlabel('Model'),ylabel('Average loglikelihood of one observation');
    
    % plot the parameters %
        subplot(1,3,3)

        %rename, graph, title
        names = categorical({'Single Term Rho','Full Jehee Rho','Full Jehee Sigma','Nu'}); names = reordercats(names,{'Single Term Rho','Full Jehee Rho','Full Jehee Sigma','Nu'});
        bar(names, [singleTermRho fullRho Sigma Nu])
        ylim([0 1]),title('parameter estimates')
        


keyboard









% END OF SCRIPT %



%~%~%~%~%~%~%~%~%~%~%~
%% HELPER FUNCTIONS %%
%~%~%~%~%~%~%~%~%~%~%~



%%%%%%%%%%%%%%%%%%%
%% structureData %%
%%%%%%%%%%%%%%%%%%%

function [xCenters yCenters stdCenters tSeries pRFtSeries pRFresidualtSeries scanCoords testpRFresidualtSeries] = structureData(trainData, testData, rois);

% make empty arrays/matrices of what you want
xCenters = []; yCenters = []; stdCenters = []; tSeries = [];  pRFtSeries  = []; pRFresidualtSeries = []; scanCoords = []; testpRFresidualtSeries = [];

%loop through the ROIs you want to include
for roi = rois

    % Training data: add things to the centralized matrices
    xCenters = [xCenters trainData.cleanRois(roi).vox.x];
    yCenters = [yCenters trainData.cleanRois(roi).vox.y];
    stdCenters = [stdCenters trainData.cleanRois(roi).vox.rfstd];
    tSeries = [tSeries trainData.cleanRois(roi).vox.tSeries] ;
    pRFtSeries = [pRFtSeries trainData.cleanRois(roi).vox.pRFtSeries];
    pRFresidualtSeries = [pRFresidualtSeries trainData.cleanRois(roi).vox.tSeries - trainData.cleanRois(roi).vox.pRFtSeries];
    scanCoords = [scanCoords trainData.cleanRois(roi).vox.scanCoords];

    % Testing data: get the time series
    testpRFresidualtSeries = [testpRFresidualtSeries testData.rois(roi).vox.baselineNoise(:,ismember(trainData.rois(roi).vox.linearCoords, trainData.cleanRois(roi).vox.linearCoords)) ];

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeVarCovarLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeVarCovarLikelihood(Rho,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean)

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeVarCovarMatrix = nearestSPD( (1-Rho)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeVarCovarMatrix));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getJeheeFullLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getJeheeFullLikelihood(params,covarianceMatrix,numVoxels,pRFresidualtSeries,pdfMean,rfOverlapMatrix)
    
    Rho = params(1); Sigma = params(2)^2;

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    jeheeFullMatrix = nearestSPD( (1-Rho-Sigma)*covarianceMatrix.*eye(numVoxels) + Rho*covarianceMatrix + Sigma*rfOverlapMatrix.*covarianceMatrix);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, jeheeFullMatrix));



%%%%%%%%%%%%%%%%%%%
%% getCorFromFit %%
%%%%%%%%%%%%%%%%%%%

function y = getCorFromFit(x,expFit)

% output of each matrix entry is the exponential fit to the data
y = expFit.a*exp(expFit.b*x) + expFit.c*exp(expFit.d*x);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getOverlapDistanceLikelihood %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = getOverlapDistanceLikelihood(params,pRFresidualtSeries,pdfMean,distanceVarCovarMatrix,rfOverlapVarCovarMatrix)
    
    Nu = params(1);

    % calc varCovar matrix as a weighted sum (by covariance Rho) of variance and covariance matrices
    overlapDistanceVarCovarMatrix = nearestSPD( Nu*rfOverlapVarCovarMatrix + (1-Nu)*distanceVarCovarMatrix);

    % calc likelihood with Jehee-inspired variance covariance matrix
    logLikelihood = 10^8-sum( logmvnpdf( pRFresidualtSeries, pdfMean, overlapDistanceVarCovarMatrix));



%%%%%%%%%%%%%
%% showFit %%
%%%%%%%%%%%%%
function showFit(x,y,fitToData,type)
figure(1);

% plot RF overlap
if type == 1; subplot(1,2,1)
scatter(x,y,1,'filled','k'); hold on;
graphFit = plot(fitToData,'predobs'); for i = 1:3, graphFit(i).Color = [0, 0.4470, 0.7410]; graphFit(i).LineWidth = 2; end
for i = 2:3, graphFit(i).LineStyle = '--'; graphFit(i).LineWidth = .75; end;
title('RF Overlap and residual correlation');xlabel('RF overlap (percent)'),ylabel('Residual Correlation'); ylim([-.5 1]);

% plot distance
elseif type == 2; subplot(1,2,2)
scatter(x,y,1,'filled','k'); hold on;
graphFit = plot(fitToData,'predobs'); for i = 1:3, graphFit(i).Color = [0, 0.4470, 0.7410]; graphFit(i).LineWidth = 2; end
for i = 2:3, graphFit(i).LineStyle = '--'; graphFit(i).LineWidth = .75; end;
title('Distance and residual correlation');xlabel('Voxel Distance (mm)'),ylabel('Residual Correlation'); ylim([-.5 1]);

end



%%%%%%%%%%%%%%%%%
%% imshowScale %%
%%%%%%%%%%%%%%%%%
% just a better version of rescale that multiplies the matrix by 1/avg(diag) so it's visible
function outputMatrix = imshowScale(inputMatrix)
outputMatrix = inputMatrix.*1/mean(diag(inputMatrix));






%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
% OUTDATED HELPER FUNCTIONS %
%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
% I don't call any of these but used them at some point. Keeping for documentation reasons, I guess.


%%%%%%%%%%%%%%%%%%%%%
%% fitVariance %%
function logLikelihood = fitVariance(params,testpRFresidualtSeries,pdfMean,numVoxels)

varianceMatrix = nearestSPD((params' * params) .* eye(numVoxels));
logLikelihood = 10^8 - sum(logmvnpdf( testpRFresidualtSeries, pdfMean, varianceMatrix));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makePositiveSemiDefinite %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Josh 8/4/22: I've replaced this with nearestSPD from John D'Errico(https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd).
% It seems to have less of an effect on the log likelihood values. I think he is making smaller adjustments than my brute forcing.
%
% Makes a matrix positive semi definite by adding the minimum eigenvalue. This occurs in matlab with very small values due to rounding errors.
% If the minmum eigenvalue is negative enough to have potentially arisen from something other than rounding error, it will offer you a warning.
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

