%% synthCor.m %%
%
%       By: Josh Wilson
%       Created: April 2022
%
% Pipeline:
%   1. Generate random pRF parameters.
%   2. Synthesize true time series from encoding model (specified in input).
%   3. Add noise to synthesized time series
%   4. Fit a simple gaussian receptive field to the noisy true time series by minimizing squared residual error.
%       You can add other receptive field models to synthesize or decode.
%   5. Plot relationship between receptive field overlap (recovered) and noise correlation between voxels.
%
% Things you can specify:
%   Visual field:
%       fieldSize: Size of the receptive field (pixels). True RF parameters are generated based on this size.
%       volumes: How fast the stimulus travels (time for full sweep across RF). For this and hdr, I'm thinking in seconds.
%       sweeps: Numer of stimulus sweeps across RF. Should be minimum ~10 (bars sweep randomly horizontal/vertical in forward/reverse).
%   Noise:
%       param.globalNoise: 1 if you want to add global Orenstein-Uhlenbeck noise (brownian motion with return to 0).
%           param.globalNoiseMag: std of time-correlated global noise.
%           param.theta: 0-1 OU process parameter; strength of pull back to 0 at each time point of brownian motion (set 1 = plain gaussian noise; set 0 = brownian motion).
%       param.gaussianNoise: 1 to add gaussian noise to individual voxels.
%           param.varStd/param.varAvg: mean/std of additive gaussian noise.
%           param.rho: Covariance of noise by receptive field overlap (0-1). 0 = no covariance (only variance). 1 = full covariance (currently broken).
%   Alternate model recovery:
%       negDoGFactor: For DoG encoding model, divide the surround RF by this value before subtracting.
%       DoGsize: Size of surround DoG receptive field relative to center.
%   Random:
%       rectify: If 1, will perform ReLU on DoG RF. If 0, won't.
%
%
% Usage:
%   synthCor(numVoxels,encodingModel)
%
%   Example:
%       [vox param] = synthCor(1000, 'gaussianDiff')


function [vox param rfOverlapRec noiseCor correlatedNoise] = synthCor(numVoxels,synthRFType,recRFType);


%% Set parameters %%
param.fieldSize = 80; 
param.volumes = 20; 
param.sweeps = 10;

param.globalNoise = 1; %1 to add correlated noise to all time series
param.globalNoiseMag = .05;
param.theta = .2; %strength at which global noise returns to 0

param.gaussianNoise = 1;
param.varAvg = .05; %mean of individual voxel variance
param.varStd = param.varAvg/10; %std of individual voxel variance
param.rho = .2; %gaussian noise individual covariance (1 = dependent)

param.negDoGFactor = 2; %magnitude of the negative difference of gaussian factor
param.DoGsize = 2; %size of the inhibitory gaussian receptive field relative to excitatory (std)
param.rectify = 1; %relu rectification for negative receptive field pixels

param.horizontal = round(rand(1,param.sweeps)); 
param.reverse = round(rand(1,param.sweeps)); 
vox = [];
param.p.canonical.type = 'gamma';
param.p.canonical.lengthInSeconds = 20;
param.p.canonical.timelag = 1;
param.p.canonical.tau = .6;
param.p.canonical.exponent = 6;
param.p.canonical.offset = 0;
param.numVols = param.sweeps*param.volumes+param.p.canonical.lengthInSeconds;
param.rfSynthType = synthRFType;
param.rfRecType = recRFType;




%% Synthesize the receptive fields %%
for voxel = 1:numVoxels
    vox = synthRF(vox,voxel,param);
end


%% Create the correlated noise structure to add to voxels %%
[correlatedNoise, param] = correlateNoise(param,vox,numVoxels);


%% Make the stimulus Movie %%
stimMovie = makeStimMovie(param);


%% Synth Time series w/ noise and fit receptive fields to voxels %% 
for voxel = 1:numVoxels
    vox = recoverRF(vox,voxel,param,stimMovie,correlatedNoise);
    [a, MSGID] = lastwarn(); warning('off', MSGID); %turn off LM warning
    if mod(voxel/numVoxels*100,20) == 0
        fprintf('Synthesizing voxels: %i percent done\n', voxel/numVoxels*100)
    end
end


%% Receptive field Overlap %%
rfOverlapRec = zeros(numVoxels);
for row = 1:numVoxels
    for column = 1:numVoxels
        mu1 = 0; s1 = vox{row}.Rparams(3); s2 = vox{column}.Rparams(3);
        mu2 = sqrt((vox{row}.Rparams(1)-vox{column}.Rparams(1))^2 + (vox{row}.Rparams(2)-vox{column}.Rparams(2))^2);
        if mu1 == mu2 & s1 == s2; 
            rfOverlapRec(row,column) = 1;
        else;
            c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
            rfOverlapRec(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end
    end
end
        
%% Noise Correlation %%
for row = 1:numVoxels
    for column = 1:numVoxels
        noiseCor(row,column) = corr2(vox{row}.noiseSeries,vox{column}.noiseSeries);
    end
end


%% plot %%
figure;hold on;
for i = 1:length(noiseCor); scatter(rfOverlapRec(i,:),noiseCor(i,:)); end;
rfOverlapArr = reshape(rfOverlapRec,[1 length(rfOverlapRec)^2]); NoiseCorArr = reshape(noiseCor,[1 length(noiseCor)^2]);
[rfOverlapArr,sortOrder] = sort(rfOverlapArr); NoiseCorArr = NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .01;
for bin = 0:step:40
    if sum( (bin < rfOverlapArr) & (rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(NoiseCorArr((bin < rfOverlapArr) & (rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);
title('Receptive Field Overlap and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j'); ylabel('Noise correlation between voxels i,j');

leg = {sprintf('Synth RF type: %s',param.rfSynthType), 
    sprintf('Recovery RF Type: %s',param.rfRecType), 
    sprintf('Global Noise: %1.0i',param.globalNoise), 
    sprintf('      Magnitude: %0.2f',param.globalNoiseMag),
    sprintf('Gaussian Noise: %1.0i', param.gaussianNoise), 
    sprintf('      Rho: %0.1f',param.rho),
    sprintf('      Magnitude: %0.1f',param.varAvg)};

text(.05,.8,leg)

%% compared recovered to group truth parameters %%
figure
subplot(2,2,1); hold on; for voxel = 1:length(vox); scatter(vox{voxel}.OGparams(1),vox{voxel}.Rparams(1)); end;
xlabel('Ground Truth X');ylabel('Recovered X'); plot([0 param.fieldSize], [0 param.fieldSize],'k');

subplot(2,2,2); hold on; for voxel = 1:length(vox); scatter(vox{voxel}.OGparams(2),vox{voxel}.Rparams(2)); end;
xlabel('Ground Truth Y');ylabel('Recovered Y'); plot([0 param.fieldSize], [0 param.fieldSize],'k');

subplot(2,2,3); hold on; for voxel = 1:length(vox); scatter(vox{voxel}.OGparams(3),vox{voxel}.Rparams(3)); end;
xlabel('Ground Truth Std');ylabel('Recovered Std'); plot([param.fieldSize*.05 param.fieldSize*.25], [param.fieldSize*.05 param.fieldSize*.25],'k');

rfOverlapRecvTrue = recTrueOverlap(vox,numVoxels);
subplot(2,2,4); histogram(rfOverlapRecvTrue); title('True and Recovered RF overlap');
%% end of program %%









%%%%%%%%%%%%%%%%%%%%
%% correlateNoise %%
%%%%%%%%%%%%%%%%%%%%
function [correlatedNoise param] = correlateNoise(param,vox,numVoxels)


%% global Ornstein Uhlenbeck noise %%
if param.globalNoise
    correlatedNoise.global(1) = 0;
    for volume = 2:param.numVols
      correlatedNoise.global(volume) = correlatedNoise.global(volume-1) - param.theta*correlatedNoise.global(volume-1) + normrnd(0,param.globalNoiseMag);
    end
elseif param.globalNoise == 0; 
    correlatedNoise.global(1:param.numVols) = 0; 
end;


%% correlated indivudual voxel noise %%
rfOverlapTrue = getTrueOverlap(vox,numVoxels);
rfOverlapTrue = rfOverlapTrue.*rfOverlapTrue; % fuck it, let's square this for now so it's not linear

param.tau = abs(normrnd(param.varAvg,param.varStd,numVoxels,1));
stdMatrix = param.tau * param.tau';

covarMatrix = param.rho * (stdMatrix .* rfOverlapTrue);
varMatrix = (1-param.rho) * (stdMatrix .* rfOverlapTrue .* eye(numVoxels));
varCovarMatrix = covarMatrix + varMatrix;


%% generate the noise %%
correlatedNoise.individual = mvnrnd(zeros(1,numVoxels), varCovarMatrix, param.numVols)


%%%%%%%%%%%%%%%%
%% synthRF %%
%%%%%%%%%%%%%%%%
function vox = synthRF(vox,voxel,param)


%% Make the receptive field %%
% Gaussian receptive field %
if strcmp(param.rfSynthType,'gaussian')
    x = param.fieldSize*(.1) + (param.fieldSize*(.8)).*rand;
    y = param.fieldSize*(.1) + (param.fieldSize*(.8)).*rand;
    sx = param.fieldSize*(.05) + (param.fieldSize*(.2)).*rand; sy=sx;
    gaussian1 = normpdf([1:param.fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:param.fieldSize]),y,sy);
    rf = gaussian2*gaussian1;
% Difference of Gaussians Receptive Field %
elseif strcmp(param.rfSynthType,'gaussianDiff')
    x = param.fieldSize*(.1) + (param.fieldSize*(.8)).*rand;
    y = param.fieldSize*(.1) + (param.fieldSize*(.8)).*rand;
    sx = param.fieldSize*(.05) + (param.fieldSize*(.2)).*rand; sy=sx; sx2 = param.DoGsize*sx; sy2 = param.DoGsize*sy;
    gaussian1 = normpdf([1:param.fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:param.fieldSize]),y,sy);
    gaussianDiff1 = normpdf([1:param.fieldSize],x,sx2);
    gaussianDiff2 = normpdf(transpose([1:param.fieldSize]),y,sy2);
    rf = gaussian2*gaussian1 - (gaussianDiff2*gaussianDiff1)/param.negDoGFactor;
    if param.rectify; % 
        rf(rf<0) = 0;
    end
end


%% grab things to return %%
vox{voxel}.OGparams = [x y sx];
vox{voxel}.OGrf = rf;



%%%%%%%%%%%%%%%
%% recoverRF %%
%%%%%%%%%%%%%%%
function vox = recoverRF(vox,voxel,param,stimMovie,correlatedNoise)

%% get signal time series %%
hrf = getCanonicalHRF(param.p.canonical,1);
tseries = getModelTSeries(stimMovie,vox{voxel}.OGrf,param,hrf);


%% add the gaussian noise, correlated by rho %%
if param.gaussianNoise;
    tempseries = tseries+transpose(correlatedNoise.individual(:,voxel));
else tempseries = tseries; end


%% add the global noise %%
noisytSeries = tempseries + correlatedNoise.global*std(tseries);


%% fit RF params to noisytSeries %%
startparams(1) = param.fieldSize/2; startparams(2) = param.fieldSize/2; startparams(3) = param.fieldSize/10; %start in center with param.fieldSize/10 std

minsearch = [param.fieldSize/20 param.fieldSize/20 param.fieldSize/50]; maxsearch = [param.fieldSize*19/20 param.fieldSize*19/20 param.fieldSize]; opts = optimset('display','off'); %constrain search to visual field and >.5 std

[params] = lsqnonlin(@getModelResidual,startparams,minsearch,maxsearch,opts,noisytSeries,hrf,stimMovie,param);


%% draw the recovered receptive field %%
if strcmp(param.rfRecType,'gaussian');
    recgaussian1 = normpdf([1:param.fieldSize],params(1),params(3));
    recgaussian2 = normpdf(transpose([1:param.fieldSize]),params(2),params(3));
    recrf = recgaussian2*recgaussian1; %recovers a simple gaussian RF
    rectSeries = getModelTSeries(stimMovie,recrf,param,hrf);
    param.rfRecoveryType = 'gaussian';
elseif strcmp(param.rfRecType,'gaussianDiff');
    gaussian1 = normpdf([1:param.fieldSize],params(1),params(3));
    gaussian2 = normpdf(transpose([1:param.fieldSize]),params(2),params(3));
    sx2 = param.DoGsize*params(3); sy2 = param.DoGsize*params(3);
    gaussianDiff1 = normpdf([1:param.fieldSize],params(1),sx2);
    gaussianDiff2 = normpdf(transpose([1:param.fieldSize]),params(2),sy2);
    recrf = gaussian2*gaussian1 - (gaussianDiff2*gaussianDiff1)/param.negDoGFactor;
    rectSeries = getModelTSeries(stimMovie,recrf,param,hrf);
    if param.rectify; % 
        recrf(recrf<0) = 0;
    end
end

%% grab things to return %%
vox{voxel}.hrf = hrf;
vox{voxel}.Rparams = params;
vox{voxel}.Rrf = recrf;
vox{voxel}.trueSeries = tseries;
vox{voxel}.noisyTrueSeries = noisytSeries;
vox{voxel}.recoveredSeries = rectSeries;
vox{voxel}.noiseSeries = noisytSeries-rectSeries;



%%%%%%%%%%%%%%%%%%%%%%
%% getModelResidual %%
%%%%%%%%%%%%%%%%%%%%%%
function residual = getModelResidual(params,noisytSeries,hrf,stimMovie,param)

%% make the receptive field from params you are fitting%%
if strcmp(param.rfRecType,'gaussian');
    recgaussian1 = normpdf([1:param.fieldSize],params(1),params(3));
    recgaussian2 = normpdf(transpose([1:param.fieldSize]),params(2),params(3));
    rf = recgaussian2*recgaussian1; %recovers a simple gaussian RF
    param.rfRecoveryType = 'gaussian';
elseif strcmp(param.rfRecType,'gaussianDiff');
    gaussian1 = normpdf([1:param.fieldSize],params(1),params(3));
    gaussian2 = normpdf(transpose([1:param.fieldSize]),params(2),params(3));
    sx2 = param.DoGsize*params(3); sy2 = param.DoGsize*params(3);
    gaussianDiff1 = normpdf([1:param.fieldSize],params(1),sx2);
    gaussianDiff2 = normpdf(transpose([1:param.fieldSize]),params(2),sy2);
    rf = gaussian2*gaussian1 - (gaussianDiff2*gaussianDiff1)/param.negDoGFactor;
    if param.rectify; % 
        rf(rf<0) = 0;
    end
end


%% signal time series %%
% neural response
for t = 1:param.volumes*param.sweeps
    ntseries(t) = sum(sum(stimMovie(:,:,t).*rf));
    time(t) = t;
end
tseries = conv(ntseries,hrf.hrf);
residual = sumsqr(noisytSeries - tseries);



%%%%%%%%%%%%%%%%%%%%%
%% getModelTSeries %%
%%%%%%%%%%%%%%%%%%%%%
function tseries = getModelTSeries(stimMovie,rf,param,hrf)
% neural response
for t = 1:param.volumes*param.sweeps
    ntseries(t) = sum(sum(stimMovie(:,:,t).*rf));
    time(t) = t;
end

% convolve neural and hrf
tseries = conv(ntseries,hrf.hrf);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);



%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit



%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);



%%%%%%%%%%%%%%%%%%%
%% makeStimMovie %%
%%%%%%%%%%%%%%%%%%%
function stimMovie = makeStimMovie(param)

stimMovie = zeros(param.fieldSize,param.fieldSize,param.volumes*param.sweeps);

for sweep = 1:param.sweeps
    for volume = 1:(param.volumes)
        
        if param.horizontal(sweep)
            if param.reverse(sweep)
                stimMovie((param.fieldSize/param.volumes*(param.volumes-volume+1)-param.fieldSize/param.volumes+1):(param.fieldSize/param.volumes*(param.volumes-volume+1)), : ,(sweep-1)*param.volumes+volume) = 1;
            else
                stimMovie((param.fieldSize/param.volumes*(volume)-param.fieldSize/param.volumes+1):param.fieldSize/param.volumes*(volume), : ,(sweep-1)*param.volumes+volume) = 1;
            end
        else
            if param.reverse(sweep)
                stimMovie(:, (param.fieldSize/param.volumes*(param.volumes-volume+1)-param.fieldSize/param.volumes+1):(param.fieldSize/param.volumes*(param.volumes-volume+1)), (sweep-1)*param.volumes+volume) = 1;
            else
                stimMovie(:, (param.fieldSize/param.volumes*(volume)-param.fieldSize/param.volumes+1):param.fieldSize/param.volumes*(volume), (sweep-1)*param.volumes+volume) = 1;
            end
        end

    end
end



%%%%%%%%%%%%%%%%%%%%%
%% getTrueOverlap %%
%%%%%%%%%%%%%%%%%%%%
function rfOverlapTrue = getTrueOverlap(vox,numVoxels)

rfOverlap = zeros(numVoxels);
for row = 1:numVoxels
    for column = 1:numVoxels
        mu1 = 0; s1 = vox{row}.OGparams(3); s2 = vox{column}.OGparams(3);
        mu2 = sqrt((vox{row}.OGparams(1)-vox{column}.OGparams(1))^2 + (vox{row}.OGparams(2)-vox{column}.OGparams(2))^2);
        if mu1 == mu2 & s1 == s2; 
            rfOverlap(row,column) = 1;
        else;
            c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
            rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end
    end
end

rfOverlapTrue = 0.5 * (rfOverlap + rfOverlap');

%n = size(rfOverlap,1);
%cvx_begin
%variable rfOverlapTrue(n,n)
%minimize(norm(rfOverlapTrue-rfOverlap,'fro'))
%rfOverlapTrue -m *eye(n) == semidefinite(n)
%cvx_end

%%%%%%%%%%%%%%%%%%%%
%% recTrueOverlap %%
%%%%%%%%%%%%%%%%%%%%
%Get overlap between recovered and true parameters
function rfOverlapRecvTrue = recTrueOverlap(vox,numVoxels)

rfOverlapRecvTrue = zeros(1,numVoxels);
for voxel = 1:numVoxels;
    mu1 = 0; s1 = vox{voxel}.OGparams(3); s2 = vox{voxel}.Rparams(3);
    mu2 = sqrt((vox{voxel}.OGparams(1)-vox{voxel}.Rparams(1))^2 + (vox{voxel}.OGparams(2)-vox{voxel}.Rparams(2))^2);
    if mu1 == mu2 & s1 == s2; 
        rfOverlapRecvTrue(voxel) = 1;
    else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        rfOverlapRecvTrue(voxel) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end
end




function showSynthVoxels
figure;hold on;
scatter(1:length(vox{v}.noisyTrueSeries),vox{v}.noisyTrueSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.recoveredSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.trueSeries,'black')


%% convex semi-definite program solution to make varCovar positive semidefinite %%
%n = size(vcMatrix,1);
%cvx_begin
%variable varCovarMatrix(n,n)
%minimize(norm(varCovarMatrix-vcMatrix,'fro'))
%(varCovarMatrix -0*eye(n)) == semidefinite(n)
%cvx_end