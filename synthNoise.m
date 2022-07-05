%% synthNoise.m %%
%
%       By: Josh Wilson
%       Created: June 2022
%
% Pipeline:
%   First part does the same thing as synthCor.m - creates receptive fields and time series, recovers, etc.
%
%   After that, creates several more time series for each voxel from the same ground truth parameters
%   (can specify number of extra scans with param.numScans). Then, averages those scans together to use 
%   as a mean model of voxel activity. Subtracts the mean model from each voxels and shows you the noise
%   characteristics (which you can set the same way you would in synthCor). The interesting part is adding
%   multiplicative noise (param.multi). 
%
%
%   Example:
%       [vox param rfOverlapRec noiseCor correlatedNoise stimMovie] = synthNoise(100,'gaussian','gaussian');


function [vox param rfOverlapRec noiseCor correlatedNoise, stimMovie] = synthCor(synthRFType,varargin);

getArgs(varargin);


%% Set parameters %%
param.fieldSize = 80; 
param.volumes = 20; 
param.sweeps = 10;
%param.numScans = 3;

param.gaussianNoise = 1;
param.varAvg = .1; %mean of individual voxel variance
param.varStd = param.varAvg/10; %std of individual voxel variance
param.rho = 0; %global covariance
param.sigma = 0; %channel dependent noise contribution
param.multi = 0; %set 1 for multiplicative noise
param.multiDegree = .2; % exponent for multiplicative noise

param.globalNoise = 0; %1 to add correlated noise to all time series
param.globalNoiseMag = .03;
param.theta = 0; %strength at which global noise returns to 0

param.gaussianNoiseTime = 0;
param.varAvgT = .05; %mean of individual voxel variance
param.varStdT = param.varAvgT/10; %std of individual voxel variance
param.rhoT = .4; %gaussian noise individual covariance (1 = dependent)
param.thetaT = .2;

param.rectify = 0; %relu for negative RF values
param.DoGsize = 2;
param.negDoGFactor = 2;

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
%param.rfRecType = recRFType;

getArgs(varargin);
param.numScans = numScans;
param.multi = multi;
param.multiDegree = degree;



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



%% make the extra time series %%
for scan = 1:param.numScans
    [correlatedNoise, param] = correlateNoise(param,vox,numVoxels);
    for voxel = 1:numVoxels
        tseries = vox{voxel}.trueSeries;
        %% add the gaussian noise, correlated by rho %%
        if param.gaussianNoise;
        if param.multi %multiplicative noise
            tempnoise = transpose(correlatedNoise.individual(:,voxel));
            noise = tempnoise.*((tseries/max(tseries)+1).^param.multiDegree);
            noise = noise/(std(noise)/std(tempnoise)); %standardize
            tempseries = tseries+noise;
        else
            tempseries = tseries+transpose(correlatedNoise.individual(:,voxel));
        end
        else tempseries = tseries; end
        if param.gaussianNoiseTime
            tempseries2 = tempseries + transpose(correlatedNoise.individualTime(:,voxel));
        else tempseries2 = tempseries; end

        %% add the global noise %%
        vox{voxel}.dupeScans{scan}.noisyTrueSeries  = tempseries2 + correlatedNoise.global*std(tseries);

    end
end

%% plot the noise %%
groupsize = param.numScans;
noisyTrueSeries = []; avg = []; modelSeries = [];
for group = 1:numVoxels;
    avgseries = zeros(1,length(vox{1}.trueSeries));
    for scan = 1:groupsize, avgseries = avgseries+vox{group}.dupeScans{scan}.noisyTrueSeries;end
avgseries = avgseries/groupsize; vox{group}.avgSeries = avgseries;

avg = [avg zscore(avgseries)];
noisyTrueSeries = [noisyTrueSeries zscore(vox{group}.noisyTrueSeries)];
modelSeries = [modelSeries zscore(vox{group}.trueSeries)];
residuals = noisyTrueSeries-avg;
end

figure,subplot(2,2,1),scatter(modelSeries,noisyTrueSeries-avg,1,'filled','k'); xlabel('model time series % signal'),ylabel('residual (true - avg)'); title('Residuals by pRF predictions')
subplot(2,2,2),scatter(noisyTrueSeries,noisyTrueSeries-avg,1,'filled','k'); xlabel('true time series % signal'); ylabel('residual (true - avg)'); title('Residuals by time series activity')
subplot(2,2,3),scatter(avg,noisyTrueSeries-avg,1,'filled','k'); xlabel('avg time series % signal'); ylabel('residual (true - avg)'); title('Residuals by mean model prediction')
subplot(2,2,4),scatter(noisyTrueSeries,avg,1,'filled','k'); xlabel('Single Scan % signal'); ylabel('Avg % signal'); title('Single scan and average activity'); hold on;


startparams(1) = 0; startparams(2) = 1; startparams(3) = 1; startparams(4) = 0;

minsearch = [-1 -inf -inf -inf]; maxsearch = [1 inf inf inf]; opts = optimset('display','off'); opts = optimset('maxIter',2500);

[params, logLikelihood] = lsqnonlin(@fitNoiseParameters,startparams,minsearch,maxsearch,opts,residuals,numScans,avg,noisyTrueSeries)










keyboard





%%%%%%%%%%%%%%%%%%%
%% END OF SCRIPT %%
%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%% fitNoiseParameters %%
%%%%%%%%%%%%%%%%%%%%%%%%
function logLikelihood = fitNoiseParameters(params,residuals,numScans,avg,noisyTrueSeries)
mu = params(1); compstd = params(2); scale = params(3); beta = params(4);


std = ((normcdf(noisyTrueSeries)*scale+beta)+compstd).^2+(compstd/sqrt(numScans))^2;
logLikelihood = sum(log(normpdf(residuals,mu,std)));





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


%% gaussian time-independent voxel noise %%
rfOverlapTrue = getTrueOverlap(vox,numVoxels);
rfOverlapTrue = rfOverlapTrue.*rfOverlapTrue; % fuck it, let's square this - looks more like real data

param.tau = abs(normrnd(param.varAvg,param.varStd,numVoxels,1));
stdMatrix = param.tau * param.tau';

covarMatrix = param.rho * (stdMatrix);
varMatrix = (1-param.rho) * (stdMatrix .* eye(numVoxels));

channelMatrix = param.sigma^2 * (stdMatrix .* rfOverlapTrue);

varCovarMatrix = covarMatrix + varMatrix + channelMatrix;

correlatedNoise.individual = mvnrnd(zeros(1,numVoxels), varCovarMatrix, param.numVols);


%% individual Ornstein Uhlenbeck noise %%
if param.gaussianNoiseTime
    %% regenerate correlated noise %%
    covarMatrix = param.rhoT * (stdMatrix .* rfOverlapTrue);
    varMatrix = (1-param.rhoT) * (stdMatrix .* rfOverlapTrue .* eye(numVoxels));
    varCovarMatrix = covarMatrix + varMatrix;
    correlatedNoiseTemp = mvnrnd(zeros(1,numVoxels), varCovarMatrix, param.numVols);
    
    for voxel = 1:numVoxels
        correlatedNoise.individualTime(1,voxel) = correlatedNoiseTemp(1,voxel);
        for time = 2:param.numVols
            correlatedNoise.individualTime(time,voxel) = correlatedNoise.individualTime(time-1,voxel) - param.thetaT*correlatedNoise.individualTime(time-1,voxel) + correlatedNoiseTemp(time,voxel);
        end
    end

end
    


%%%%%%%%%%%%%%%%
%% synthRF %%
%%%%%%%%%%%%%%%%
function vox = synthRF(vox,voxel,param)


%% Make the receptive field %%
% Gaussian receptive field %
if strcmp(param.rfSynthType,'gaussian')
    x = param.fieldSize*(.3) + (param.fieldSize*(.4)).*rand;
    y = param.fieldSize*(.3) + (param.fieldSize*(.4)).*rand;
    sx = (param.fieldSize*(.025) + (param.fieldSize*(.025)).*rand) * 2 * (1 + sqrt((x-param.fieldSize/2)^2+(y-param.fieldSize/2)^2) / sqrt((param.fieldSize/2)^2 + (param.fieldSize/2)^2) ); sy=sx;
    %x=40;y=40;sx=5;sy=5;
    gaussian1 = normpdf([1:param.fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:param.fieldSize]),y,sy);
    rf = gaussian2*gaussian1;
% Difference of Gaussians Receptive Field %
elseif strcmp(param.rfSynthType,'gaussianDiff')
    x = param.fieldSize*(.2) + (param.fieldSize*(.6)).*rand;
    y = param.fieldSize*(.2) + (param.fieldSize*(.6)).*rand;
    sx = (param.fieldSize*(.04) + (param.fieldSize*(.04)).*rand) * 2 * (1 + sqrt((x-param.fieldSize/2)^2+(y-param.fieldSize/2)^2) / sqrt((param.fieldSize/2)^2 + (param.fieldSize/2)^2) ); sy=sx; sx2 = param.DoGsize*sx; sy2 = param.DoGsize*sy;
    %x=40;y=40;sx=8;sy=8;sx2 = sx*param.DoGsize; sy2 = sy*param.DoGsize;
    gaussian1 = normpdf([1:param.fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:param.fieldSize]),y,sy);
    gaussianDiff1 = normpdf([1:param.fieldSize],x,sx2);
    gaussianDiff2 = normpdf(transpose([1:param.fieldSize]),y,sy2);
    rf = gaussian2*gaussian1 - (gaussianDiff2*gaussianDiff1)/param.negDoGFactor;
    vox{voxel}.g1 = gaussian2*gaussian1; vox{voxel}.g2 = (gaussianDiff2*gaussianDiff1)/param.negDoGFactor;
end

    if param.rectify; % 
        rf(rf<0) = 0;
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
[tseries ntseries] = getModelTSeries(stimMovie,vox{voxel}.OGrf,param,hrf);


%% add the gaussian noise, correlated by rho %%
if param.gaussianNoise;
 
    if param.multi %multiplicative noise
    tempnoise = transpose(correlatedNoise.individual(:,voxel));
    noise = tempnoise.*((tseries/max(tseries)+1).^param.multiDegree);
    noise = noise/(std(noise)/std(tempnoise)); %
    tempseries = tseries+noise;
    else
    tempseries = tseries+transpose(correlatedNoise.individual(:,voxel));
    end

else tempseries = tseries; end

if param.gaussianNoiseTime
    tempseries2 = tempseries + transpose(correlatedNoise.individualTime(:,voxel));
else tempseries2 = tempseries; end


%% add the global noise %%
noisytSeries = tempseries2 + correlatedNoise.global*std(tseries);

%if you want to induce hrf recovery noise
%param.p.canonical.tau = .75;
%hrf = getCanonicalHRF(param.p.canonical,1);

%% fit RF params to noisytSeries %%
startparams(1) = param.fieldSize/2; startparams(2) = param.fieldSize/2; startparams(3) = param.fieldSize/10; %start in center with param.fieldSize/10 std

minsearch = [param.fieldSize/20 param.fieldSize/20 param.fieldSize/50]; maxsearch = [param.fieldSize*19/20 param.fieldSize*19/20 param.fieldSize]; opts = optimset('display','off'); %constrain search to visual field and >.5 std

recover = 0; if recover
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
end
%% grab things to return %%
vox{voxel}.hrf = hrf;
%vox{voxel}.Rparams = params;
%vox{voxel}.Rrf = recrf;
vox{voxel}.trueSeries = tseries;
vox{voxel}.noisyTrueSeries = noisytSeries;
%vox{voxel}.recoveredSeries = rectSeries;
%vox{voxel}.noiseSeries = noisytSeries-rectSeries;
%vox{voxel}.ntseries = ntseries;

%param.p.canonical.timelag = 3;
%hrf = getCanonicalHRF(param.p.canonical,1);
%hrftseries = getModelTSeries(stimMovie,vox{voxel}.OGrf,param,hrf);
%vox{voxel}.trueSerieshrf = hrftseries;


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
function [tseries ntseries] = getModelTSeries(stimMovie,rf,param,hrf)
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
figure;
subplot(1,2,1);hold on;
scatter(1:length(vox{v}.noisyTrueSeries),vox{v}.noisyTrueSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.recoveredSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.trueSeries,'black');
subplot(2,2,2);imshow(rescale(vox{v}.OGrf)); title('True RF')
subplot(2,2,4);imshow(rescale(vox{v}.Rrf)); title('Recovered RF')


%% convex semi-definite program solution to make varCovar positive semidefinite %%
%n = size(vcMatrix,1);
%cvx_begin
%variable varCovarMatrix(n,n)
%minimize(norm(varCovarMatrix-vcMatrix,'fro'))
%(varCovarMatrix -0*eye(n)) == semidefinite(n)
%cvx_end