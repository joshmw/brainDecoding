%%%%%%%%%%%%%%%%
%% synthCor.m %%
%%%%%%%%%%%%%%%%

function synthCor(numVoxels,rfType)

fieldSize = 100;
volumes = 20;
sweeps = 10;
vox = [];

for voxel = 1:numVoxels
    vox = synthVoxel(vox,voxel,rfType,fieldSize,volumes,sweeps);
    if mod(voxel/numVoxels*100,20) == 0
        disp(sprintf('Synthesizing voxels: %i percent done', voxel/numVoxels*100))
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
keyboard

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
title('V1 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j');


keyboard






%%%%%%%%%%%%%%%%
%% synthVoxel %%
%%%%%%%%%%%%%%%%
function [vox, stimMovie] = synthVoxel(vox,voxel,rfType,fieldSize,volumes,sweeps)


%% set parameters %%
horizontal = round(rand(1,sweeps)); 
reverse = round(rand(1,sweeps)); 
fitParams.framePeriod = 1; %seconds
noiseMag = 5; %divide std of tseries by noisemag to generate gaussian noise: lower = more noise
negDoGFactor = 2; %magnitude of the negative difference of gaussian factor

p.canonical.type = 'gamma';
p.canonical.lengthInSeconds = 20;
p.canonical.timelag = 1;
p.canonical.tau = .6;
p.canonical.exponent = 6;
p.canonical.offset = 0;


%% Make the stimulus movie %%

stimMovie = zeros(fieldSize,fieldSize,volumes*sweeps);
for sweep = 1:sweeps
    for volume = 1:(volumes)
        
        if horizontal(sweep)
            if reverse(sweep)
                stimMovie((fieldSize/volumes*(volumes-volume+1)-fieldSize/volumes+1):(fieldSize/volumes*(volumes-volume+1)), : ,(sweep-1)*volumes+volume) = 1;
            else
                stimMovie((fieldSize/volumes*(volume)-fieldSize/volumes+1):fieldSize/volumes*(volume), : ,(sweep-1)*volumes+volume) = 1;
            end
        else
            if reverse(sweep)
                stimMovie(:, (fieldSize/volumes*(volumes-volume+1)-fieldSize/volumes+1):(fieldSize/volumes*(volumes-volume+1)), (sweep-1)*volumes+volume) = 1;
            else
                stimMovie(:, (fieldSize/volumes*(volume)-fieldSize/volumes+1):fieldSize/volumes*(volume), (sweep-1)*volumes+volume) = 1;
            end
        end

    end
end

%% Make the receptive field %%
%gaussian receptive field
if strcmp(rfType,'gaussian')
    x = round(fieldSize*(.1)) + round((fieldSize*(.8)).*rand);
    y = round(fieldSize*(.1)) + round((fieldSize*(.8)).*rand);
    sx = round(fieldSize*(.05)) + round((fieldSize*(.2)).*rand); sy=sx;
    gaussian1 = normpdf([1:fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:fieldSize]),y,sy);
    rf = gaussian2*gaussian1;
elseif strcmp(rfType,'gaussianDiff')
    x = round(fieldSize*(.1)) + round((fieldSize*(.8)).*rand);
    y = round(fieldSize*(.1)) + round((fieldSize*(.8)).*rand);
    sx = round(fieldSize*(.05)) + round((fieldSize*(.2)).*rand); sy=sx; sx2 = 2*sx; sy2 = 2*sy;
    gaussian1 = normpdf([1:fieldSize],x,sx);
    gaussian2 = normpdf(transpose([1:fieldSize]),y,sy);
    gaussianDiff1 = normpdf([1:fieldSize],x,sx2);
    gaussianDiff2 = normpdf(transpose([1:fieldSize]),y,sy2);
    rf = gaussian2*gaussian1 - (gaussianDiff2*gaussianDiff1)/negDoGFactor;
end


%% signal time series %%
hrf = getCanonicalHRF(p.canonical,fitParams.framePeriod);
tseries = getModelTSeries(stimMovie,rf,sweeps,volumes,hrf);
noisytSeries = tseries + normrnd(0,std(tseries)/noiseMag,[1 length(tseries)]);


%% recover time series from noisytSeries %%
startparams(1) = fieldSize/2; startparams(2) = fieldSize/2; startparams(3) = fieldSize/10;
[params] = fminsearch(@getModelResidual,startparams,[],noisytSeries,hrf,stimMovie,fieldSize,volumes,sweeps);

recgaussian1 = normpdf([1:fieldSize],params(1),params(3));
recgaussian2 = normpdf(transpose([1:fieldSize]),params(2),params(3));
recrf = recgaussian2*recgaussian1;
rectSeries = getModelTSeries(stimMovie,recrf,sweeps,volumes,hrf);

% grab things to return
vox{voxel}.hrf = hrf;
vox{voxel}.OGparams = [x y sx];
vox{voxel}.Rparams = params;
vox{voxel}.OGrf = rf;
vox{voxel}.Rrf = recrf;
vox{voxel}.trueSeries = tseries;
vox{voxel}.noisyTrueSeries = noisytSeries;
vox{voxel}.recoveredSeries = rectSeries;
vox{voxel}.noiseSeries = noisytSeries-rectSeries;
vox{voxel}.rfType = rfType;


%%%%%%%%%%%%%%%%%%%%%%
%% getModelResidual %%
%%%%%%%%%%%%%%%%%%%%%%
function residual = getModelResidual(params,noisytSeries,hrf,stimMovie,fieldSize,volumes,sweeps)
x = params(1); y = params(2); sx = params(3); sy = sx;
gaussian1 = normpdf([1:fieldSize],x,sx);
gaussian2 = normpdf(transpose([1:fieldSize]),y,sy);
rf = gaussian2*gaussian1;


%% signal time series %%
% neural response
for t = 1:volumes*sweeps
    ntseries(t) = sum(sum(stimMovie(:,:,t).*rf));
    time(t) = t;
end
tseries = conv(ntseries,hrf.hrf);
residual = sumsqr(noisytSeries - tseries);


%%%%%%%%%%%%%%%%%%%%%
%% getModelTSeries %%
%%%%%%%%%%%%%%%%%%%%%
function tseries = getModelTSeries(stimMovie,rf,sweeps,volumes,hrf)
% neural response
for t = 1:volumes*sweeps
    ntseries(t) = sum(sum(stimMovie(:,:,t).*rf));
    time(t) = t;
end

% convolve neural and hrf; add noise
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





function showSynthVoxels
figure;hold on;
scatter(1:length(vox{v}.noisyTrueSeries),vox{v}.noisyTrueSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.recoveredSeries);
plot(1:length(vox{v}.recoveredSeries),vox{v}.trueSeries,'black')


