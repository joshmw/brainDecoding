%%%%%%%%%%%%%%%%%%%%%%%
%% saveSumarryFigs.m %%
%%%%%%%%%%%%%%%%%%%%%%%
%
% Uses the core noise correlation functions (getpRFTSeries, avgModelNoise, shufflecor) to compute summary statistics
% for all of the subjects for which we have data. The main comparisons are:
%
%     Noise correlations using the prf model and the mean moddle
%
%     Noise correlations between different scans using the pRF model
%
%     Same comparisons but will different levels of filtering



%Get the list of subjects - set this manually, stored in onedrive with getpRFTSeries format and naming convention
subjectList = {'s0350', 's0399', 's0401', 's0403', 's0404', 's0405', 's0406', 's0416', 's0419', 's0420', 's0421' 's0423'}
% Excluded subjects: 's0414',


%% Get the exponents of the base rfOverlap/NoiseCor fit with the prf model %%
%initialize empty arrays
prfExponents = [];
prfExponentConfLow = [];
prfExponentConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr =  strcat('data=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    [v1v3Exponent expConfInt] = getpRFTSeries('loadData=1',dataStr{1},'graphStuff=0'); close all;
    prfExponents = [prfExponents v1v3Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];
end


%% Get the exponents using another scan as the model of scan 1 %
% initialize empty arrays
MMExponents = [];
MMExponentConfLow = [];
MMExponentConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr1 = strcat('data1=', subjectList(subject), 'pRFGaussianHdrNMp1.mat');
    dataStr2 = strcat('data2=', subjectList(subject), 'pRFGaussianHdrNMp2.mat');
    baseStr = strcat('base=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    [v1v3Exponent expConfInt] = avgModelNoise(dataStr1{1}, dataStr2{1}, baseStr{1},'graphStuff=0');
    MMExponents = [MMExponents v1v3Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all
end


%% plot the data %%
figure, hold on

errorbar(MMExponents,prfExponents,prfExponentConfLow-prfExponents,prfExponentConfHigh-prfExponents,'LineStyle','none','Color','k','CapSize',2) 
errorbar(MMExponents,prfExponents,MMExponentConfLow-MMExponents,MMExponentConfHigh-MMExponents,'horizontal','LineStyle','none','Color','k','CapSize',2) 

maxVal = max(prfExponents); xlim([0 maxVal*1.2]); ylim([0 maxVal*1.2]);
plot([0 maxVal*1.05],[0 maxVal*1.05],'--','Color','k')

xlabel('Mean Model Correlation Fit Exponent');
ylabel('pRF Model Correlation Fit Exponent')
title('Overlap/Correlation Relationship Exponents')
drawPublishAxis('labelFontSize=14','yAxisOffset=-1/25');legend('off');



%% look at noise correlations between different scans with shufflecor %%
%initialize empty arrays
shuffleV1V3Exponents = [];
shuffleV1V3ExponentConfLow = [];
shuffleV1V3ExponentConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr1 = strcat('data1=', subjectList(subject), 'pRFGaussianHdrNMp1.mat');
    dataStr2 = strcat('data2=', subjectList(subject), 'pRFGaussianHdrNMp2.mat');
    baseStr = strcat('base=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    [v1v3Exponent v1v3ExponentConfInt] = shufflecor(dataStr1{1},dataStr2{1},baseStr{1},'graphStuff=0');
    shuffleV1V3Exponents = [shuffleV1V3Exponents v1v3Exponent]; shuffleV1V3ExponentConfLow = [shuffleV1V3ExponentConfLow v1v3ExponentConfInt(1)]; shuffleV1V3ExponentConfHigh = [shuffleV1V3ExponentConfHigh v1v3ExponentConfInt(2)]; close all
end

%plot
figure; hold on
errorbar([1:length(shuffleV1V3Exponents)],shuffleV1V3Exponents,shuffleV1V3ExponentConfLow-shuffleV1V3Exponents,shuffleV1V3ExponentConfHigh-shuffleV1V3Exponents,'LineStyle','none','Color','k','CapSize',2)
