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
prfExponentsv1 = []; prfExponentv1ConfLow = []; prfExponentv1ConfHigh = [];
prfExponentsv2 = []; prfExponentv2ConfLow = []; prfExponentv2ConfHigh = [];
prfExponentsv3 = []; prfExponentv3ConfLow = []; prfExponentv3ConfHigh = [];
prfExponentsv1v3 = []; prfExponentv1v3ConfLow = []; prfExponentv1v3ConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr =  strcat('data=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    %get the exponents
    [v1Exponent, v1ExponentConfInt, v2Exponent, v2ExponentConfInt, v3Exponent, v3ExponentConfInt, v1v2Exponent, v1v2ExponentConfInt, v1v3Exponent, v1v3ExponentConfInt] = getpRFTSeries('loadData=1',dataStr{1},'graphStuff=0'); close all;
    prfExponentsv1 = [prfExponentsv1 v1Exponent]; prfExponentv1ConfLow = [prfExponentv1ConfLow v1ExponentConfInt(1)]; prfExponentv1ConfHigh = [prfExponentv1ConfHigh v1ExponentConfInt(2)];
    prfExponentsv2 = [prfExponentsv2 v2Exponent]; prfExponentv2ConfLow = [prfExponentv2ConfLow v2ExponentConfInt(1)]; prfExponentv2ConfHigh = [prfExponentv2ConfHigh v2ExponentConfInt(2)];
    prfExponentsv3 = [prfExponentsv3 v3Exponent]; prfExponentv3ConfLow = [prfExponentv3ConfLow v3ExponentConfInt(1)]; prfExponentv3ConfHigh = [prfExponentv3ConfHigh v3ExponentConfInt(2)];
    prfExponentsv1v3 = [prfExponentsv1v3 v1v3Exponent]; prfExponentv1v3ConfLow = [prfExponentv1v3ConfLow v1v3ExponentConfInt(1)]; prfExponentv1v3ConfHigh = [prfExponentv1v3ConfHigh v1v3ExponentConfInt(2)];
end 


%% Get the exponents using another scan as the model of scan 1 %
% initialize empty arrays
MMExponentsv1 = []; MMExponentv1ConfLow = []; MMExponentv1ConfHigh = [];
MMExponentsv2 = []; MMExponentv2ConfLow = []; MMExponentv2ConfHigh = [];
MMExponentsv3 = []; MMExponentv3ConfLow = []; MMExponentv3ConfHigh = [];
MMExponentsv1v3 = []; MMExponentv1v3ConfLow = []; MMExponentv1v3ConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr1 = strcat('data1=', subjectList(subject), 'pRFGaussianHdrNMp1.mat');
    dataStr2 = strcat('data2=', subjectList(subject), 'pRFGaussianHdrNMp2.mat');
    baseStr = strcat('base=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    [v1Exponent, v1ExponentConfInt, v2Exponent, v2ExponentConfInt, v3Exponent, v3ExponentConfInt, v1v3Exponent, v1v3ExponentConfInt] = avgModelNoise(dataStr1{1}, dataStr2{1}, baseStr{1},'graphStuff=0'); close all
    MMExponentsv1 = [MMExponentsv1 v1Exponent]; MMExponentv1ConfLow = [MMExponentv1ConfLow v1ExponentConfInt(1)]; MMExponentv1ConfHigh = [MMExponentv1ConfHigh v1ExponentConfInt(2)];
    MMExponentsv2 = [MMExponentsv2 v2Exponent]; MMExponentv2ConfLow = [MMExponentv2ConfLow v2ExponentConfInt(1)]; MMExponentv2ConfHigh = [MMExponentv2ConfHigh v2ExponentConfInt(2)];
    MMExponentsv3 = [MMExponentsv3 v3Exponent]; MMExponentv3ConfLow = [MMExponentv3ConfLow v3ExponentConfInt(1)]; MMExponentv3ConfHigh = [MMExponentv3ConfHigh v3ExponentConfInt(2)];
    MMExponentsv1v3 = [MMExponentsv1v3 v1v3Exponent]; MMExponentv1v3ConfLow = [MMExponentv1v3ConfLow v1v3ExponentConfInt(1)]; MMExponentv1v3ConfHigh = [MMExponentv1v3ConfHigh v1v3ExponentConfInt(2)];
end


%% get exponents of correlations between scans %%
%initialize empty arrays
shuffleExponentsv1 = []; shuffleExponentv1ConfLow = []; shuffleExponentv1ConfHigh = [];
shuffleExponentsv2 = []; shuffleExponentv2ConfLow = []; shuffleExponentv2ConfHigh = [];
shuffleExponentsv3 = []; shuffleExponentv3ConfLow = []; shuffleExponentv3ConfHigh = [];
shuffleExponentsv1v3 = []; shuffleExponentv1v3ConfLow = []; shuffleExponentv1v3ConfHigh = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr1 = strcat('data1=', subjectList(subject), 'pRFGaussianHdrNMp1.mat');
    dataStr2 = strcat('data2=', subjectList(subject), 'pRFGaussianHdrNMp2.mat');
    baseStr = strcat('base=', subjectList(subject), 'pRFGaussianHdrNMFull.mat');
    [v1Exponent, v1ExponentConfInt, v2Exponent, v2ExponentConfInt, v3Exponent, v3ExponentConfInt, v1v3Exponent, v1v3ExponentConfInt] = shufflecor(dataStr1{1},dataStr2{1},baseStr{1},'graphStuff=0'); close all
    shuffleExponentsv1 = [shuffleExponentsv1 v1Exponent]; shuffleExponentv1ConfLow = [shuffleExponentv1ConfLow v1ExponentConfInt(1)]; shuffleExponentv1ConfHigh = [shuffleExponentv1ConfHigh v1ExponentConfInt(2)];
    shuffleExponentsv2 = [shuffleExponentsv2 v2Exponent]; shuffleExponentv2ConfLow = [shuffleExponentv2ConfLow v2ExponentConfInt(1)]; shuffleExponentv2ConfHigh = [shuffleExponentv2ConfHigh v2ExponentConfInt(2)];
    shuffleExponentsv3 = [shuffleExponentsv3 v3Exponent]; shuffleExponentv3ConfLow = [shuffleExponentv3ConfLow v3ExponentConfInt(1)]; shuffleExponentv3ConfHigh = [shuffleExponentv3ConfHigh v3ExponentConfInt(2)];
    shuffleExponentsv1v3 = [shuffleExponentsv1v3 v1v3Exponent]; shuffleExponentv1v3ConfLow = [shuffleExponentv1v3ConfLow v1v3ExponentConfInt(1)]; shuffleExponentv1v3ConfHigh = [shuffleExponentv1v3ConfHigh v1v3ExponentConfInt(2)];
end

%% get exponents from unfiltered data %%
%initialize empty arrays
prfExponentsv1UF = []; prfExponentv1ConfLowUF = []; prfExponentv1ConfHighUF = [];
prfExponentsv2UF = []; prfExponentv2ConfLowUF = []; prfExponentv2ConfHighUF = [];
prfExponentsv3UF = []; prfExponentv3ConfLowUF = []; prfExponentv3ConfHighUF = [];
prfExponentsv1v3UF = []; prfExponentv1v3ConfLowUF = []; prfExponentv1v3ConfHighUF = [];

%for each subject, get the exponent of the exponential fit between v1 and v3
for subject = 1:length(subjectList);
    sprintf('Doing subject %i of %i',subject,length(subjectList))
    dataStr =  strcat('data=', subjectList(subject), 'pRFGaussianHdrNMFullUF.mat');
    %get the exponents
    [v1Exponent, v1ExponentConfInt, v2Exponent, v2ExponentConfInt, v3Exponent, v3ExponentConfInt, v1v2Exponent, v1v2ExponentConfInt, v1v3Exponent, v1v3ExponentConfInt] = getpRFTSeries('loadData=1',dataStr{1},'graphStuff=0'); close all;
    prfExponentsv1UF = [prfExponentsv1UF v1Exponent]; prfExponentv1ConfLowUF = [prfExponentv1ConfLowUF v1ExponentConfInt(1)]; prfExponentv1ConfHighUF = [prfExponentv1ConfHighUF v1ExponentConfInt(2)];
    prfExponentsv2UF = [prfExponentsv2UF v2Exponent]; prfExponentv2ConfLowUF = [prfExponentv2ConfLowUF v2ExponentConfInt(1)]; prfExponentv2ConfHighUF = [prfExponentv2ConfHighUF v2ExponentConfInt(2)];
    prfExponentsv3UF = [prfExponentsv3UF v3Exponent]; prfExponentv3ConfLowUF = [prfExponentv3ConfLowUF v3ExponentConfInt(1)]; prfExponentv3ConfHighUF = [prfExponentv3ConfHighUF v3ExponentConfInt(2)];
    prfExponentsv1v3UF = [prfExponentsv1v3UF v1v3Exponent]; prfExponentv1v3ConfLowUF = [prfExponentv1v3ConfLowUF v1v3ExponentConfInt(1)]; prfExponentv1v3ConfHighUF = [prfExponentv1v3ConfHighUF v1v3ExponentConfInt(2)];
end 





%% plot the data %%
color1 = "#0072BD"; color2 = "#77AC30"; color3 = "#EDB120"; color4 = "#7E2F8E";

% plot prf vs mean model comparison
figure, hold on

errorbar(MMExponentsv1,prfExponentsv1,prfExponentv1ConfLow-prfExponentsv1,prfExponentv1ConfHigh-prfExponentsv1,'.','markerSize',14,'LineStyle','none','Color',color1,'CapSize',2,'lineWidth',1) 
errorbar(MMExponentsv1,prfExponentsv1,MMExponentv1ConfLow-MMExponentsv1,MMExponentv1ConfHigh-MMExponentsv1,'.','horizontal','markerSize',14,'LineStyle','none','Color',color1,'CapSize',2,'lineWidth',1) 

errorbar(MMExponentsv2,prfExponentsv2,prfExponentv2ConfLow-prfExponentsv2,prfExponentv2ConfHigh-prfExponentsv2,'.','markerSize',14,'LineStyle','none','Color',color2,'CapSize',2,'lineWidth',1) 
errorbar(MMExponentsv2,prfExponentsv2,MMExponentv2ConfLow-MMExponentsv2,MMExponentv2ConfHigh-MMExponentsv2,'.','horizontal','markerSize',14,'LineStyle','none','Color',color2,'CapSize',2,'lineWidth',1) 

errorbar(MMExponentsv3,prfExponentsv3,prfExponentv3ConfLow-prfExponentsv3,prfExponentv3ConfHigh-prfExponentsv3,'.','markerSize',14,'LineStyle','none','Color',color3,'CapSize',2,'lineWidth',1) 
errorbar(MMExponentsv3,prfExponentsv3,MMExponentv3ConfLow-MMExponentsv3,MMExponentv3ConfHigh-MMExponentsv3,'.','horizontal','markerSize',14,'LineStyle','none','Color',color3,'CapSize',2,'lineWidth',1) 

errorbar(MMExponentsv1v3,prfExponentsv1v3,prfExponentv1v3ConfLow-prfExponentsv1v3,prfExponentv1v3ConfHigh-prfExponentsv1v3,'.','markerSize',14,'LineStyle','none','Color',color4,'CapSize',2,'lineWidth',1) 
errorbar(MMExponentsv1v3,prfExponentsv1v3,MMExponentv1v3ConfLow-MMExponentsv1v3,MMExponentv1v3ConfHigh-MMExponentsv1v3,'.','horizontal','markerSize',14,'LineStyle','none','Color',color4,'CapSize',2,'lineWidth',1) 


maxVal = max([prfExponentsv1v3 prfExponentsv1 prfExponentsv2 prfExponentsv3]); xlim([0 maxVal*1.2]); ylim([0 maxVal*1.2]);
plot([0 maxVal*1.05],[0 maxVal*1.05],'--','Color','k')

xlabel('Mean Model Correlation Fit Exponent');
ylabel('pRF Model Correlation Fit Exponent')
title('Mean model and pRF model Exponents')
legend('v1','','','v2','','','v3','','','v1v3','','','')
%drawPublishAxis('labelFontSize=14','yAxisOffset=-1/25');legend('off');

% plot shufflecor exponents
figure; hold on
errorbar([1:length(shuffleExponentsv1)],shuffleExponentsv1,shuffleExponentv1ConfLow-shuffleExponentsv1,shuffleExponentv1ConfHigh-shuffleExponentsv1,'.','markerSize',14,'LineStyle','none','Color',color1,'CapSize',2)

errorbar([1:length(shuffleExponentsv2)],shuffleExponentsv2,shuffleExponentv2ConfLow-shuffleExponentsv2,shuffleExponentv2ConfHigh-shuffleExponentsv2,'.','markerSize',14,'LineStyle','none','Color',color2,'CapSize',2)

errorbar([1:length(shuffleExponentsv3)],shuffleExponentsv3,shuffleExponentv3ConfLow-shuffleExponentsv3,shuffleExponentv3ConfHigh-shuffleExponentsv3,'.','markerSize',14,'LineStyle','none','Color',color3,'CapSize',2)

errorbar([1:length(shuffleExponentsv1v3)],shuffleExponentsv1v3,shuffleExponentv1v3ConfLow-shuffleExponentsv1v3,shuffleExponentv1v3ConfHigh-shuffleExponentsv1v3,'.','markerSize',14,'LineStyle','none','Color',color4,'CapSize',2)

xlim([0 length(subjectList)+1])
xlabel('Subject');
ylabel('Exponent of correlations between scans')
title('Residual Correlations between scans')
legend('v1','','v2','','v3','','v1v3','')
%drawPublishAxis('labelFontSize=14','yAxisOffset=-1/25');legend('off');

%filtered vs unfiltered data
% plot prf vs mean model comparison
figure, hold on

errorbar(prfExponentsv1UF,prfExponentsv1,prfExponentv1ConfLow-prfExponentsv1,prfExponentv1ConfHigh-prfExponentsv1,'.','markerSize',14,'LineStyle','none','Color',color1,'CapSize',2) 
errorbar(prfExponentsv1UF,prfExponentsv1,prfExponentv1ConfLowUF-prfExponentsv1UF,prfExponentv1ConfHighUF-prfExponentsv1UF,'.','horizontal','markerSize',14,'LineStyle','none','Color',color1,'CapSize',2) 

errorbar(prfExponentsv2UF,prfExponentsv2,prfExponentv2ConfLow-prfExponentsv2,prfExponentv2ConfHigh-prfExponentsv2,'.','markerSize',14,'LineStyle','none','Color',color2,'CapSize',2) 
errorbar(prfExponentsv2UF,prfExponentsv2,prfExponentv2ConfLowUF-prfExponentsv2UF,prfExponentv2ConfHighUF-prfExponentsv2UF,'.','horizontal','markerSize',14,'LineStyle','none','Color',color2,'CapSize',2) 

errorbar(prfExponentsv3UF,prfExponentsv3,prfExponentv3ConfLow-prfExponentsv3,prfExponentv3ConfHigh-prfExponentsv3,'.','markerSize',14,'LineStyle','none','Color',color3,'CapSize',2) 
errorbar(prfExponentsv3UF,prfExponentsv3,prfExponentv3ConfLowUF-prfExponentsv3UF,prfExponentv3ConfHighUF-prfExponentsv3UF,'.','horizontal','markerSize',14,'LineStyle','none','Color',color3,'CapSize',2) 

errorbar(prfExponentsv1v3UF,prfExponentsv1v3,prfExponentv1v3ConfLow-prfExponentsv1v3,prfExponentv1v3ConfHigh-prfExponentsv1v3,'.','markerSize',14,'LineStyle','none','Color',color4,'CapSize',2) 
errorbar(prfExponentsv1v3UF,prfExponentsv1v3,prfExponentv1v3ConfLowUF-prfExponentsv1v3UF,prfExponentv1v3ConfHighUF-prfExponentsv1v3UF,'.','horizontal','markerSize',14,'LineStyle','none','Color',color4,'CapSize',2) 

maxVal = max([prfExponentsv1v3 prfExponentsv1 prfExponentsv2 prfExponentsv3]); xlim([0 maxVal*1.2]); ylim([0 maxVal*1.2]);
plot([0 maxVal*1.05],[0 maxVal*1.05],'--','Color','k')

xlabel('Unfiltered data pRF Model Correlation Fit Exponent');
ylabel('pRF Model Correlation Fit Exponent')
title('Filtered and Unfiltered Fit Exponents')
legend('v1','','','v2','','','v3','','','v1v3','','','')
%drawPublishAxis('labelFontSize=14','yAxisOffset=-1/25');legend('off');



keyboard


