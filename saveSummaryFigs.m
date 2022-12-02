prfExponents = [];
prfExponentConfLow = [];
prfExponentConfHigh = [];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0415mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0416mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0405mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0406mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0419mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0350mc1348GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0420mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0422mc1234GaussianNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0399prf.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0401prf.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0403prf.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0404prf.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];

[v1Exponent expConfInt] = getpRFTSeries('loadData=1','data=s0423mc12345GaussianHdrNM.mat'); close all;
prfExponents = [prfExponents v1Exponent]; prfExponentConfLow = [prfExponentConfLow expConfInt(1)]; prfExponentConfHigh = [prfExponentConfHigh expConfInt(2)];


% do the mean model
MMExponents = [];
MMExponentConfLow = [];
MMExponentConfHigh = [];

[v1Exponent expConfInt] = avgModelNoise('data1=s0415mc12GaussianNM.mat','data2=s0415mc34GaussianNM.mat','base=s0415mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0416mc12GaussianNM.mat','data2=s0416mc34GaussianNM.mat','base=s0416mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0405mc12GaussianNM.mat','data2=s0405mc34GaussianNM.mat','base=s0405mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0406mc12GaussianNM.mat','data2=s0406mc34GaussianNM.mat','base=s0406mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0419mc12GaussianNM.mat','data2=s0419mc34GaussianNM.mat','base=s0419mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0350mc13GaussianNM.mat','data2=s0350mc48GaussianNM.mat','base=s0350mc1348GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0420mc12GaussianNM.mat','data2=s0420mc34GaussianNM.mat','base=s0420mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0422mc12GaussianNM.mat','data2=s0422mc34GaussianNM.mat','base=s0422mc1234GaussianNM.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0399mc345GaussianHdrNM.mat','data2=s0399mc12GaussianHdrNM.mat','base=s0399pRF.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0401mc345GaussianHdrNM.mat','data2=s0401mc12GaussianHdrNM.mat','base=s0401pRF.mat')
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0403mc345GaussianHdrNM.mat','data2=s0403mc12GaussianHdrNM.mat','base=s0403pRF.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0404mc345GaussianHdrNM.mat','data2=s0404mc12GaussianHdrNM.mat','base=s0404pRF.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all

[v1Exponent expConfInt] = avgModelNoise('data1=s0423mc678GaussianHdrNM.mat','data2=s0423mc12345GaussianHdrNM.mat','base=s0423pRF.mat');
MMExponents = [MMExponents v1Exponent]; MMExponentConfLow = [MMExponentConfLow expConfInt(1)]; MMExponentConfHigh = [MMExponentConfHigh expConfInt(2)]; close all






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








%% shufflecor %%

shuffleV1Exponents = []; shuffleV1ExponentConfLow = []; shuffleV1ExponentConfHigh = [];
shuffleV2Exponents = []; shuffleV2ExponentConfLow = []; shuffleV2ExponentConfHigh = [];
shuffleV3Exponents = []; shuffleV3ExponentConfLow = []; shuffleV3ExponentConfHigh = [];

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0415mc12GaussianNM.mat','data2=s0415mc34GaussianNM.mat','base=s0415mc1234GaussianNM.mat');
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0416mc12GaussianNM.mat','data2=s0416mc34GaussianNM.mat','base=s0416mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0405mc12GaussianNM.mat','data2=s0405mc34GaussianNM.mat','base=s0405mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0406mc12GaussianNM.mat','data2=s0406mc34GaussianNM.mat','base=s0406mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0419mc12GaussianNM.mat','data2=s0419mc34GaussianNM.mat','base=s0419mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0350mc13GaussianNM.mat','data2=s0350mc48GaussianNM.mat','base=s0350mc1348GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0420mc12GaussianNM.mat','data2=s0420mc34GaussianNM.mat','base=s0420mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0422mc12GaussianNM.mat','data2=s0422mc34GaussianNM.mat','base=s0422mc1234GaussianNM.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0399mc12GaussianNM','data2=s0399mc345GaussianNM','base=s0399pRF.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0401mc12GaussianNM','data2=s0401mc345GaussianNM','base=s0401pRF.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0403mc12GaussianNM','data2=s0403mc345GaussianNM','base=s0403pRF.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0404mc12GaussianNM','data2=s0404mc345GaussianNM','base=s0404pRF.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all

[v1Exponent v1ExponentConfInt v2Exponent v2ExponentConfInt v3Exponent v3ExponentConfInt] = shufflecor('data1=s0423mc678GaussianNM','data2=s0423mc12345GaussianNM','base=s0423pRF.mat')
shuffleV1Exponents = [shuffleV1Exponents v1Exponent]; shuffleV1ExponentConfLow = [shuffleV1ExponentConfLow v1ExponentConfInt(1)]; shuffleV1ExponentConfHigh = [shuffleV1ExponentConfHigh v1ExponentConfInt(2)]; close all
shuffleV2Exponents = [shuffleV2Exponents v2Exponent]; shuffleV2ExponentConfLow = [shuffleV2ExponentConfLow v2ExponentConfInt(1)]; shuffleV2ExponentConfHigh = [shuffleV2ExponentConfHigh v2ExponentConfInt(2)]; close all
shuffleV3Exponents = [shuffleV3Exponents v3Exponent]; shuffleV3ExponentConfLow = [shuffleV3ExponentConfLow v3ExponentConfInt(1)]; shuffleV3ExponentConfHigh = [shuffleV3ExponentConfHigh v3ExponentConfInt(2)]; close all




figure; hold on
errorbar([1:13],shuffleV1Exponents,shuffleV1ExponentConfLow-shuffleV1Exponents,shuffleV1ExponentConfHigh-shuffleV1Exponents,'LineStyle','none','Color','k','CapSize',2)
errorbar([1:13],shuffleV2Exponents,shuffleV2ExponentConfLow-shuffleV2Exponents,shuffleV2ExponentConfHigh-shuffleV2Exponents,'LineStyle','none','Color','r','CapSize',2)
errorbar([1:13],shuffleV3Exponents,shuffleV3ExponentConfLow-shuffleV3Exponents,shuffleV3ExponentConfHigh-shuffleV3Exponents,'LineStyle','none','Color','b','CapSize',2)



keyboard

