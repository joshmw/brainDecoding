%%save figs%%

function savefigs

getpRFTSeries('loadData=1','data=s0415mc1234GaussianNM.mat');
savepdf(figure(11),'s0415v1prf')
savepdf(figure(12),'s0415v2prf')
savepdf(figure(13),'s0415v3prf')
savepdf(figure(15),'s0415v1v3prf')
close all;

getpRFTSeries('loadData=1','data=s0416mc1234GaussianNM.mat');
savepdf(figure(11),'s0416v1prf')
savepdf(figure(12),'s0416v2prf')
savepdf(figure(13),'s0416v3prf')
savepdf(figure(15),'s0416v1v3prf')
close all;



getpRFTSeries('loadData=1','data=s0405mc1234GaussianNM.mat');
savepdf(figure(11),'s0405v1prf')
savepdf(figure(12),'s0405v2prf')
savepdf(figure(13),'s0405v3prf')
savepdf(figure(15),'s0405v1v3prf')
close all;

getpRFTSeries('loadData=1','data=s0406mc1234GaussianNM.mat');
savepdf(figure(11),'s0406v1prf')
savepdf(figure(12),'s0406v2prf')
savepdf(figure(13),'s0406v3prf')
savepdf(figure(15),'s0406v1v3prf')
close all;

getpRFTSeries('loadData=1','data=s0419mc1234GaussianNM.mat');
savepdf(figure(11),'s0419v1prf')
savepdf(figure(12),'s0419v2prf')
savepdf(figure(13),'s0419v3prf')
savepdf(figure(15),'s0419v1v3prf')
close all;


getpRFTSeries('loadData=1','data=s0350mc1348GaussianNM.mat');
savepdf(figure(11),'s0350v1prf')
savepdf(figure(12),'s0350v2prf')
savepdf(figure(13),'s0350v3prf')
savepdf(figure(15),'s0350v1v3prf')
close all;


getpRFTSeries('loadData=1','data=s0420mc1234GaussianNM.mat');
savepdf(figure(11),'s0420v1prf')
savepdf(figure(12),'s0420v2prf')
savepdf(figure(13),'s0420v3prf')
savepdf(figure(15),'s0420v1v3prf')
close all;


getpRFTSeries('loadData=1','data=s0422mc1234GaussianNM.mat');
savepdf(figure(11),'s0422v1prf')
savepdf(figure(12),'s0422v2prf')
savepdf(figure(13),'s0422v3prf')
savepdf(figure(15),'s0422v1v3prf')
close all;






getpRFTSeries('loadData=1','data=s0399prf.mat');
savepdf(figure(11),'s0399v1prf')
savepdf(figure(12),'s0399v2prf')
savepdf(figure(13),'s0399v3prf')
savepdf(figure(15),'s0399v1v3prf')
close all;

getpRFTSeries('loadData=1','data=s0401prf.mat');
savepdf(figure(11),'s0401v1prf')
savepdf(figure(12),'s0401v2prf')
savepdf(figure(13),'s0401v3prf')
savepdf(figure(15),'s0401v1v3prf')
close all

getpRFTSeries('loadData=1','data=s0404prf.mat');
savepdf(figure(11),'s0404v1prf')
savepdf(figure(12),'s0404v2prf')
savepdf(figure(13),'s0404v3prf')
savepdf(figure(15),'s0404v1v3prf')
close all


getpRFTSeries('loadData=1','data=s0403prf.mat');
savepdf(figure(11),'s0403v1prf')
savepdf(figure(12),'s0403v2prf')
savepdf(figure(13),'s0403v3prf')
savepdf(figure(15),'s0403v1v3prf')
close all

getpRFTSeries('loadData=1','data=s0423mc12345GaussianHdrNM.mat');
savepdf(figure(11),'s0423v1prf')
savepdf(figure(12),'s0423v2prf')
savepdf(figure(13),'s0423v3prf')
savepdf(figure(15),'s0423v1v3prf')
close all



% avgmodelnoise

avgModelNoise('data1=s0416mc12GaussianNM.mat','data2=s0416mc34GaussianNM.mat','base=s0416mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0416v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0416v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0416v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0416v1v3MeanModel')
close all


avgModelNoise('data1=s0415mc12GaussianNM.mat','data2=s0415mc34GaussianNM.mat','base=s0415mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0415v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0415v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0415v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0415v1v3MeanModel')
close all


avgModelNoise('data1=s0405mc12GaussianNM.mat','data2=s0405mc34GaussianNM.mat','base=s0405mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0405v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0405v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0405v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0405v1v3MeanModel')
close all

avgModelNoise('data1=s0406mc12GaussianNM.mat','data2=s0406mc34GaussianNM.mat','base=s0406mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0406v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0406v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0406v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0406v1v3MeanModel')
close all

avgModelNoise('data1=s0419mc12GaussianNM.mat','data2=s0419mc34GaussianNM.mat','base=s0419mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0419v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0419v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0419v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0419v1v3MeanModel')
close all

avgModelNoise('data1=s0350mc13GaussianNM.mat','data2=s0350mc48GaussianNM.mat','base=s0350mc1348GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0350v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0350v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0350v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0350v1v3MeanModel')
close all

avgModelNoise('data1=s0420mc12GaussianNM.mat','data2=s0420mc34GaussianNM.mat','base=s0420mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0420v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0420v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0420v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0420v1v3MeanModel')
close all

avgModelNoise('data1=s0422mc12GaussianNM.mat','data2=s0422mc34GaussianNM.mat','base=s0422mc1234GaussianNM.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0422v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0422v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0422v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0422v1v3MeanModel')
close all





avgModelNoise('data1=s0399mc345GaussianHdrNM.mat','data2=s0399mc12GaussianHdrNM.mat','base=s0399pRF.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0399v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0399v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0399v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0399v1v3MeanModel')
close all


avgModelNoise('data1=s0401mc345GaussianHdrNM.mat','data2=s0401mc12GaussianHdrNM.mat','base=s0401pRF.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0401v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0401v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0401v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0401v1v3MeanModel')
close all

avgModelNoise('data1=s0403mc345GaussianHdrNM.mat','data2=s0403mc12GaussianHdrNM.mat','base=s0403pRF.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0403v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0403v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0403v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0403v1v3MeanModel')
close all


avgModelNoise('data1=s0404mc345GaussianHdrNM.mat','data2=s0404mc12GaussianHdrNM.mat','base=s0404pRF.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0404v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0404v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0404v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0404v1v3MeanModel')
close all


avgModelNoise('data1=s0423mc678GaussianHdrNM.mat','data2=s0423mc12345GaussianHdrNM.mat','base=s0423pRF.mat')
figure(11);set(gcf,'renderer','Painters'),savepdf(figure(11),'s0423v1MeanModel')
figure(12);set(gcf,'renderer','Painters'),savepdf(figure(12),'s0423v2MeanModel')
figure(13);set(gcf,'renderer','Painters'),savepdf(figure(13),'s0423v3MeanModel')
figure(15);set(gcf,'renderer','Painters'),savepdf(figure(15),'s0423v1v3MeanModel')
close all


%shufflecor

shufflecor('data1=s0416mc12GaussianNM.mat','data2=s0416mc34GaussianNM.mat','base=s0416mc1234GaussianNM.mat')
savepdf(figure(11),'s0416shuffleCorV1')
savepdf(figure(12),'s0416shuffleCorV2')
savepdf(figure(13),'s0416shuffleCorV3')
close all

shufflecor('data1=s0415mc12GaussianNM.mat','data2=s0415mc34GaussianNM.mat','base=s0415mc1234GaussianNM.mat')
savepdf(figure(11),'s0415shuffleCorV1')
savepdf(figure(12),'s0415shuffleCorV2')
savepdf(figure(13),'s0415shuffleCorV3')
close all


shufflecor('data1=s0405mc12GaussianNM.mat','data2=s0405mc34GaussianNM.mat','base=s0405mc1234GaussianNM.mat')
savepdf(figure(11),'s0405shuffleCorV1')
savepdf(figure(12),'s0405shuffleCorV2')
savepdf(figure(13),'s0405shuffleCorV3')
close all

shufflecor('data1=s0406mc12GaussianNM.mat','data2=s0406mc34GaussianNM.mat','base=s0406mc1234GaussianNM.mat')
savepdf(figure(11),'s0406shuffleCorV1')
savepdf(figure(12),'s0406shuffleCorV2')
savepdf(figure(13),'s0406shuffleCorV3')
close all

shufflecor('data1=s0419mc12GaussianNM.mat','data2=s0419mc34GaussianNM.mat','base=s0419mc1234GaussianNM.mat')
savepdf(figure(11),'s0419shuffleCorV1')
savepdf(figure(12),'s0419shuffleCorV2')
savepdf(figure(13),'s0419shuffleCorV3')
close all

shufflecor('data1=s0350mc13GaussianNM.mat','data2=s0350mc48GaussianNM.mat','base=s0350mc1348GaussianNM.mat')
savepdf(figure(11),'s0350shuffleCorV1')
savepdf(figure(12),'s0350shuffleCorV2')
savepdf(figure(13),'s0350shuffleCorV3')
close all

shufflecor('data1=s0420mc12GaussianNM.mat','data2=s0420mc34GaussianNM.mat','base=s0420mc1234GaussianNM.mat')
savepdf(figure(11),'s0420shuffleCorV1')
savepdf(figure(12),'s0420shuffleCorV2')
savepdf(figure(13),'s0420shuffleCorV3')
close all

shufflecor('data1=s0422mc12GaussianNM.mat','data2=s0422mc34GaussianNM.mat','base=s0422mc1234GaussianNM.mat')
savepdf(figure(11),'s0422shuffleCorV1')
savepdf(figure(12),'s0422shuffleCorV2')
savepdf(figure(13),'s0422shuffleCorV3')
close all



shufflecor('data1=s0399mc12GaussianNM','data2=s0399mc345GaussianNM','base=s0399pRF.mat')
savepdf(figure(11),'s0399shuffleCorV1')
savepdf(figure(12),'s0399shuffleCorV2')
savepdf(figure(13),'s0399shuffleCorV3')
close all


shufflecor('data1=s0401mc12GaussianNM','data2=s0401mc345GaussianNM','base=s0401pRF.mat')
savepdf(figure(11),'s0401shuffleCorV1')
savepdf(figure(12),'s0401shuffleCorV2')
savepdf(figure(13),'s0401shuffleCorV3')
close all


shufflecor('data1=s0403mc12GaussianNM','data2=s0403mc345GaussianNM','base=s0403pRF.mat')
savepdf(figure(11),'s0403shuffleCorV1')
savepdf(figure(12),'s0403shuffleCorV2')
savepdf(figure(13),'s0403shuffleCorV3')
close all


shufflecor('data1=s0404mc12GaussianNM','data2=s0404mc345GaussianNM','base=s0404pRF.mat')
savepdf(figure(11),'s0404shuffleCorV1')
savepdf(figure(12),'s0404shuffleCorV2')
savepdf(figure(13),'s0404shuffleCorV3')
close all


shufflecor('data1=s0423mc678GaussianNM','data2=s0423mc12345GaussianNM','base=s0423pRF.mat')
savepdf(figure(11),'s0423shuffleCorV1')
savepdf(figure(12),'s0423shuffleCorV2')
savepdf(figure(13),'s0423shuffleCorV3')
close all












% multiplicative figures

[additiveFit scanA scanB] = avgModelNoise('data1=s0399mc345GaussianHdrNM.mat','data2=s0399mc12GaussianHdrNM.mat','base=s0399pRF.mat')
savepdf(figure(1),'s0399multiCurvesZscores')
savepdf(figure(2),'s0399multiHistZscores')
save('s0399graphFitsZscores','additiveFit','scanA','scanB')
close all

[additiveFit scanA scanB] = avgModelNoise('data1=s0401mc345GaussianHdrNM.mat','data2=s0401mc12GaussianHdrNM.mat','base=s0401pRF.mat')
savepdf(figure(1),'s0401multiCurvesZscores')
savepdf(figure(2),'s0401multiHistZscores')
save('s0401graphFitsZscores','additiveFit','scanA','scanB')
close all

[additiveFit scanA scanB] = avgModelNoise('data1=s0403mc345GaussianHdrNM.mat','data2=s0403mc12GaussianHdrNM.mat','base=s0403pRF.mat')
savepdf(figure(1),'s0403multiCurvesZscores')
savepdf(figure(2),'s0403multiHistZscores')
save('s0403graphFitsZscores','additiveFit','scanA','scanB')
close all

[additiveFit scanA scanB] = avgModelNoise('data1=s0404mc345GaussianHdrNM.mat','data2=s0404mc12GaussianHdrNM.mat','base=s0404pRF.mat')
savepdf(figure(1),'s0404multiCurvesZscores')
savepdf(figure(2),'s0404multiHistZscores')
save('s0404graphFitsZscores','additiveFit','scanA','scanB')
close all

[additiveFit scanA scanB] = avgModelNoise('data1=s0423mc678GaussianHdrNM.mat','data2=s0423mc12345GaussianHdrNM.mat','base=s0423pRF.mat')
savepdf(figure(1),'s0423multiCurvesZscores')
savepdf(figure(2),'s0423multiHistZscores')
save('s0423graphFitsZscores','additiveFit','scanA','scanB')
close all


















