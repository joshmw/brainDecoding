function graphCor(filename)

load(filename)

figure(9);
for roi = 1:length(cleanRois)
    x = []; y = []; avg = [];
for bin = 0:.1:.9
    z = []; avgs = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    if mean(cleanRois(roi).vox.pRFtSeries(:,voxel)) > .5 %this fixes an issue where some model responses are centered around 0
    z = [z std(cleanRois(roi).vox.baselineNoise(quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin+.025)>cleanRois(roi).vox.pRFtSeries(:,voxel)&cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin),voxel))]; 
    avgs = [avgs mean(cleanRois(roi).vox.baselineNoise(quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin+.025)>cleanRois(roi).vox.pRFtSeries(:,voxel)&cleanRois(roi).vox.pRFtSeries(:,voxel)>quantile(cleanRois(roi).vox.pRFtSeries(:,voxel),bin),voxel))];
    end
end
x = [x bin];
y = [y mean(z)];
avg = [avg mean(avgs)];
end
    
subplot(2,length(cleanRois),roi); hold on;
scatter(x,y*100,'black');
title(sprintf('%s Residual Std by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Noise (% signal change)');

subplot(2,length(cleanRois),roi+3); hold on;
scatter(x,avg*100,'black');plot([0,1],[0,0],'black--');
title(sprintf('%s Average Residual by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Average Residual (% signal change)');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);subplot(1,2,1);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1kldCor(i,:),v1NoiseCor(i,:)); end;

v1kldCorArr = reshape(v1kldCor,[1 length(v1kldCor)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1kldCorArr,sortOrder] = sort(v1kldCorArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1kldCorArr) & (v1kldCorArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1kldCorArr) & (v1kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V1 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxel i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');

% second, do v2 %
figure(12);hold on;subplot(1,2,1);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2kldCor(i,:),v2NoiseCor(i,:)); end;

v2kldCorArr = reshape(v2kldCor,[1 length(v2kldCor)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2kldCorArr,sortOrder] = sort(v2kldCorArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldCorArr) & (v2kldCorArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2kldCorArr) & (v2kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V2 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


% third, do v3 %
figure(13);hold on;subplot(1,2,1);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3kldCor(i,:),v3NoiseCor(i,:)); end;

v3kldCorArr = reshape(v3kldCor,[1 length(v3kldCor)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3kldCorArr,sortOrder] = sort(v3kldCorArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldCorArr) & (v3kldCorArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3kldCorArr) & (v3kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,20]);set(gca,'XDir','reverse');
title('V3 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


%between v1 and v2
figure(14);hold on;subplot(1,2,1);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2kldCor(i,:),v1v2NoiseCor(i,:)); end;

v1v2kldCorArr = reshape(v1v2kldCor,[1 min(size(v1v2kldCor))*max(size(v1v2kldCor))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2kldCorArr,sortOrder] = sort(v1v2kldCorArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2kldCorArr) & (v1v2kldCorArr < bin+step) ) > 15
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2kldCorArr) & (v1v2kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,20]);set(gca,'XDir','reverse');
title('Receptive Field and Noise Correlations between v1 and v2'); xlabel('KL Divergence between voxel V1i, V2j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V2j (mm)');


%between v1 and v3
figure(15);hold on;subplot(1,2,1);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3kldCor(i,:),v1v3NoiseCor(i,:)); end;

v1v3kldCorArr = reshape(v1v3kldCor,[1 min(size(v1v3kldCor))*max(size(v1v3kldCor))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3kldCorArr,sortOrder] = sort(v1v3kldCorArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3kldCorArr) & (v1v3kldCorArr < bin+step) ) > 15
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3kldCorArr) & (v1v3kldCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,20]);set(gca,'XDir','reverse');
title('Receptive Field and Noise Correlations between v1 and v3'); xlabel('KL Divergence between voxel V1i, V3j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V3j (mm)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distance and noise correlation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);hold on;subplot(1,2,2);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1dist(i,:),v1NoiseCor(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);set(gca,'XDir','reverse');
title('V1 Distance and Noise Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');

% second, do v2 %
figure(12);hold on;subplot(1,2,2);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2dist(i,:),v2NoiseCor(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);set(gca,'XDir','reverse');
title('V2 Distance and Noise Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');


% third, do v3 %
figure(13);hold on;subplot(1,2,2);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3dist(i,:),v3NoiseCor(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('v3 Distance and Noise correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('Noise correlation between voxels i,j');


%between v1 and v2
figure(14);hold on;subplot(1,2,2);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2dist(i,:),v1v2NoiseCor(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 15
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('Distance and Noise Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('Noise correlation between voxels V1i, V2j');


%between v1 and v3
figure(15);hold on;subplot(1,2,2);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3dist(i,:),v1v3NoiseCor(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 15
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); set(gca,'XDir','reverse');
title('Distance and Noise Correlations between v1 and v3'); xlabel('Distance between voxels V1i, V3j'); ylabel('Noise correlation between voxels V1i, V2j');






%%%%%%%%%%%%%%%%%%%%%%%%
%% distance vs rf cor %%
%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(25);hold on;subplot(2,3,1);hold on;

for i = 1:length(v1kldCor); scatter(v1dist(i,:),v1kldCor(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1kldCorArr = reshape(v1kldCor,[1 length(v1kldCor)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1kldCorArr = v1kldCorArr(sortOrder);
bins = [];kldCorAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 15
        kldCorAvg = median(v1kldCorArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V1 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v2 %
figure(25);hold on;subplot(2,3,2);hold on;

for i = 1:length(v2kldCor); scatter(v2dist(i,:),v2kldCor(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2kldCorArr = reshape(v2kldCor,[1 length(v2kldCor)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2kldCorArr = v2kldCorArr(sortOrder);
bins = [];kldCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 15
        kldCorAvg = median(v2kldCorArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V2 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v3 %
figure(25);hold on;subplot(2,3,3);hold on;

for i = 1:length(v3kldCor); scatter(v3dist(i,:),v3kldCor(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3kldCorArr = reshape(v3kldCor,[1 length(v3kldCor)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3kldCorArr = v3kldCorArr(sortOrder);
bins = [];kldCorAvgs = []; 
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 15
        kldCorAvg = median(v3kldCorArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end
plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('V3 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


%between v1 and v2
figure(25);hold on;subplot(2,3,4);hold on;

for i = 1:min(size(v1v2kldCor)); scatter(v1v2dist(i,:),v1v2kldCor(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2kldCorArr = reshape(v1v2kldCor,[1 min(size(v1v2kldCor))*max(size(v1v2kldCor))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2kldCorArr = v1v2kldCorArr(sortOrder);

bins = [];kldCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 15
        kldCorAvg = median(v1v2kldCorArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end

plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('Distance and RF Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('KL divergence between voxel RFs V1i, V2j (distance)');


%between v1 and v3
figure(25);hold on;subplot(2,3,5);hold on;

for i = 1:min(size(v1v3kldCor)); scatter(v1v3dist(i,:),v1v3kldCor(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3kldCorArr = reshape(v1v3kldCor,[1 min(size(v1v3kldCor))*max(size(v1v3kldCor))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3kldCorArr = v1v3kldCorArr(sortOrder);

bins = [];kldCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 15
        kldCorAvg = median(v1v3kldCorArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; kldCorAvgs = [kldCorAvgs kldCorAvg];
    end
end

plot(bins,kldCorAvgs,'black','LineWidth',8);ylim([0,25]);
title('Distance and RF Correlations between v1 and v3'); xlabel('Distance between voxels V1i, V3j'); ylabel('KL divergence between voxel RFs RFs V1i, V2j (distance)');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time series and noise correlation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% V1
figure(26);hold on;subplot(1,3,1);hold on;

v1tSeriesCor = corrcoef(cleanRois(1).vox.pRFtSeries);

for i = 1:length(v1NoiseCor); scatter(v1tSeriesCor(i,:),v1NoiseCor(i,:)); end;

v1tSeriesCorArr = reshape(v1tSeriesCor,[1 length(v1tSeriesCor)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1tSeriesCorArr,sortOrder] = sort(v1tSeriesCorArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; step = .05;
for bin = -.5:step:1
    if sum( (bin < v1tSeriesCorArr) & (v1tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1tSeriesCorArr) & (v1tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V1 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');


%% V2
figure(26);hold on;subplot(1,3,2);hold on;

v2tSeriesCor = corrcoef(cleanRois(2).vox.pRFtSeries);

for i = 1:length(v2NoiseCor); scatter(v2tSeriesCor(i,:),v2NoiseCor(i,:)); end;

v2tSeriesCorArr = reshape(v2tSeriesCor,[1 length(v2tSeriesCor)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2tSeriesCorArr,sortOrder] = sort(v2tSeriesCorArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = []; 
for bin = -.5:step:1
    if sum( (bin < v2tSeriesCorArr) & (v2tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2tSeriesCorArr) & (v2tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V2 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');


%% V3
figure(26);hold on;subplot(1,3,3);hold on;

v3tSeriesCor = corrcoef(cleanRois(3).vox.pRFtSeries);

for i = 1:length(v3NoiseCor); scatter(v3tSeriesCor(i,:),v3NoiseCor(i,:)); end;

v3tSeriesCorArr = reshape(v3tSeriesCor,[1 length(v3tSeriesCor)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3tSeriesCorArr,sortOrder] = sort(v3tSeriesCorArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = -.5:step:1
    if sum( (bin < v3tSeriesCorArr) & (v3tSeriesCorArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3tSeriesCorArr) & (v3tSeriesCorArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([-.5,1]);ylim([-.5,1]);
title('V3 Time Series and Noise Correlation'); xlabel('Time Series Correlation between voxels i,j'); ylabel('Noise correlation between voxels i,j');

