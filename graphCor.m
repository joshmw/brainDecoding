function graphCor(filename)

load(filename)
   


% receptive field overlaps % delete this later once i save it in the environment files
%%%%%%%%%%%%%%%%
%% RF Overlap %%
%%%%%%%%%%%%%%%%

% receptive field overlaps %
sprintf('Correlation RFs and calculating distances (takes about a minute)...')

%v1
roi = 1; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)   
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v1rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end

% v2 %
roi = 2; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)      
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v2rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v2rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end


roi = 3; for row = 1:length(cleanRois(roi).vox.linearCoords)
    for column = 1:length(cleanRois(roi).vox.linearCoords)       
        mu1 = 0; s1 = cleanRois(roi).vox.rfstd(column); s2 = cleanRois(roi).vox.rfstd(row);
        mu2 = sqrt((cleanRois(roi).vox.x(column)-cleanRois(roi).vox.x(row))^2 + (cleanRois(roi).vox.y(column)-cleanRois(roi).vox.y(row))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(row,column) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v3rfOverlap(row,column) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end




% v1/v2 
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(2).vox.linearCoords)       
        mu1 = 0; s1 = cleanRois(1).vox.rfstd(roi1vox); s2 = cleanRois(2).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRois(1).vox.x(roi1vox)-cleanRois(2).vox.x(roi2vox))^2 + (cleanRois(1).vox.y(roi1vox)-cleanRois(2).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1v2rfOverlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
    v1v2rfOverlap = transpose(v1v2rfOverlap);

% v1/v3
for roi1vox = 1:length(cleanRois(1).vox.linearCoords)
    for roi2vox = 1:length(cleanRois(3).vox.linearCoords)
        mu1 = 0; s1 = cleanRois(1).vox.rfstd(roi1vox); s2 = cleanRois(3).vox.rfstd(roi2vox);
        mu2 = sqrt((cleanRois(1).vox.x(roi1vox)-cleanRois(3).vox.x(roi2vox))^2 + (cleanRois(1).vox.y(roi1vox)-cleanRois(3).vox.y(roi2vox))^2);
        if mu1 == mu2 & s1 == s2; v3rfOverlap(roi1vox,roi2vox) = 1;else;
        c = (mu2*(s1^2) - s2*(mu1*s2 + s1 * sqrt((mu1 - mu2)^2 + 2*(s1^2 - s2^2)*log(s1/s2))))/(s1^2-s2^2);
        v1v3rfOverlap(roi1vox,roi2vox) = 1 - normcdf(c,mu1,s1) + normcdf(c,mu2,s2); end;
    end
end
    v1v3rfOverlap = transpose(v1v3rfOverlap);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% graph %%
%%%%%%%%%%%
    
   
sprintf('Graphing things (also takes about a minute)...')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise std in signal absent and present time frames %%

figure(9);
for roi = 1:length(cleanRois)
    x = []; y = []; avg = []; errhigh = []; errlow = []; errmeanhigh = []; errmeanlow = []; yAll = []; avgAll = [];
seg = .1; for bin = 0:seg:(1-seg)
    z = []; avgs = [];
for voxel = 1:length(cleanRois(roi).vox.linearCoords)
    if mean(cleanRois(roi).vox.pRFtSeries(:,voxel)) > .5 %this fixes an issue where some model responses are centered around 0
    lowbound = min(cleanRois(roi).vox.pRFtSeries(:,voxel))+bin*(max(cleanRois(roi).vox.pRFtSeries(:,voxel))-min(cleanRois(roi).vox.pRFtSeries(:,voxel)));    
    highbound = min(cleanRois(roi).vox.pRFtSeries(:,voxel))+(bin+seg)*(max(cleanRois(roi).vox.pRFtSeries(:,voxel))-min(cleanRois(roi).vox.pRFtSeries(:,voxel)));
    z = [z std(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & cleanRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    avgs = [avgs mean(cleanRois(roi).vox.baselineNoise(cleanRois(roi).vox.pRFtSeries(:,voxel)>=lowbound & cleanRois(roi).vox.pRFtSeries(:,voxel)<=highbound ,voxel))]; 
    end
end

x = [x bin];
y = [y mean(z)];
yAll{length(x)} = z;
avg = [avg mean(avgs)];
avgAll{length(x)} = avgs;

% errorbars %
for i = 1:1000; err(i) = mean(datasample(z,length(cleanRois(roi).vox.linearCoords)));end; err = sort(err);
errhigh = [errhigh quantile(err,.95)-mean(z)];
errlow = [errlow mean(z)-quantile(err,.05)];

for i = 1:1000; err(i) = mean(datasample(avgs,length(cleanRois(roi).vox.linearCoords)));end; err = sort(err);
errmeanhigh = [errmeanhigh quantile(err,.95)-mean(avgs)];
errmeanlow = [errmeanlow mean(avgs)-quantile(err,.05)];
end

% plot %
subplot(2,length(cleanRois),roi); hold on;
%for i = 1:length(x); scatter(repmat(x(i),1,length(cleanRois(roi).vox.linearCoords)),yAll{i}*100); hold on; end;
scatter(x,y*100,60,'black','filled'); 
eb = errorbar(x,y*100,errlow*100,errhigh*100,'o'); eb.Color = 'black';

title(sprintf('%s Residual Std by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Noise (% signal change)');

subplot(2,length(cleanRois),roi+3); hold on;
%for i = 1:length(x); scatter(repmat(x(i),1,length(cleanRois(roi).vox.linearCoords)),avgAll{i}*100); hold on; end;
scatter(x,avg*100,60,'black','filled');plot([0,1],[0,0],'black--');
eb = errorbar(x,avg*100,errmeanlow*100,errmeanhigh*100,'o'); eb.Color = 'black';


title(sprintf('%s Average Residual by Activity Quantile',cleanRois(roi).name))
xlabel('Quantile of Voxel Activity');
ylabel('Average Residual (% signal change)');

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field overlaps cutoff by distance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distanceHigh = 10; distanceLow = 3;

% v1 %
figure(30); v1NoiseArrClose = v1NoiseCor(v1dist<distanceLow & v1dist>0); v1NoiseArrFar = v1NoiseCor(v1dist>distanceHigh);
v1kldArrClose = v1rfOverlap(v1dist<distanceLow & v1dist>0); v1kldArrFar = v1rfOverlap(v1dist>distanceHigh);

subplot(1,2,1);scatter(v1kldArrClose,v1NoiseArrClose);hold on;
bins = [];noiseCorAvgs = []; step = .01;
for bin = 0:step:40
    if sum( (bin < v1kldArrClose) & (v1kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v1NoiseArrClose((bin < v1kldArrClose) & (v1kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('V1 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v1kldArrFar,v1NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v1kldArrFar) & (v1kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v1NoiseArrFar((bin < v1kldArrFar) & (v1kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('V1 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


% v2 %
figure(31); v2NoiseArrClose = v2NoiseCor(v2dist<distanceLow & v2dist>0); v2NoiseArrFar = v2NoiseCor(v2dist>distanceHigh);
v2kldArrClose = v2rfOverlap(v2dist<distanceLow & v2dist>0); v2kldArrFar = v2rfOverlap(v2dist>distanceHigh);

subplot(1,2,1);scatter(v2kldArrClose,v2NoiseArrClose);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldArrClose) & (v2kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v2NoiseArrClose((bin < v2kldArrClose) & (v2kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v2 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v2kldArrFar,v2NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2kldArrFar) & (v2kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v2NoiseArrFar((bin < v2kldArrFar) & (v2kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v2 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


% v3 %
figure(32); v3NoiseArrClose = v3NoiseCor(v3dist<distanceLow & v3dist>0); v3NoiseArrFar = v3NoiseCor(v3dist>distanceHigh);
v3kldArrClose = v3rfOverlap(v3dist<distanceLow & v3dist>0); v3kldArrFar = v3rfOverlap(v3dist>distanceHigh);

subplot(1,2,1);scatter(v3kldArrClose,v3NoiseArrClose);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldArrClose) & (v3kldArrClose < bin+step) ) > 30
        noiseCorAvg = median(v3NoiseArrClose((bin < v3kldArrClose) & (v3kldArrClose < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v3 Receptive field overlap and Noise Correlations (close voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);


subplot(1,2,2);scatter(v3kldArrFar,v3NoiseArrFar);hold on;
bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3kldArrFar) & (v3kldArrFar < bin+step) ) > 30
        noiseCorAvg = median(v3NoiseArrFar((bin < v3kldArrFar) & (v3kldArrFar < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',5);
title('v3 Receptive field overlap and Noise Correlations (far voxels)'); xlabel('Receptive field overlap (percent)'); ylabel('Noise correlation between voxels i,j'); xlim([0,1]); ylim([-.4,1]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise and receptive field correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, do v1 %
figure(11);subplot(1,2,1);hold on;

v1NoiseCor = corrcoef(cleanRois(1).vox.baselineNoise);

for i = 1:length(v1NoiseCor); scatter(v1rfOverlap(i,:),v1NoiseCor(i,:)); end;

v1rfOverlapArr = reshape(v1rfOverlap,[1 length(v1rfOverlap)^2]); v1NoiseCorArr = reshape(v1NoiseCor,[1 length(v1NoiseCor)^2]);
[v1rfOverlapArr,sortOrder] = sort(v1rfOverlapArr); v1NoiseCorArr = v1NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v1rfOverlapArr) & (v1rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1NoiseCorArr((bin < v1rfOverlapArr) & (v1rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V1 Receptive Field and Noise Correlations'); xlabel('Receptive field overlap between voxels i,j (percent)'); ylabel('Noise correlation between voxels i,j (mm)');

% second, do v2 %
figure(12);hold on;subplot(1,2,1);hold on;

v2NoiseCor = corrcoef(cleanRois(2).vox.baselineNoise);   

for i = 1:length(v2NoiseCor); scatter(v2rfOverlap(i,:),v2NoiseCor(i,:)); end;

v2rfOverlapArr = reshape(v2rfOverlap,[1 length(v2rfOverlap)^2]); v2NoiseCorArr = reshape(v2NoiseCor,[1 length(v2NoiseCor)^2]);
[v2rfOverlapArr,sortOrder] = sort(v2rfOverlapArr); v2NoiseCorArr = v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2rfOverlapArr) & (v2rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v2NoiseCorArr((bin < v2rfOverlapArr) & (v2rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V2 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


% third, do v3 %
figure(13);hold on;subplot(1,2,1);hold on;

v3NoiseCor = corrcoef(cleanRois(3).vox.baselineNoise);   

for i = 1:length(v3NoiseCor); scatter(v3rfOverlap(i,:),v3NoiseCor(i,:)); end;

v3rfOverlapArr = reshape(v3rfOverlap,[1 length(v3rfOverlap)^2]); v3NoiseCorArr = reshape(v3NoiseCor,[1 length(v3NoiseCor)^2]);
[v3rfOverlapArr,sortOrder] = sort(v3rfOverlapArr); v3NoiseCorArr = v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:40
    if sum( (bin < v3rfOverlapArr) & (v3rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v3NoiseCorArr((bin < v3rfOverlapArr) & (v3rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8);xlim([0,1]);
title('V3 Receptive Field and Noise Correlations'); xlabel('KL Divergence between voxels i,j RFs (distance)'); ylabel('Noise correlation between voxels i,j (mm)');


%between v1 and v2
figure(14);hold on;subplot(1,2,1);hold on;

v1v2NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(2).vox.baselineNoise));        

for i = 1:min(size(v1v2NoiseCor)); scatter(v1v2rfOverlap(i,:),v1v2NoiseCor(i,:)); end;

v1v2rfOverlapArr = reshape(v1v2rfOverlap,[1 min(size(v1v2rfOverlap))*max(size(v1v2rfOverlap))]);
v1v2NoiseCorArr = reshape(v1v2NoiseCor,[1 min(size(v1v2NoiseCor))*max(size(v1v2NoiseCor))]);
[v1v2rfOverlapArr,sortOrder] = sort(v1v2rfOverlapArr); v1v2NoiseCorArr = v1v2NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2rfOverlapArr) & (v1v2rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1v2NoiseCorArr((bin < v1v2rfOverlapArr) & (v1v2rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,1]);
title('Receptive Field and Noise Correlations between v1 and v2'); xlabel('KL Divergence between voxel V1i, V2j RFs (distance)'); ylabel('Noise correlation between voxel V1i, V2j (mm)');


%between v1 and v3
figure(15);hold on;subplot(1,2,1);hold on;

v1v3NoiseCor = transpose(corr(cleanRois(1).vox.baselineNoise,cleanRois(3).vox.baselineNoise));        

for i = 1:min(size(v1v3NoiseCor)); scatter(v1v3rfOverlap(i,:),v1v3NoiseCor(i,:)); end;

v1v3rfOverlapArr = reshape(v1v3rfOverlap,[1 min(size(v1v3rfOverlap))*max(size(v1v3rfOverlap))]);
v1v3NoiseCorArr = reshape(v1v3NoiseCor,[1 min(size(v1v3NoiseCor))*max(size(v1v3NoiseCor))]);
[v1v3rfOverlapArr,sortOrder] = sort(v1v3rfOverlapArr); v1v3NoiseCorArr = v1v3NoiseCorArr(sortOrder);

bins = [];noiseCorAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3rfOverlapArr) & (v1v3rfOverlapArr < bin+step) ) > 15
        noiseCorAvg = median(v1v3NoiseCorArr((bin < v1v3rfOverlapArr) & (v1v3rfOverlapArr < bin+step)));
        bins = [bins bin]; noiseCorAvgs = [noiseCorAvgs noiseCorAvg];
    end
end

plot(bins,noiseCorAvgs,'black','LineWidth',8); xlim([0,1]);
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

for i = 1:length(v1rfOverlap); scatter(v1dist(i,:),v1rfOverlap(i,:)); end;

v1distArr = reshape(v1dist,[1 length(v1dist)^2]); v1rfOverlapArr = reshape(v1rfOverlap,[1 length(v1rfOverlap)^2]);
[v1distArr,sortOrder] = sort(v1distArr); v1rfOverlapArr = v1rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = []; step = .5
for bin = 0:step:40
    if sum( (bin < v1distArr) & (v1distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1rfOverlapArr((bin < v1distArr) & (v1distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V1 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v2 %
figure(25);hold on;subplot(2,3,2);hold on;

for i = 1:length(v2rfOverlap); scatter(v2dist(i,:),v2rfOverlap(i,:)); end;

v2distArr = reshape(v2dist,[1 length(v2dist)^2]); v2rfOverlapArr = reshape(v2rfOverlap,[1 length(v2rfOverlap)^2]);
[v2distArr,sortOrder] = sort(v2distArr); v2rfOverlapArr = v2rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = [];
for bin = 0:step:40
    if sum( (bin < v2distArr) & (v2distArr < bin+step) ) > 15
        rfOverlapAvg = median(v2rfOverlapArr((bin < v2distArr) & (v2distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V2 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


% do v3 %
figure(25);hold on;subplot(2,3,3);hold on;

for i = 1:length(v3rfOverlap); scatter(v3dist(i,:),v3rfOverlap(i,:)); end;

v3distArr = reshape(v3dist,[1 length(v3dist)^2]); v3rfOverlapArr = reshape(v3rfOverlap,[1 length(v3rfOverlap)^2]);
[v3distArr,sortOrder] = sort(v3distArr); v3rfOverlapArr = v3rfOverlapArr(sortOrder);
bins = [];rfOverlapAvgs = []; 
for bin = 0:step:40
    if sum( (bin < v3distArr) & (v3distArr < bin+step) ) > 15
        rfOverlapAvg = median(v3rfOverlapArr((bin < v3distArr) & (v3distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end
plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('V3 Distance and RF Correlations'); xlabel('Distance between voxels i,j (mm)'); ylabel('KL divergence between voxel RFs i,j');


%between v1 and v2
figure(25);hold on;subplot(2,3,4);hold on;

for i = 1:min(size(v1v2rfOverlap)); scatter(v1v2dist(i,:),v1v2rfOverlap(i,:)); end;

v1v2distArr = reshape(v1v2dist,[1 min(size(v1v2dist))*max(size(v1v2dist))]);
v1v2rfOverlapArr = reshape(v1v2rfOverlap,[1 min(size(v1v2rfOverlap))*max(size(v1v2rfOverlap))]);
[v1v2distArr,sortOrder] = sort(v1v2distArr); v1v2rfOverlapArr = v1v2rfOverlapArr(sortOrder);

bins = [];rfOverlapAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v2distArr) & (v1v2distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1v2rfOverlapArr((bin < v1v2distArr) & (v1v2distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end

plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
title('Distance and RF Correlations between v1 and v2'); xlabel('Distance between voxels V1i, V2j (mm)'); ylabel('KL divergence between voxel RFs V1i, V2j (distance)');


%between v1 and v3
figure(25);hold on;subplot(2,3,5);hold on;

for i = 1:min(size(v1v3rfOverlap)); scatter(v1v3dist(i,:),v1v3rfOverlap(i,:)); end;

v1v3distArr = reshape(v1v3dist,[1 min(size(v1v3dist))*max(size(v1v3dist))]);
v1v3rfOverlapArr = reshape(v1v3rfOverlap,[1 min(size(v1v3rfOverlap))*max(size(v1v3rfOverlap))]);
[v1v3distArr,sortOrder] = sort(v1v3distArr); v1v3rfOverlapArr = v1v3rfOverlapArr(sortOrder);

bins = [];rfOverlapAvgs = [];
for bin = 0:step:45
    if sum( (bin < v1v3distArr) & (v1v3distArr < bin+step) ) > 15
        rfOverlapAvg = median(v1v3rfOverlapArr((bin < v1v3distArr) & (v1v3distArr < bin+step)));
        bins = [bins bin]; rfOverlapAvgs = [rfOverlapAvgs rfOverlapAvg];
    end
end

plot(bins,rfOverlapAvgs,'black','LineWidth',8);ylim([0,1]);
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






