% Collective Migration Figure: Correlations
% using 5 high temporal resolution birds. 

%==========  USER INPUT 
%spec_birds = [5,11,3,4,61]; %all birds!
spec_birds = [4,11,3,61]; %something weird happens when add in lime5L because of some NaNs
%============

clearvars -except spec_birds
addpath('/Users/naomishvedov/Documents/GitHub/In-vivo-Migratory-Dynamics/CollectiveMigrationFig/')
load('birds_headCorrs_table.mat')
load('birds_speedCorr_table.mat')

birds_headCorrs_table = birds_headCorrs_table(ismember(birds_headCorrs_table.ID,spec_birds),:);
birds_speedCorr_table = birds_speedCorr_table(ismember(birds_speedCorr_table.ID, spec_birds),:);

%% Heading Correlations 

%remove correlations that are with self
birds_headCorrs_table(birds_headCorrs_table.DistDiffs < 5,:) = [];

%are correlations in the population above chance?
figure()
h = histogram(birds_headCorrs_table.AngleDiff,36,'FaceColor',"black",'EdgeColor','none','FaceAlpha',0.4); %DATA
hold on
%histogram(birds_headCorrs_table.ShuffAngleDiff,'FaceColor',"yellow",'FaceAlpha',0.55)
[counts,centers] = hist(birds_headCorrs_table.ShuffAngleDiff,h.NumBins); %shuffled data
plot(centers,counts,'LineWidth',3,'Color','blue') %SHUFFLED DATA LINE plotted
set(gca,'FontSize',16), set(gcf,'color','w');
legend('Data','Shuffled')
[h,p] = kstest2(birds_headCorrs_table.AngleDiff, birds_headCorrs_table.ShuffAngleDiff); %K-S test, 2 sample
title(['Angle Differences, n = ' num2str(height(birds_headCorrs_table)) ' pairs, p = ' num2str(p)])
ylabel('Pair Count')
xlabel('Degrees')

% is there a relationship between correlation and distance?
angmdl = fitlm(reshape(birds_headCorrs_table.DistDiffs,1,[]), reshape(birds_headCorrs_table.AngleDiff,1,[]))
angmdl_shuff = fitlm(reshape(birds_headCorrs_table.ShuffDistDiffs,1,[]), reshape(birds_headCorrs_table.ShuffAngleDiff,1,[]))


%plot the scatter plot as density heatmap scatter plot
figure()
hold on
C = ksdensity([birds_headCorrs_table.DistDiffs(:),birds_headCorrs_table.AngleDiff(:)], [birds_headCorrs_table.DistDiffs(:),birds_headCorrs_table.AngleDiff(:)]); %compute density of each point
h = scatter(birds_headCorrs_table.DistDiffs(:),birds_headCorrs_table.AngleDiff(:), [], C,'filled');
xlabel('Distance Difference (um)'),ylabel('Degrees')
colorbar
title('Angle Difference vs. Distance'),set(gca,'FontSize',16), set(gcf,'color','w')
%refline(angmdl.Coefficients.Estimate(2),angmdl.Coefficients.Estimate(1)) %reference line with model fit
hold off

%% Speed Correlations

%remove correlations that are with self
birds_speedCorr_table(birds_speedCorr_table.DistDiffs < 5,:) = [];

figure()
h = histogram(birds_speedCorr_table.SpeedR,'FaceColor',"black",'EdgeColor','none','FaceAlpha',0.4); %DATA
hold on
%histogram(birds_speedCorr_table.ShuffSpeedR)
[counts,centers] = hist(birds_speedCorr_table.ShuffSpeedR,h.NumBins); %shuffled data
plot(centers,counts,'LineWidth',3,'Color','blue') %SHUFFLED DATA LINE plotted
set(gca,'FontSize',16), set(gcf,'color','w');
legend('Data','Shuffled')
[h,p] = kstest2(birds_speedCorr_table.SpeedR, birds_speedCorr_table.ShuffSpeedR); %K-S test, 2 sample
title(['Speed Correlation R, n = ' num2str(height(birds_speedCorr_table)) ' pairs, p =' num2str(p)])
ylabel('Pair Count')
xlabel('Speed Correlation')

% is there a relationship between correlation and distance?
spmdl = fitlm(reshape(birds_speedCorr_table.DistDiffs,1,[]), reshape(birds_speedCorr_table.SpeedR,1,[]))
spmdl_shuff = fitlm(reshape(birds_speedCorr_table.ShuffDistDiffs,1,[]), reshape(birds_speedCorr_table.ShuffSpeedR,1,[]))

%plot the scatter plot as density heatmap scatter plot
figure()
hold on
C = ksdensity([birds_speedCorr_table.DistDiffs(:),birds_speedCorr_table.SpeedR(:)], [birds_speedCorr_table.DistDiffs(:),birds_speedCorr_table.SpeedR(:)]); %compute density of each point
h = scatter(birds_speedCorr_table.DistDiffs(:),birds_speedCorr_table.SpeedR(:), [], C,'filled');
xlabel('Distance Difference (um)'),ylabel('Speed Correlation')
colorbar
title('Speed Correlation vs. Distance'),set(gca,'FontSize',16), set(gcf,'color','w')
hline = refline(spmdl.Coefficients.Estimate(2),spmdl.Coefficients.Estimate(1)); %reference line with model fit
hline.Color = 'red';
hold off
