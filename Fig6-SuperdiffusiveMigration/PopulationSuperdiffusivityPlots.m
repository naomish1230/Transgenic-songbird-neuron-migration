%PopulationSuperdiffusivity
% this script creates all figures for superdiffusivity (figure 6).

%execute from the root directory of cloned repository!

clear 
addpath ../ProcessedData 
addpath ../Functions

%load in ALL birds' combined data
load('combined_ALL_dynamic_stats.mat');
current = combined_ALL_dynamic_stats;

% Filter outliers:
current = FilterCells(current,5,7); %input: data, displacement threshold (5), tortuosity threshold (7)
%% PLOT 1: Male HVC, including grey3L, green4R, whiteband, and pink61 HVC. 
%lime12 not shown because of hourly temporal acquisition

clear mdl61 mdl3 mdl0 mdl4 %clear variables if not re-running only the section
hvc = figure(); %make figure
hold on

%for pink61 HVC: (dark orange)
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(current(current.Bird_ID == 61 & current.R == 0,:));
mdl61.mdl = mdl; %save output of the model
mdl61.msd = msd;
mdl61.times = times;
mdl61.rsq = rsq;
scatter([0 mdl61.times],[0 mdl61.msd],100,[0.8500 0.3250 0.0980],'filled') %,'LineWidth',2) %,'Color','black') %size is 100
plot([0 mdl61.times],[0 (mdl61.mdl.a*mdl61.times .^mdl61.mdl.b)],'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)

%for grey3L HVC: (light orange)
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(current(current.Bird_ID == 3 & current.R == 0,:));
mdl3.mdl = mdl; %save output of the model
mdl3.msd = msd;
mdl3.times = times;
mdl3.rsq = rsq;
scatter([0 mdl3.times],[0 mdl3.msd],100,[0.9290 0.6940 0.1250],'filled') %,'LineWidth',2) %,'Color','black') %size is 100
plot([0 mdl3.times],[0 (mdl3.mdl.a*mdl3.times .^mdl3.mdl.b)],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)


%for green4R HVC: (purple)
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(current(current.Bird_ID == 4,:));
mdl4.mdl = mdl; %save output of the model
mdl4.msd = msd;
mdl4.times = times;
mdl4.rsq = rsq;
scatter([0 mdl4.times],[0 mdl4.msd],100,[0.4940 0.1840 0.5560],'filled') %,'LineWidth',2) %,'Color','black') %size is 100
plot([0 mdl4.times],[0 (mdl4.mdl.a*mdl4.times .^mdl4.mdl.b)],'--','Color',[0.4940 0.1840 0.5560],'LineWidth',1.5)


%for whiteband S1 HVC: (blue)
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(current(current.Bird_ID == 0 & current.S == 1,:));
mdl0.mdl = mdl; %save output of the model
mdl0.msd = msd;
mdl0.times = times;
mdl0.rsq = rsq;
scatter([0 mdl0.times],[0 mdl0.msd],100,[0 0.4470 0.7410],'filled') %,'LineWidth',2) %,'Color','black') %size is 100
plot([0 mdl0.times],[0 (mdl0.mdl.a*mdl0.times .^mdl0.mdl.b)],'--','Color',[0 0.4470 0.7410],'LineWidth',1.5)

%aesthetics
set(gca,'FontSize',16), set(gcf,'color','w')
ylabel('MSD (um^2)')
xlabel('Time lag (minutes)')
title('Single Bird HVC fits')
hold off

%% across bird groups' fits

grps = figure();
hold on


%male hyperpallium (HP) - grey3L and pink61 HP plotted in blue
maleHP = current(current.Bird_ID == 3 | current.Bird_ID == 61,:);
maleHP = maleHP(maleHP.R == 1,:);
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(maleHP);
mdlhp.mdl = mdl;
mdlhp.msd = msd;
mdlhp.times = times;
mdlhp.rsq = rsq;
scatter([0 mdlhp.times],[0 mdlhp.msd],100,[0.3010 0.7450 0.9330],'filled')
plot([0 mdlhp.times],[0 (mdlhp.mdl.a*mdlhp.times .^mdlhp.mdl.b)], '--','Color',[0.3010 0.7450 0.9330],'LineWidth',1.5)


%female NP (blue53, black11L)
femNP = [current(current.Bird_ID == 53 & current.S == 1,:) ; current(current.Bird_ID == 11,:)] ; %current(current.Bird_ID == 5 & current.S ~= 3,:)];
%femNP = current(current.Bird_ID == 11,:);
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(femNP);
mdlnp.mdl = mdl;
mdlnp.msd = msd;
mdlnp.times = times;
mdlnp.rsq = rsq;
scatter([0 mdlnp.times],[0 mdlnp.msd],100,[0.4660 0.6740 0.1880],'filled')
plot([0 mdlnp.times],[0 (mdlnp.mdl.a*mdlnp.times .^mdlnp.mdl.b)], '--','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)

%male HVC (green4R, grey3L, pink61)... plotted in black
maleHVC = current(current.Bird_ID == 4 | current.Bird_ID == 3 | current.Bird_ID == 61,:);
maleHVC = maleHVC(maleHVC.R == 0,:);
[mdl,msd,times,alpha,rsq] = ComputeMSD_function(maleHVC);
mdlhvc.mdl = mdl; %save output of the model
mdlhvc.msd = msd;
mdlhvc.times = times;
mdlhvc.rsq = rsq;
scatter([0 mdlhvc.times],[0 mdlhvc.msd],100,'filled','black') %size is 100
plot([0 mdlhvc.times],[0 (mdlhvc.mdl.a*mdlhvc.times .^mdlhvc.mdl.b)],'--','Color',"black",'LineWidth',1.5)

%aesthetics
set(gca,'FontSize',16), set(gcf,'color','w')
ylabel('MSD (um^2)')
xlabel('Time lag (minutes)')
title('Combined Bird Fits')
hold off
%% Individual cell fit

ind = figure;
hold on

m = 115; %chosen cell ID from grey3L
indcell = RenumberCells(current(current.Bird_ID == 3,:));
[indalpha, indrsq, indbeta, tracks] = AlphaVsDepth(indcell,2);
times = 10*tracks{m}.frames;
[mdl,gof] = fit(times(1:length(tracks{m}.msd)),tracks{m}.msd','power1');
scatter([0 times(1:length(tracks{m}.msd))'],[0 tracks{m}.msd],100,[0.6350 0.0780 0.1840],'filled') %size is 100
plot([0 times(1:length(tracks{m}.msd))'],[0 mdl.a*times(1:length(tracks{m}.msd)).^mdl.b'],'--','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5)

%aesthetics
set(gca,'FontSize',16), set(gcf,'color','w')
ylabel('MSD (um^2)')
xlabel('Time lag (minutes)')
title(['Cell ' num2str(m) ' Fit, alpha = ' num2str(mdl.b) ', R2 = ' num2str(gof.rsquare)])
hold off

%% fit to all cells 

[allalpha, allrsq, allbeta, tracks, alldepths, alltorts] = AlphaVsDepth(current,2);


    h = histfit(allalpha(allalpha < 4 & allrsq>0.95),round(max(allalpha(allalpha < 4 & allrsq>0.95))/0.2),'Loglogistic') %filtering coeffs by being positive and having good Rsquared value
    ndl = fitdist(allalpha(allalpha < 4 & allrsq>0.95)','Loglogistic')
    [testdec,p,stats] = chi2gof(allalpha(allalpha < 4 & allrsq>0.95)','CDF',ndl)


legend('All cells','Log-logistic Distribution Fit')
set(gcf,'color','w'),set(gca,'FontSize',16)
xlabel('Alpha coefficient')
ylabel('Cell count')
title(['Mean alpha (R2>0.95): ' num2str(mean(allalpha(allalpha < 4 & allrsq>0.95)))]);