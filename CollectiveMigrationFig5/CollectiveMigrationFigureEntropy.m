%% Position Entropy computation in 2D, written initially by B.S. and N.S.
%see also: https://www.mathworks.com/matlabcentral/fileexchange/12857-entropy

%execute from the root directory of cloned repository!

%load in data and set data as "current"
addpath ../ProcessedData 
addpath ../Functions
load('grey3LHVC_dynamic_stats.mat')  %grey3LHVC in manuscript, but can change to different bird/region
current = grey3LHVC_dynamic_stats; %grey3LHVC in manuscript, but can change to different bird/region

clearvars -except current
current = RenumberCells(current); %makes sure all cells going in order with no skips in cell number even if missing/removed cells

%create figure
figure();
subplot(1,2,1)
hold on

%set your microns per pixel scale, number of timepoints and number of bins
microns_per_FOV=current.MaxX(1)*2; timepoints=current.MaxF(1); number_of_bins=200;
% (MaxX*2 is width of FOV, because MaxX is width of one quadrant, which is 1/2 total width)

%% compute entropy for each number of bins and loop over number of timepoints
for t=1:timepoints
 clear XandY
XandY(:,1) = current.X(current.F == t); %X positions for  time point = t
XandY(:,2) = current.Y(current.F == t); %Y positions for  time point = t
numcells(t)=length(XandY);
for nbins=3:number_of_bins
bincounts2d = hist3(XandY,[nbins nbins]); %use hist3 to count the number of cells in each bin
r=hist(bincounts2d(:)); %turn the bin counts into a histogram 
p = r/sum(r); %compute the probability distribution by normalizing p to number of cells
p(p==0)=[]; %remove zeros
bent(nbins,t) = -sum(p.*log(p)); %calculate entropy
scale(nbins)=microns_per_FOV/nbins;
end
end

%plot positions 
plot(XandY(:,1),XandY(:,2),'o','Color',[0.2 0.2 0.2],'LineWidth',2,'MarkerSize',4) 
 ylabel ('Y (microns)')
 xlabel ('X (microns)')
 set(gca,'YDir','reverse')
 title ('Final Cell Position')

 %plot mean entropy[positions for all timepoints] as a function of bin size
subplot(1,2,2)
hold on
 errorbar(log10(scale),mean(bent,2),std(bent')','.k','MarkerSize',10)
 ylabel ('Entropy')
 xlabel ('log10 Bin size (microns)')
 drawnow
 set(gca,'FontSize',12), set(gcf,'color','w')
 
%% compute "maximum entropy" simulation using the multinomial distribution 
for i=1:1000 %choose number of timepoints
for nbins=3:number_of_bins
numbins=nbins^2; %count the total number of bins
probs=ones(numbins,1)/numbins; %compute the probability that a cell gets added to each bin
%now compute "maximum entropy" entropy of a multinomial distribution
simulated_data=mnrnd(ceil(mean(numcells)),probs); %min number cells
r=histcounts(simulated_data);
p=r./sum(r);
ment(nbins,i)= -nansum(p.*log(p)); %calculate entropy
end
end
mment = mean(ment'); %average across timepoints
sment = std(ment');
scale(3:1:number_of_bins) = microns_per_FOV./(3:1:number_of_bins);

% plot the mean entropy of the simulation for each bin size, with standard deviation as confidence interval
plot(log10(scale),mment,'-r','Linewidth',1.5)
plot(log10(scale),sment+mment,'--','Linewidth',1,'Color',[1, 0, 0, 0.4]) %upper std bound
plot(log10(scale),mment - sment,'--','Linewidth',1,'Color',[1, 0, 0, 0.4]) %lower std bound
x2 = [log10(scale), fliplr(log10(scale))];
inBetween = [mment - sment, fliplr(sment+mment)];
x2(~isfinite(x2)) = x2(3);
fill(x2(inBetween>0), inBetween(inBetween>0), 'r','FaceAlpha',0.1,'EdgeColor','none');
ylabel ('Entropy')
xlabel ('log10 Bin size (microns)')
set(gca,'FontSize',12), set(gcf,'color','w')


%print percent error between data and entropy curve
PE = (trapz(mean(bent,2)) - trapz(mment))/trapz(mean(bent,2));
disp(['Percent error: ' num2str((PE*100)) '%'])


%% Heading Entropy 
clearvars -except current
%=========USER INPUT (keep default for replicating paper figure)
timebin = max(current.F); %setting for the figure in Shvedov et al.: timebin = max(current.F)
arrowwidth= 0.8; %'0.8' in paper
colorarrows = 0; %if 1, will color according to colorvec. '0' in paper
    colorvec = current.Z; %set the vector you want to color arrows by. 'current.Z' or 'current.F' recommended

microns_per_FOV=current.MaxX(1)*2; %multiplies the maxX by 2 (a quadrant diameter *2 is the diameter of FOV)
current = RenumberCells(current); %makes sure all cells going in order with no skips
 
%find heading difference across time by finding heading between current and first position over time
figure();
subplot(1,2,1)
hold on

%set up variables
vector = [];
numcells=max(current.N);

if colorarrows == 1
%for coloring the vectors
c = parula(ceil(max(colorvec) - min(colorvec)));
else
end

%plot vectors
for i=1:numcells
    thiscellXdata=current.X(current.N == i);
    thiscellYdata=current.Y(current.N == i);
    vector(i,:)=[thiscellXdata(end)-thiscellXdata(1) thiscellYdata(end)-thiscellYdata(1)]; 
    if colorarrows == 1
    quiver(vector(i,1),vector(i,2),'Color',c(ceil(mean(colorvec(current.N == i)) - min(colorvec)),:), 'AutoScale', 'off', 'LineWidth', arrowwidth)
    elseif colorarrows == 0
        quiver(vector(i,1),vector(i,2),'Color',[0.1 0.1 0.1],'LineWidth',arrowwidth)
    end
    drawnow
end
title('Cell Heading')
ylabel ('Y (microns)')
xlabel ('X (microns)')
set(gca,'FontSize',12), set(gcf,'color','w')
set(gca,'YDir','reverse')
if colorarrows == 1
    colorbar
    caxis([min(colorvec) max(colorvec)])
else
end


%compute cumulative heading entropy over time, set up vectors
vectors = []; %will be 3D matrix with 1st dim: cell, 2nd dim: X and Y vector, 3rd dim: increasing time since start
cells = unique(current.N);
tick = 0;
for t = 2+1:timebin %all possible frames in the data, starting at 2 because many don't have 1 in grey3L
    count = 0; %will count # of cells that are not available for this timepoint
    tick = tick+1; %counting # of timepoints
for i=1:length(cells); %loop through cells
    thiscell=current(current.N == cells(i),:); %isolate cell
    if ~ismember(t,unique(thiscell.F)) % if cell doesn't have this timepoint
        vectors(tick,:,i) = [NaN NaN NaN]; %just fill in nans
        count = count+1; %count this occurence to be subtracted from cell counts later 
        continue
    else
    vectors(tick,:,i)=[thiscell.X(thiscell.F == t)-thiscell.X(1) thiscell.Y(thiscell.F == t)-thiscell.Y(1) t]; %find increasing vector over time, timepoint that subtracting first cell position from
    end
end
numcells(tick) = length(cells) - count; %number of cells used for comparisons
end

%compute entropy
for t = 1:size(vectors,1)
for nbins=2:100
    numbins=nbins^2;
    bincounts2d = hist3(squeeze(vectors(t,1:2,:))',[nbins nbins]); %use hist3 to count the number of cells in each bin
    r=histcounts(bincounts2d(:));
    p = r/sum(r);%normalize p to number of cells
    p(p==0)=[];
    bent(nbins,t) = -sum(p.*log(p)); %calculate entropy
end
end
scale(2:1:100)=microns_per_FOV./(2:1:100);

subplot(1,2,2)
hold on
errorbar(log10(scale),mean(bent'),std(bent'),'.k','MarkerSize',10)
ylabel ('Entropy')
xlabel ('log10 Bin size')
drawnow
set(gca,'FontSize',12), set(gcf,'color','w')

%% SIMULATION
%choose number of bins and loop
for i=1:1000 %choose number of timepoints, i=1:1000 for paper
for nbins=2:100
numbins=nbins^2; %count the total number of bins
probs=ones(numbins,1)/numbins; %compute the probability that a cell gets added to each bin
%now compute "maximum entropy" entropy of a multinomial distribution
simulated_data=mnrnd(ceil(mean(numcells)),probs); %min numb cells
r=histcounts(simulated_data);
p=r./sum(r);
ment(nbins,i)= -nansum(p.*log(p)); %calculate entropy
end
end
mment = mean(ment'); %average across timepoints
sment = std(ment');
scale(2:1:100) = microns_per_FOV./(2:1:100);

% plot the mean entropy of the simulation for each bin size, with standard deviation as confidence interval
plot(log10(scale),mment,'-r','Linewidth',1.5)
plot(log10(scale),sment+mment,'--','Linewidth',1,'Color',[1, 0, 0, 0.4]) %upper std bound
plot(log10(scale),mment - sment,'--','Linewidth',1,'Color',[1, 0, 0, 0.4]) %lower std bound
x2 = [log10(scale), fliplr(log10(scale))];
inBetween = [mment - sment, fliplr(sment+mment)];
x2(~isfinite(x2)) = x2(3);
fill(x2(inBetween>0), inBetween(inBetween>0), 'r','FaceAlpha',0.1,'EdgeColor','none');
ylabel ('Entropy')
xlabel ('log10 Bin size (microns)')
set(gca,'FontSize',12), set(gcf,'color','w')

%print percent error between data and entropy curve
PE = (trapz(mean(bent,2)) - trapz(mment))/trapz(mean(bent,2));
disp(['Percent error: ' num2str((PE*100)) '%'])