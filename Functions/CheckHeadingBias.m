function CheckHeadingBias(current,nbin, dots, colors) %current is neuroblast struct, nbin is value, colors is string

%pathh =
%'/Users/naomishvedov/Desktop/Scott_Lab/Useful_Scripts/A_FINAL_PIPELINE/CircStat2012a'; %laptop
pathh = '/Users/naomishvedov/Documents/GitHub/In-vivo-Migratory-Dynamics/CircStat2012a'; %for lab PC
%nbin = number of bins
%current = make current neuroblaststats
%colors = color string 
%dots = 1, plot endpoint dot. dots = 0, just the polar histogram
%do Rayleigh test and polarhistogram. calculate angle relative to reference
%frame which is center of the volume

if ~exist('colors','var')
    colors = 'auto';
end

%FUNCTION ASSUMES:
%ALL QUADRANTS ARE THE SAME SIZE 
%FINDING CENTER OF TRACKED DATA, not WHOLE VOLUME
cells = unique([current.N]);
azims = 0;
elevs = 0;
rho = 0;
n = 0;

for n = 1:length(current);
    if isempty(current(n).dispXYangle) == 1
        continue
    else
    azims(n) = current(n).dispXYangle;
    elevs(n) = current(n).dispXZangle;
    rhos(n) = current(n).disptot/10;
    end
end

addpath(pathh)
rA = circ_rtest(azims); %Rayleigh's test from Circular Statistics Toolbox 
rE = circ_rtest(elevs);
wA = circ_symtest(azims);
wE = circ_symtest(elevs);

figure()
polarhistogram(azims,nbin,'FaceColor',colors) 
if dots == 1
    hold on
    polarscatter(azims, rhos, 'filled','red')
else
end
set(gcf,'color','w'), set(gca,'FontSize',16);
set(gca,'ThetaZeroLocation','bottom')
title('Neuroblast XY Angles')
subtitle(['Rayleighs p value: ' num2str(rA) ' , Wilcoxon signed-rank p value: ' num2str(wA)])

figure()
polarhistogram(elevs,nbin,'FaceColor',colors) 
if dots == 1
    hold on
    polarscatter(azims, rhos, 'filled', 'red')
else
end
set(gcf,'color','w'), set(gca,'FontSize',16);
set(gca,'ThetaZeroLocation','bottom')
title('Neuroblast XZ Angles')
subtitle(['Rayleighs p value: ' num2str(rE) ' , Wilcoxon signed-rank p value: ' num2str(wE)])

 figure()
 subplot(1,2,1)
 title('XY Displacement Vectors of all Cells')
 hold on
for n = 1:length(current)
     plot([0 current(n).dispX],[0 current(n).dispY],'LineWidth',1)
end
xlabel('X displacement (um)')
ylabel('Y displacement (um)')
hold off
 subplot(1,2,2)
 title('XZ Displacement Vectors of all Cells')
 hold on
for n = 1:length(current)
     plot([0 current(n).dispY],[0 current(n).dispZ],'LineWidth',1)
end
xlabel('X displacement (um)')
ylabel('Z displacement (um)')
hold off
set(gcf,'color','w');

end
