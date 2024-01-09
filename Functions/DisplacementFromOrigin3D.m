function [f] = DisplacementFromOrigin3D(current, plotdots)

%set current = dynamic_stats;
    
cells = unique(current.N);
f = figure();
title('Relative XYZ Displacement of all Cells')
hold on
for c = 1:length(cells)
        nID = cells(c);
        coor = find(current.N(:) == nID);
        if length(coor) < 3
            continue
        else
        patch([current.currentdispX(coor)' nan], [current.currentdispY(coor)' nan], [current.currentdispZ(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2)
        if plotdots == 1
            scatter3(current.currentdispX(coor(end)), current.currentdispY(coor(end)),current.currentdispZ(coor(end)), 'filled','red')
        else
        end
        end
end
xlabel('X displacement (um)')
xlim([-60 60])
ylabel('Y displacement (um)')
ylim([-60 60])
zlabel('Z displacement (um)')
zlim([-60 60])
hold off
set(gcf,'color','w');
view(3)
grid on
colormap(jet)
set(gca, 'Zdir', 'reverse')
colorbar

end