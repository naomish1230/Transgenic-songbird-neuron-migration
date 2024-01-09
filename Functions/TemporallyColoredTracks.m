function [current] = TemporallyColoredTracks(current);     
%input is data, output is figure and quadrant adjusted data
%code will only run if not ran before this session

if length(unique(current.Q)) == 1; %if data only has 1 quadrant tracked, use this to plot: 
cells = unique(current.N);
figure()
colormap(jet)
    hold on
    for n = cells(1):cells(end); %patch allows temporal color coding
        coor = find(current.N(:) == n);
        if length(current.F(coor)) < 3
        continue
        else
        patch([current.X(coor)' nan], [current.Y(coor)' nan], [current.Z(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2)
        end
    end
    h = colorbar;
    ylabel(h, 'Timepoint')
    title(['Neuroblasts 1 - ' num2str(length(cells))])
    xlabel('X Position (um)')
    xlim([1 current.MaxX(1)*2]) 
    ylabel('Y Position (um)')
    ylim([1 current.MaxY(1)*2])
    zlabel('Z Position (um)')
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
    zlim([1 (current.MaxZ(1))*0.5])
    view(3)
    set(gcf,'color','w');
    grid
    hold off
elseif length(unique(current.Q)) >1 & ismember('fixedQ',current.Properties.VariableNames) == 0; %if there are multiple quadrants and has not already been adjusted
    current.fixedQ(:) = 1; %column indicating quadrant adjustment, skipped if already adjusted
    if find(current.Q(:) == 2) > 0 %for upper right quadrant, add to X but not Y
    current.X((find(current.Q(:) == 2))) = current.X((find(current.Q(:) == 2))) + current.MaxX((find(current.Q(:) == 2)));
    end
    if find(current.Q(:) == 3) > 0 %for lower left quadrant, add to Y but not X
    current.Y((find(current.Q(:) == 3))) = current.Y((find(current.Q(:) == 3))) + current.MaxY((find(current.Q(:) == 3)));
    end
    if find(current.Q(:) == 4) > 0 %for lower right quadrant, add to X and Y
    current.X((find(current.Q(:) == 4))) = current.X((find(current.Q(:) == 4))) + current.MaxX((find(current.Q(:) == 4)));  
    current.Y((find(current.Q(:) == 4))) = current.Y((find(current.Q(:) == 4))) + current.MaxY((find(current.Q(:) == 4)));
    end
    cells = unique(current.N);
figure()
colormap(jet)
    hold on
    for n = cells(1):cells(end); %patch allows temporal color coding
        coor = find(current.N(:) == n);
        if length(current.F(coor)) < 3
            continue
        else
        patch([current.X(coor)' nan], [current.Y(coor)' nan], [current.Z(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2)
        end
   end
    h = colorbar;
    ylabel(h, 'Timepoint')
    title(['Neuroblasts 1 - ' num2str(length(cells))])
    xlabel('X Position (um)')
    xlim([1 current.MaxX(1)*2]) 
    ylabel('Y Position (um)')
    ylim([1 current.MaxY(1)*2])
    zlabel('Z Position (um)')
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
    zlim([1 (current.MaxZ(1))*0.5])
    view(3)
    set(gcf,'color','w');
    grid
    hold off
elseif ismember('fixedQ',current.Properties.VariableNames) == 1;
      cells = unique(current.N);
figure() 
colormap(jet)
    hold on
    for n = cells(1):cells(end); %patch allows temporal color coding
        coor = find(current.N(:) == n);
        if length(current.F(coor)) < 3
            continue
        else
        patch([current.X(coor)' nan], [current.Y(coor)' nan], [current.Z(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2)
        end
   end
    h = colorbar;
    ylabel(h, 'Timepoint')
    title(['Neuroblasts 1 - ' num2str(length(cells))])
    xlabel('X Position (um)')
    xlim([1 current.MaxX(1)*2]) 
    ylabel('Y Position (um)')
    ylim([1 current.MaxY(1)*2])
    zlabel('Z Position (um)')
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
    zlim([1 (current.MaxZ(1))*0.5])
    view(3)
    set(gcf,'color','w');
    grid
    hold off
end

end 
