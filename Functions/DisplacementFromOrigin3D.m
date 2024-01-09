function DisplacementFromOrigin3D(current)

%run workflow_for_compiled_tracks.m first

%set current = dynamic_stats;

separatetorts = 0;
    tortthresh = 1.7;
    uppertortthresh = 5;
    
plotdots = 1;
    



if separatetorts == 0
cells = unique(current.N);
figure()
title('XYZ Displacement Vectors of all Cells')
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

elseif separatetorts == 1
    radialN = [];
    wandN = [];
    cells = unique(current.N);
    r = 0;
    w = 0;
      for c = 1:length(cells)
          nID = cells(c);
           track = [];
           track = current(find(current.N(:) == nID),:);
                if (sum(abs(track.pathdiffs(2:end)))/(track.currentdisp(end))) < tortthresh %if tortuosity is < tortthresh (radial)
                    r = r+1;
                    radialN(r) = track.N(1);
                elseif (sum(abs(track.pathdiffs(2:end)))/(track.currentdisp(end))) > tortthresh & (sum(abs(track.pathdiffs(2:end)))/(track.currentdisp(end))) < uppertortthresh
                    w = w+1;
                    wandN(w) = track.N(1);
                end
      end
%figure()
title(['XYZ Displacement Vectors of Radial Tracks, T < ' num2str(tortthresh) ', n = ' num2str(length(radialN))])
hold on
for n = 1:length(radialN)
        coor = find(current.N(:) == radialN(n));
       if length(coor) < 3
            continue
        else
        patch([current.currentdispX(coor)' nan], [current.currentdispY(coor)' nan], [current.currentdispZ(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2,'EdgeAlpha',0.5)
       if plotdots == 1
        scatter3(current.currentdispX(coor(end)), current.currentdispY(coor(end)),current.currentdispZ(coor(end)), 'filled','red')
       else
       end
    end
end
xlabel('X displacement (um)')
ylabel('Y displacement (um)')
zlabel('Z displacement (um)')
xlim([-60 60])
ylim([-60 60])
zlim([-60 60])
hold off
set(gcf,'color','w');
view(3) 
grid on

figure()
title(['XYZ Displacement Vectors of Wandering Tracks, T > ' num2str(tortthresh) ', n = ' num2str(length(wandN))])
hold on
for n = 1:length(wandN)
        coor = find(current.N(:) == wandN(n));
        if length(coor) < 3
            continue
        else
        patch([current.currentdispX(coor)' nan], [current.currentdispY(coor)' nan], [current.currentdispZ(coor)' nan], [current.F(coor)' nan], 'FaceColor', 'none','EdgeColor','interp','LineWidth',2, 'EdgeAlpha',0.5)
       if plotdots == 1
        scatter3(current.currentdispX(coor(end)), current.currentdispY(coor(end)),current.currentdispZ(coor(end)), 'filled','red')
       else
       end
       end
end
xlabel('X displacement (um)')
ylabel('Y displacement (um)')
zlabel('Z displacement (um)')
xlim([-60 60])
ylim([-60 60])
zlim([-60 60])
hold off
set(gcf,'color','w');
view(3)
grid on
end

colormap(jet)
set(gca, 'Zdir', 'reverse')
colorbar
%set(gca, 'Ydir', 'reverse')

end