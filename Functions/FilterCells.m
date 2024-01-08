% filter cells based on total displacement and tortuosity. 

function [current] = FilterCells(current, dispthresh, tortthresh)
%inputs are data, displacement threshold (5) and tortuosity threshold (7)

current = RenumberCells(current);
cells = unique(current.N);
remove = [];
rc = 0; %removed cells counter
    for ii = 1:length(cells)
        nID = cells(ii);
        track = current(current.N == nID,:);
        clear pathX pathY pathZ paths_3d_cell hourdiffs currentdispXYZ currentdisp
        tc = 0;
        for t = 2:length(track.F);
           tc = tc+1;
           pathX(tc) = (track.X(t) - track.X(t-1));
           pathY(tc) = (track.Y(t) - track.Y(t-1));
            pathZ(tc) = (track.Z(t) - track.Z(t-1));
            paths_3d_cell(tc) = sqrt(pathX(tc)^2 + pathY(tc)^2 + pathZ(tc)^2);
            currentdispXYZ(tc,:) = [(track.X(t) - track.X(1)), (track.Y(t) - track.Y(1)),(track.Z(t) - track.Z(1))];
            currentdisp(tc) = sqrt(sum(currentdispXYZ(tc,1)^2 + currentdispXYZ(tc,2)^2 + currentdispXYZ(tc,3)^2));
           hourdiffs(tc) = ((track.H(t))-(track.H(t-1))); %hourly difference
        end
        if (sum(paths_3d_cell)/currentdisp(end)) > tortthresh || currentdisp(end) < dispthresh
            rc = rc+1;
            remove(rc) = nID;
        end
    end

current(ismember(current.N, remove),:) = []; %removes all filtered cells
end