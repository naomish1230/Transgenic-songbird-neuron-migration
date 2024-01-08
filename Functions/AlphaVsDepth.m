%compute MSD coefficients of all cells and their depths

function [alphas, rsqs, betas, tracks, depths, torts, speeds, disps, ventralmovement] = AlphaVsDepth(current,dimensions) 

if  nargin < 2
    dimensions = 2; %default is 2 dimensions
end

%% set framerate per bird
current = RenumberCells(current); %just in case out of order or deleted cells
current.framerate = nan(height(current),1);
    for h = 1:height(current)
        if current.Bird_ID(h) == 61 || current.Bird_ID(h) == 3 || current.Bird_ID(h) == 4; %pink61, grey3L, green4Rprism2
             current.framerate(h) = 0.1;
        elseif current.Bird_ID(h) == 12; %lime12
            current.framerate(h) = 0.01666667;
        elseif current.Bird_ID(h) == 0 & (current.S(h) == 1 | current.S(h) == 0); %whiteband S1
            current.framerate(h) = 0.047962; % once every 20 minutes
        elseif current.Bird_ID(h) == 0 & current.S(h) == 2; %whiteband S2
            current.framerate(h) =  0.04902; % once every 20 minutes
        elseif current.Bird_ID(h) == 53; %blue53 S1
            current.framerate(h) = 0.02; %once every 48 min
        elseif current.Bird_ID(h) == 5 & current.S(h) == 3; %lime5Lprism4 high temp res
            current.framerate(h) = 0.22124; %1 frame every 4 minutes and 52 seconds.
        elseif current.Bird_ID(h) == 5 & current.S(h) ~= 3; %lime5Lprism4 normal
            current.framerate(h) = 0.083333; %once every 12 minutes
        elseif current.Bird_ID(h) == 11; %1 frame every 12 minutes
            current.framerate(h) = 0.083333;
        end
    end

    %% compute alpha of all cells
longestTrack = 0;
cellsToSkip = [];
cells = unique(current.N);
numberOfCells = length(cells); %robust to length ~= cellIDs (if missing cells)
tracks = cell(1,numberOfCells);
framerates = [];
    for ii = 1:numberOfCells
        nID = cells(ii);
        track = current(current.N == nID,:);
        clear pathX pathY pathZ paths_3d_cell hourdiffs currentdispXYZ currentdisp
        tc = 0;
        for t = 2:length(track.F);
           tc = tc+1;
           pathX(tc) = (track.X(t) - track.X(t-1));
           pathY(tc) = (track.Y(t) - track.Y(t-1));
           if dimensions == 3
            pathZ(tc) = (track.Z(t) - track.Z(t-1));
            paths_3d_cell(tc) = sqrt(pathX(tc)^2 + pathY(tc)^2 + pathZ(tc)^2);
            currentdispXYZ(tc,:) = [(track.X(t) - track.X(1)), (track.Y(t) - track.Y(1)),(track.Z(t) - track.Z(1))];
            currentdisp(tc) = sqrt(sum(currentdispXYZ(tc,1)^2 + currentdispXYZ(tc,2)^2 + currentdispXYZ(tc,3)^2));
           elseif dimensions == 2
               paths_3d_cell(tc) = sqrt(pathX(tc)^2 + pathY(tc)^2);
               currentdispXYZ(tc,:) = [(track.X(t) - track.X(1)), (track.Y(t) - track.Y(1))];
               currentdisp(tc) = sqrt(sum(currentdispXYZ(tc,1)^2 + currentdispXYZ(tc,2)^2));
           end
           hourdiffs(tc) = ((track.H(t))-(track.H(t-1))); %hourly difference
       end
        numberOfSamples = size(track,1);
        tracks{ii} = struct();
        numberOfMsds = round(numberOfSamples-1);
        tracks{ii}.msd = zeros(1, numberOfMsds);
        tracks{ii}.msdUncertainty = zeros(1, numberOfMsds);
        tracks{ii}.cellID = nID; %exact cell ID 
        tracks{ii}.framerate = current.framerate(current.N == nID);
        tracks{ii}.frames = track.F;
        tracks{ii}.torts = sum(paths_3d_cell)/currentdisp(end); %tortuosity
        tracks{ii}.speed = sum(paths_3d_cell)/(max([track.H]) - min([track.H]));
        tracks{ii}.stepsize = paths_3d_cell(:);
        tracks{ii}.currentdisp = currentdisp(:);
        if ~ismember('DVDist',current.Properties.VariableNames)
            tracks{ii}.depth = track.Z;
        elseif unique(track.Bird_ID) == 5 |unique(track.Bird_ID) == 11 |unique(track.Bird_ID) == 4 %the prism birds
            tracks{ii}.depth = track.DVDist;
            tracks{ii}.ventralmovement = track.DVDist(end) - track.DVDist(1);
        end
        
        framerates(ii) = tracks{ii}.framerate(1); %index framerates to find minimum later

        if numberOfMsds==0
            cellsToSkip = [cellsToSkip ii];
        else
            track = table2array(track(:,3:5)); %turns track into array following code can manipulate
            for tau = 1:numberOfMsds
                dx = track(1:1:(end-tau), 1) - track((tau+1):1:end, 1);
                dy = track(1:1:(end-tau), 2) - track((tau+1):1:end, 2);
                if dimensions == 3
                   dz = track(1:1:(end-tau), 3) - track((tau+1):1:end, 3);
                   drSquared = dx.^2+dy.^2 + dz.^2;
                elseif dimensions == 2
                   drSquared = dx.^2+dy.^2;
                end
                msd = mean( drSquared );

                tracks{ii}.msd(tau) = msd;
                tracks{ii}.msdUncertainty(tau) = std( drSquared ) / sqrt(length( drSquared ));
            end
if isnan(tracks{ii}.framerate(1))
    error('Need to provide approximate frame rate for this Bird_ID manually.')
end
            tracks{ii}.diffusionCoefficient = tracks{ii}.msd(1) / (2 * dimensions * tracks{ii}.framerate(1)^-1); %2*dimension*time
            
            longestTrack = max( longestTrack, numberOfSamples );  
        end
    end

times = (1:longestTrack)/min(framerates);
clear torts speeds alphas betas depths rsqs
    for m = 1:length(tracks)
        maxlength = length(tracks{m}.msd);
        
        if maxlength > 1
            [mdl,gof] = fit(times(1:maxlength)',tracks{m}.msd','power1');
            betas(m) = mdl.a; 
            alphas(m) = mdl.b; 
            rsqs(m) = gof.rsquare;
        elseif maxlength <= 1
            betas(m) = nan;
            alphas(m) = nan;
            rsqs(m) = nan;
        end
        torts(m) = tracks{m}.torts;
        speeds(m) = tracks{m}.speed;
        depths(m) = mean(tracks{m}.depth); %average depth of the track
        disps(m) = tracks{m}.currentdisp(end);
        if isfield(tracks, 'ventralmovement') == 1
        ventralmovement(m) = tracks{m}.ventralmovement;
        else
        end
    end
    
end