function [mdl, avgMSD, times, alpha, rsq] = ComputeMSD_function(current)
%USED in PopulationSuperdiffusivity.m

%calculate MSD of population using online formula found here:
%https://measurebiology.org/wiki/Calculating_MSD_and_Diffusion_Coefficients
%to find MSD and diffusion coefficient
%% ======== OPTIONS

dimensions = 2; %set to number of dimensions you want to analyze. if 2, will look at XY. if 3, will look at XYZ; '2' used for paper due to limited axial resolution in Z.

format shortg %preference; 'shortg' for paper
%% ASSIGN APPROPRIATE FRAME RATES TO DIFFERENT BIRDS

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
%% COMPUTE MSD OF EACH TRACK if same frame rate
longestTrack = 0;
cellsToSkip = [];
cells = unique(current.N);


numberOfCells = length(cells); 
tracks = cell(1,numberOfCells);
framerates = [];
    for ii = 1:numberOfCells %cycle through cells 
        nID = cells(ii);
        track = current(current.N == nID,:); %isolate one cell's entire track 
        clear pathX pathY pathZ paths_3d_cell hourdiffs currentdispXYZ currentdisp
        tc = 0; %set counter
        for t = 2:length(track.F); %cycle through that cell's track over time
           tc = tc+1; %add to counter each loop
           pathX(tc) = (track.X(t) - track.X(t-1)); %differences between current and previous X position
           pathY(tc) = (track.Y(t) - track.Y(t-1)); %differences between current and previous Y position
           if dimensions == 3
            pathZ(tc) = (track.Z(t) - track.Z(t-1)); %differences between current and previous Z position
            paths_3d_cell(tc) = sqrt(pathX(tc)^2 + pathY(tc)^2 + pathZ(tc)^2); %Euclidean 3D distance between current and previous position
            currentdispXYZ(tc,:) = [(track.X(t) - track.X(1)), (track.Y(t) - track.Y(1)),(track.Z(t) - track.Z(1))]; %XYZ distance between current and first position
            currentdisp(tc) = sqrt(sum(currentdispXYZ(tc,1)^2 + currentdispXYZ(tc,2)^2 + currentdispXYZ(tc,3)^2)); %Euclidean 3D distance between current and first position
           elseif dimensions == 2
               paths_3d_cell(tc) = sqrt(pathX(tc)^2 + pathY(tc)^2); %Euclidean 2D distance between current and previous position
               currentdispXYZ(tc,:) = [(track.X(t) - track.X(1)), (track.Y(t) - track.Y(1))]; %XYZ distance between current and first position
               currentdisp(tc) = sqrt(sum(currentdispXYZ(tc,1)^2 + currentdispXYZ(tc,2)^2)); %Euclidean 2D distance between current and first position
           end
           hourdiffs(tc) = ((track.H(t))-(track.H(t-1))); %hourly difference
       end
        numberOfSamples = size(track,1);
        tracks{ii} = struct(); 
        numberOfMsds = round(numberOfSamples-1);
        tracks{ii}.msd = zeros(1, numberOfMsds);
        tracks{ii}.msdUncertainty = zeros(1, numberOfMsds);
        tracks{ii}.number = numberOfMsds;
        tracks{ii}.cellID = nID; 
        tracks{ii}.framerate = current.framerate(current.N == nID);
        tracks{ii}.frames = track.F;
        tracks{ii}.torts = sum(paths_3d_cell)./currentdisp(end); %tortuosity
        tracks{ii}.speed = sum(paths_3d_cell)/(max([track.H]) - min([track.H])); %speed
        tracks{ii}.stepsize = paths_3d_cell(:); %step distances
        if ~ismember(unique(track.Bird_ID), [4,5,11]) %if not one of the prism birds, then...
            tracks{ii}.depth = track.Z;
        elseif unique(track.Bird_ID) == 4
            tracks{ii}.depth = track.X + 305.848; %hard coded because computed from raw image data and distance to brain surface. do not change.
        elseif unique(track.Bird_ID) == 5
            tracks{ii}.depth = track.Y + 167.201; %hard coded because computed from raw image data and distance to brain surface. do not change.
        elseif unique(track.Bird_ID) == 11
            tracks{ii}.depth = track.Y + 81.568; %hard coded because computed from raw image data and distance to brain surface. do not change.
        end
        
        framerates(ii) = tracks{ii}.framerate(1); %index framerates to find minimum later

        if numberOfMsds==0
            cellsToSkip = [cellsToSkip ii];
        else
            track = table2array(track(:,3:5)); %turns track into array
            for tau = 1:numberOfMsds %finds the MSD across different time intervals tau
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

            tracks{ii}.diffusionCoefficient = tracks{ii}.msd(1) / (2 * dimensions * tracks{ii}.framerate(1)^-1); %2*dimension*time
            
            longestTrack = max( longestTrack, numberOfSamples );
        end
    end

% compute the average value of D and the ensemble average of msd
d = zeros(1,numberOfCells);
averageMsd = zeros( 1, longestTrack );
averageMsdCount = zeros( 1, longestTrack);
 
for ii = 1:numberOfCells
    if sum(cellsToSkip==ii)==0 
        numberOfMsds = length(tracks{ii}.msd);
        d(ii) = tracks{ii}.diffusionCoefficient;
        averageMsd(1:numberOfMsds) = averageMsd(1:numberOfMsds) + tracks{ii}.msd;
        averageMsdCount(1:numberOfMsds) = averageMsdCount(1:numberOfMsds) + 1;
    end
end
 
minFrameRate = min(framerates); %not optimal if combining birds across vastly different timescales, not all birds can be compared with each other
averageMsd = averageMsd ./ averageMsdCount;


%% finding the equations of best fit 

%shorten ensemble MSD to before the peak, it drops off because MSD is still averaged for missing and short cell tracks
maxlength = find(averageMsd(:) == max(averageMsd)) - 1; 

times = (1:longestTrack)/minFrameRate;

%Fit power model to alpha
[mdl,gof] = fit(times(1:maxlength)',averageMsd(1:maxlength)','power1');

alpha = mdl.b;
rsq = gof.adjrsquare;


    avgMSD = averageMsd(1:maxlength);
    times = times(1:maxlength);

end
