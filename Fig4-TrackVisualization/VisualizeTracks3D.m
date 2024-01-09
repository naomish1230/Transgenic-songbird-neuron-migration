%Cell track visualization script

%load in files from individual birds and individual sessions (check BirdID and S columns of data in Processed Data or Processed Data/IndividualBirdData folder)

addpath ../ProcessedData
addpath ../Functions
load('grey3L_dynamic_stats.mat'); %load dynamic stats file for bird, rename it to data
data = grey3L_dynamic_stats; %need to check this, change to 1 session at a time if length(unique(data.S))>1 [ex: (pink61_dynamic_stats(pink61_dynamic_stats.S == 1,:)]

%make 3D spatial plot of all cells' positions over time, color coded by time
t = TemporallyColoredTracks(data);

%make 3D plot of relative displacements from origin, color coded by time
d = DisplacementFromOrigin3D(data,1); %input2: 1 = plot red dots at ends of cells, 0 = no dot plotting