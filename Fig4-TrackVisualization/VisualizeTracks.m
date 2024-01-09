%Cell track visualization script

%load in files from individual birds and individual sessions (check BirdID and S columns of data in Processed Data or Processed Data/IndividualBirdData folder)

addpath ../ProcessedData
addpath ../Functions

load(); %input here
current = (); %input here

t = TemporallyColoredTracks(current)
d = DisplacementFromOrigin3D(current)
c = CheckHeadingBias(current)