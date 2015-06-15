
% gruRetinotopyStanford.m
% 
%                 date: 06/04/2015
%              Purpose: fMRI mux8 sequences for retinotopy
%                Usage: 
%
%                      copy paste each line in matlab
%
%Sequence info
%-------------
%
%MUX 8 ARC 1
%TR 0.5s
%7 slices per vol
%
%wedges and Rings
%----------------
%
%"phase per location": 496 vols
% duration: 4:07 min
%
%Bars
%----
%
%"phase per location": 472 vols
% duration: 3:55 min



%% Wedges and rings

%details:
%10 cycles, 24s/cycle, 48 vols/cycle, 2.4 mm iso vox.res, 0.5s TR, 480 vols + mux8*nummux = 496 vols
%10*48 + 16 = 496 vols

%CCW wedges
%mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% CW wedges
%mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% expanding rings
%mglRetinotopy('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% contracting rings
%mglRetinotopy('displayName=fMRIproj16','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% CCW wedges
%mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% CW wedges
%mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%% BARS

%details:
%for mux 8 arc 1 we set: same with 0.5s TR , 48 vols/cycle , 8*48 + 3*24 = 456 vols + mux8*nummux = 472 vols (3:55 min)
%mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');








