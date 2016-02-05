
% gruRetinotopyStanford.m
% 
%                 date: 09/09/2015
%               author: steeve
%              Purpose: fMRI mux8 sequences for retinotopy
%                Usage: 
%
%                      run each command (times X) in matlab
%
%Sequence info
%-------------
%
%32 channels coil
%MUX 8 ARC 1
%TR 0.5s
%7 slices per vol
%
%Bars
%----
%
%"phase per location": 493 vols
% duration: 4:08 min
%
%MT localizer
%------------
%
%"phase per location" : 548 vols
% duration: 4:33 min

%% BARS
%Stimulus
%- always starts with a 12s blank (half 'stimulusPeriod=24' period)
%- interleave 3 other half period blanks ('blanks=3')
%- always interleave 8 24s period bar cycles

%mux 8 arc 1: 8 cycle * 48 vols/cycle + 3(+1 initial) blanks * 24vols = 480 vols + 8(mux)*2(nummux) = 493 vols (4:06min)
%% Scan 1 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 2 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 3 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 4 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 5 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 6 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 7 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 8 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 9 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 10 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 11 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 12 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 13 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 14 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 15 (BarsTask)
mglRetinotopy('displayName=fMRIproj32','barsTask=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% Scan 16 (BarsTaskFixation)
mglRetinotopy('displayName=fMRIproj32','barsFix=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% ...





%%
% %% MT localizer
% % Mtloc('0%',.5)






