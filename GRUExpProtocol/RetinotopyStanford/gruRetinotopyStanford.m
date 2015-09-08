
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
%"phase per location": 500 vols
% duration: 4:09 min
%
%MT localizer
%------------
%
%"phase per location" : 548 vols
% duration: 4:33 min

%% BARS
%for mux 8 arc 1 we should set: 48 vols/cycle, 8*48 vols + 3 blanks * 24 s = 456 vols + 8(mux)*2(nummux) = 472 vols (3:55 min)
%but 500 actually works best (need to findout why)
mglRetinotopy('displayName=fMRIproj32','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% MT localizer
Mtloc('0%',.5)






