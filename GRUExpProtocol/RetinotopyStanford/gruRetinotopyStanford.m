
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
mglRetinotopy('displayName=fMRIprojFlex','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%% MT localizer`
Mtloc('0%',.5)






