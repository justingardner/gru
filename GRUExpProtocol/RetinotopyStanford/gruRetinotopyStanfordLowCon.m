
% gruRetinotopyStanfordLowCon.m
% 
%                 date: 08/25/2017
%               author: minyoung lee
%              Purpose: fMRI mux8 sequences for retinotopy with low vs high
%                       contrast bars
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
%"phase per location": 493 vols (scanner will end slightly first)
% duration: 4:08 min
%

%% BARS
%Stimulus
%- always starts with a 12s blank (half 'stimulusPeriod=24' period)
%- interleave 3 other half period blanks ('blanks=3')
%- always interleave 8 24s period bar cycles

%mux 8 arc 1: 8 cycle * 48 vols/cycle + 3(+1 initial) blanks * 24vols = 480 vols + 8(mux)*2(nummux) = 493 vols (4:06min)
% Low contrast (15%)
mglRetinotopy('displayName=fMRIprojFlex','bars=1','barContrast=.035','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

% High contrast (100%)
mglRetinotopy('displayName=fMRIprojFlex','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');





