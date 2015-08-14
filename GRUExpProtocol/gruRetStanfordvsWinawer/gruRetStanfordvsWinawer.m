
% gruRetStanfordvsWinawer.m
% 
%                 date: 08/07/2015
%               author: steeve
%              Purpose: fMRI mux8 sequences. Compare Winawer retinotopy
%                       with ours.
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
%
%Bars
%----
%
%"phase per location": 500 vols
% duration : 4:09 min


%% BARS

%details:
%for mux 8 arc 1 we set: same with 0.5s TR , 48 vols/cycle , 8*48 + 3*24 = 456 vols + mux8*nummux = 472 vols, but stim seems to last 
%but a bit longer so 500 vols (4 min)

%% width 1
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=1','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=1','doEyeCalib=0');

%% width 2
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=2','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=2','doEyeCalib=0');

%% width 3
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=3','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=3','doEyeCalib=0');

%repeat the whole


%% width 1
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=1','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=1','doEyeCalib=0');

%% width 2
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=2','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=2','doEyeCalib=0');

%% width 3
mglRetinotopy('displayName=fMRIproj16','bars','barsTaskDefault','fixedRandom=1','barWidth=3','doEyeCalib=0');

%%
mglRetinotopy('displayName=fMRIproj16','barsTask','barsTaskDefault','fixedRandom=1','barWidth=3','doEyeCalib=0');



