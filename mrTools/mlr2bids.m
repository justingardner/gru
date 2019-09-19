function mlr2bids(group,subj)
%% MLR2BIDS Convert an mrTools folder into a BIDS folder
%
% INPUTS
%   folder: 'group/subj', for example: 'mglRetinotopy/s035920180808'
% OUTPUTS
%   folder/bids/soft-links and folder/bids/mriqc
%
% The mrTools format is:
% s0###YYYYMMDD
%  - Anatomy/
%    - *.nii
%  - Averages/
%  - Concatenation/
%  - Etc/
%  - MotionComp/
%  - Pre/
%  - ROIs/
%  - Raw/
%    - *.nii
%
% We don't care about anything except the Anatomy and Raw folders
%
% The BIDS format is:
% sub-s0###YYYYMMDD
%  - anat/
%    - *.nii (or *.nii.gz)
%  - func/
%    - *.nii (or *.nii.gz)
%
% To do this we will create symbolic links and then call the mriqc
% functions. 

prefix = '~/data/';
bsubj = 'sub-S01';

mFolder = fullfile(prefix,group,subj);
bFolder = fullfile(prefix,group,subj,'bids',bsubj,'ses-1');

if ~isdir(bFolder)
    mkdir(bFolder);
end
if ~isdir(fullfile(bFolder,'anat'))
    mkdir(fullfile(bFolder,'anat'));
end
if ~isdir(fullfile(bFolder,'func'))
    mkdir(fullfile(bFolder,'func'));
end

%% Collect information about the relevant files
anat = dir(fullfile(mFolder,'Anatomy','*.nii*'));
bold = dir(fullfile(mFolder,'Raw','*.nii*'));

%% Make new MRIQC files
for ai = 1:length(anat)
    symbolic(fullfile(mFolder,'Anatomy',anat(ai).name),fullfile(bFolder,'anat',sprintf('%s_%s',bsubj,anat(ai).name)));
end

for bi = 1:length(bold)
    symbolic(fullfile(mFolder,'Raw',anat(ai).name),fullfile(bFolder,'func',sprintf('%s_%s',bsubj,anat(ai).name)));
end

function succ = symbolic(ofile,nfile)

command = sprintf('ln -s %s %s',ofile,nfile);
system(command);

succ = isfile(nfile);