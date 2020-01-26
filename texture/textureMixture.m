% texmix.m
%
%      usage: texmix()
%         by: justin gardner
%       date: 11/13/18
%    purpose: 
%
function retval = texmix(varargin)

% get arguments
getArgs(varargin,{'nScales=4','nOrientations=4','nSpatialNeighborhood=5','nIteration=50'});

% path to images
texturePath = '/Users/justin/Box Sync/Original_Images/bw_orig';
texture1 = 'rocks.jpg';
texture2 = 'drops.jpg';

texturePath = '/Users/justin/Box Sync/Original_Images/FS_orig';
texture1 = 'im13.png';
texture2 = 'im71.png';
texture1 = 'im56.png';
texture2 = 'im60.png';
texture1 = 'im18.png';
texture2 = 'im38.png';

%mixtures = [0 0.25 0.5 0.75 1];
mixtures = [0 0.25 0.5 0.75 1];
% load the two textures
texture1Data = imread(fullfile(texturePath,texture1));
texture2Data = imread(fullfile(texturePath,texture2));

% make grayscale
texture1Data = mean(texture1Data,3);
texture2Data = mean(texture2Data,3);

% get texture statistic
texture1Stats = textureAnalysis(texture1Data,nScales,nOrientations,nSpatialNeighborhood);
texture2Stats = textureAnalysis(texture2Data,nScales,nOrientations,nSpatialNeighborhood);

f2 = mlrSmartfig('textSynth','reuse');clf;

for iMix = 1:length(mixtures)
  % mix the two
  mixStats = mixTextures(mixtures(iMix),texture1Stats,texture2Stats);
  % synthesize image
  mixSynth(:,:,iMix) = textureSynthesis(mixStats,size(texture1Data),nIteration);
end

f = mlrSmartfig('texmix','reuse');clf;
for iMix = 1:length(mixtures)
  % display final image
  figure(f);
  subplot(1,length(mixtures),iMix);
  imagesc(mixSynth(:,:,iMix));
  colormap(gray);
  if mixtures(iMix) == 1
    title(sprintf('%s',texture1));
  elseif mixtures(iMix) == 0
    title(sprintf('%s',texture2));
  else
    title(sprintf('Mixture: %0.2f',mixtures(iMix)));
  end
end

keyboard
%%%%%%%%%%%%%%%%%%%%%
%    mixTextures    %
%%%%%%%%%%%%%%%%%%%%%
function mixStats = mixTextures(mixPercent,texture1Stats,texture2Stats)

statNames = fieldnames(texture1Stats);
for iStat = 1:length(statNames);
  % get combination of the two stats
  mixStats.(statNames{iStat}) = mixPercent * texture1Stats.(statNames{iStat}) + (1-mixPercent) * texture2Stats.(statNames{iStat});
end


