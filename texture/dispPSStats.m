% dispPSStats.m
%
%      usage: dispPSStats(imageFileName)
%         by: justin gardner
%       date: 10/19/18
%    purpose: Program to display Portilla and Simoncelli stats of a particular image. Requires
%             steerable pyramid and texture synthesis code from: http://www.cns.nyu.edu/~lcv/texture/
%
%             github repositories from Portilla & Simoncelli:
%             git clone https://github.com/LabForComputationalVision/textureSynth.git
%             git clone https://github.com/LabForComputationalVision/matlabPyrTools.git
%
%     usage: dispPSStats imageName.jpg
%
function retval = dispPSStats(imageFileName,varargin)

% check arguments
if nargin < 1
  help dispPSStats
  return
end

% get arguments
getArgs(varargin,{'nScales=4','nOrientations=4','nSpatialNeighborhood=7'});

% check for texture synth library
if isempty(which('textureAnalysis'))
  disp(sprintf('(dispPSStats) Portilla & Simoncelli textureSynth library missing'));
  disp(sprintf('              git clone https://github.com/LabForComputationalVision/textureSynth.git'));
  return
end

% check for matlab pyramid tools
if isempty(which('matlabPyrTools'))
  disp(sprintf('(dispPSStats) Simoncelli steerable pyramid library missing'));
  disp(sprintf('              git clone https://github.com/LabForComputationalVision/matlabPyrTools.git'));
  return
end

% check the file exist
if ~isfile(imageFileName)
  disp(sprintf('(dispPSStats) Could not find file: %s',imageFileName));
  return
end

% load file
imageData = imread(imageFileName);
if isempty(imageData)
  disp(sprintf('(dsipPSStats) Image file %s could not be read',imageFileName));
  return
end

% get image width and height
imageSize = size(imageData);
imageWidth = imageSize(2);
imageHeight = imageSize(1);
imageColorPlanes = size(imageData,3);

% if color image, then make it grayscale
imageDataGray = mean(imageData,3);

% run textureAnalysis to get Portilla and Simoncelli statistics
imageStats = textureAnalysis(imageDataGray,nScales,nOrientations,nSpatialNeighborhood);

% build pyramid representation of image
[imagePyramid,imagePyramidIndices] = buildSCFpyr(imageDataGray,nScales,nOrientations-1);
%[imagePyramid,imagePyramidIndices] = buildSFpyr(imageDataGray,nScales,nOrientations-1);

% initialize figure
figHist = mlrSmartfig('dispPSStats','reuse');clf;
subplotRowsHist = 3;subplotColsHist = 2;

% display original texture imaage
subplot(subplotRowsHist,subplotColsHist,1);
imagesc(imageData);
title(sprintf('Texture image: %s (%i x %i)',imageFileName,imageWidth,imageHeight));

% display pyramid bands
figPyramid = mlrSmartfig('dispPSStats_pyramid','reuse');clf;
subplotRows = nScales;subplotCols = nOrientations+1;
bandIndex = 2;
for iScale = 1:nScales
  scaleValues{iScale} = [];
  for iOrientation = 1:nOrientations
    % get the image band for the scale / orientation
    imageBand = pyrBand(imagePyramid, imagePyramidIndices, bandIndex);
    % displlay its magnitude
    subplot(subplotRows,subplotCols,(iScale-1)*subplotCols+iOrientation);
    imagesc(abs(imageBand));
    title(sprintf('Scale: %i Orientation: %i',iScale,iOrientation));
    % update band number
    bandIndex = bandIndex + 1;
    % keep all of the orientation images so that we can get histogram
    scaleValues{iScale} = [scaleValues{iScale} imageBand(:)];
  end  
end

% show high-pass residual
imageBand = pyrBand(imagePyramid, imagePyramidIndices, 1);
subplot(subplotRows,subplotCols,subplotCols);
imagesc(abs(imageBand));
title('High pass residual');

% show low-pass residual
imageBand = pyrBand(imagePyramid, imagePyramidIndices, nScales*nOrientations+2);
subplot(subplotRows,subplotCols,(subplotRows-1)*subplotCols+subplotCols);
imagesc(abs(imageBand));
title('Low pass residual');

colormap(gray);

% display histogram
figure(figHist);
subplot(subplotRowsHist,subplotColsHist,2);
dispHist(imageData(:),imageStats.pixelStats,'');

for iScale = 1:nScales
  subplot(subplotRowsHist,subplotColsHist,2+iScale);
  dispHist(abs(scaleValues{iScale}(:)),imageStats.pixelLPStats(iScale,:),sprintf('Scale %i: ',iScale));
end

keyboard

%%%%%%%%%%%%%%%%%%
%    dispHist    %
%%%%%%%%%%%%%%%%%%
function dispHist(histData,pixelStats,titleStr)

% calculate min and max of distribution
minHist = double(min(abs(histData)));
maxHist = double(max(abs(histData)));

binSize = (maxHist-minHist)/32;
if (binSize>1),binSize = round(binSize);end

bins = minHist:binSize:maxHist;

% display histogram
myhist(histData,bins);

% set to 0 to 255 if that is the image range
if (maxHist > 128) && (maxHist <= 255) && (minHist >=0)
  xaxis(0,255);
  set(gca,'XTick',[0 64 128 192 255]);
end

% set title to display pixel statistics
if length(pixelStats) == 6
  title(sprintf('%sMean: %0.2f Var: %0.2f Skew: %0.2f Kurtosis: %0.2f Minmax: [%0.1f %0.1f]',titleStr,pixelStats(1),pixelStats(2),pixelStats(3),pixelStats(4),pixelStats(5),pixelStats(6)));
else
  title(sprintf('%sSkew: %0.2f Kurtosis: %0.2f',titleStr,pixelStats(1),pixelStats(2)));
end  
ylabel('n');
xlabel('Pixel value');

