% pmMLR.m
%
%      usage: results = pmMLR(homedir, stimfile, datafile, stimradius)
%         by: justin gardner
%       date: 01/25/20
%    purpose: Wrapper for Gari's validation call
%       e.g.: 
%r = pmMLR('~/Desktop/temp','~/Desktop/stim.nii.gz','~/Desktop/bold.nii.gz',10);
%
function results = pmMLR(homedir, stimfile, datafile, stimradius)

% some debugging displays
doDispHDR = 1;
doDispFits = 1;

% also for debugging, do a quickfit
doQuickfit = 0;

% set up return value, default to failure
results.status = -1;

% check arguments
if ~any(nargin == [4]) || (nargout ~= 1)
  help pmMLR
  results.errorString = '(pmMLR) Invalid number of arguments';
  return
end

% make sure that MLR is in a default quit state
mrQuit;

% check for proper homedir
if ~isdir(homedir)
  mkdir(homedir);
end

% name for temporary MLR directory to process data
results.dataDir = fullfile(homedir,'pmMLR');

% remove any existing driectory
if isdir(results.dataDir) & ~isempty(results.dataDir);
  disp(sprintf('(pmMLR) Removing existing directory: %s',results.dataDir));
  system(sprintf('rm -rf %s',results.dataDir));
end

% now make that as a temporary directory
makeEmptyMLRDir(results.dataDir,'description=Temporary directory for pmMLR to process pRF data','defaultParams=1');

% get a view from this directory
cd(results.dataDir);
v = newView;

% check that bold data is valid
[d h] = mlrImageLoad(datafile);
if isempty(d)
  % error string
  results.errorString= sprintf('(pmMLR) Could not read %s',datafile);
  % do clean up and return
  results = cleanUp(results);
  return
end

% set qform and sfrom to 1 and identity
h.hdr.qformcode = 1;
h.hdr.sformcode = 1;
h.hdr.qform = eye(4);
h.hdr.sform = eye(4);

% save file as a Raw/TSeries
v = saveNewTSeries(v,d,[],h.hdr);

% check to see if we got a scan
if viewGet(v,'nScans') ~= 1
  % set error string
  results.errorString = sprintf('(pmMLR) Found %i scans in %s. Scanfile %s may not have loaded properly',viewGet(v,'nScans'),results.dataDir,datafile);
  % do clean up and return
  results = cleanUp(results);
  return
end

% set scanNum
scanNum = 1;


% load the stimulus image
[stim.im stim.hdr] = mlrImageLoad(stimfile);
if isempty(stim.im)
  % set error string
  results.errorString('(pmMLR) Could not read %s',stimfile);
  % do clean up and return
  results = cleanUp(results);
  return
end

% squeeze out any singleton dimensions from stim
stim.im = squeeze(stim.im);

% set the stimulus file, needs to be a structure with field
% stim.im, stim.x, stim.y, stim.t where stim.im is and w x h x t
% volume of stim images, stim.x and stim.y are w x h matrices
% with the x and y position in degrees and stim.t is a t array
% with the times of each frame

% the following are not necessary for the pRF code, but
% setting them in the struct for possible future ref
stim.size = size(stim.im);
stim.radius = stimradius;
% compute the X and Y of every position in the image
stimX = -stimradius:(2*stimradius)/(stim.size(1)-1):stimradius;
stimY = -stimradius:(2*stimradius)/(stim.size(2)-1):stimradius;
[stim.x stim.y] = meshgrid(stimX,stimY);
% set the t
framePeriod = viewGet(v,'framePeriod');
stim.t = 0:framePeriod:framePeriod*(stim.size(3)-1);

% setup the pRF analysis
[v params] = pRF(v,[],'justGetParams=1','defaultParams=1');

% add stimulus to params
params.pRFFit.stim{1} = stim;

if doQuickfit 
  params.pRFFit.quickPrefit = 1;
  params.pRFFit.prefitOnly = 1;
else
  params.pRFFit.quickPrefit = 0;
  params.pRFFit.prefitOnly = 0;
end
  
% do not try to split data (for running on sherlock)
params.pRFFit.splitData = 0;

% display verbose info
params.pRFFit.verbose = 1;

% parameters for canonical HRF using vistasoft code
params.pRFFit.canonicalFunction = 'rmHrfTwoGammas';
params.pRFFit.canonicalParams = [5.4 5.2 10.8 7.35 0.35];

% old way of generating difference of gammas
if 0
  params.pRFFit.canonicalFunction = 'getGammaHrf';
  params.pRFFit.timelag = 1;
  params.pRFFit.tau = 0.8;
  params.pRFFit.exponent = 6;
  params.pRFFit.diffOfGamma = 1;
  params.pRFFit.amplitudeRatio = 0.2;
  params.pRFFit.timelag2 = 4;
  params.pRFFit.tau2 = 2;
  params.pRFFit.exponent2 = 4;
end

% run it
v = pRF(v,params);

% get the overlays
overlayNames = {'r2','polarAngle','eccentricity','rfHalfWidth'};
% cycle over each one
for iOverlay = 1:length(overlayNames)
  % get which overlayNum it is
  overlayNum = viewGet(v,'overlayNum',overlayNames{iOverlay});
  % then read the data into results
  results.(overlayNames{iOverlay}) = viewGet(v,'overlayData',scanNum,overlayNum);
end

% convert to cartesian
results.x = cos(results.polarAngle) .* results.eccentricity;
results.y = sin(results.polarAngle) .* results.eccentricity;

% set status to good
results.status = 1;

% do clean up
results = cleanUp(results);

% display the hdr
if doDispHDR,dispHDR(v);end

% display fits 
if doDispFits,dispFits(results,stimradius);end

%%%%%%%%%%%%%%%%%
%    dispHDR    %
%%%%%%%%%%%%%%%%%
function dispHDR(v)

% get analysis strucuture
a = viewGet(v,'Analysis');
d = a.d{1};

% make figure with canonical plot
mlrSmartfig('pmMLR_hdr','reuse');clf;
plot(d.canonicalModel.time,d.canonicalModel.hrf);
xlabel('Time (s)');
ylabel('Response amplitude (% signal change)');
title('(pmMLR) Canonical hemodynamic response');

  
%%%%%%%%%%%%%%%%%%
%    dispFits    %
%%%%%%%%%%%%%%%%%%
function dispFits(results,stimradius)

mlrSmartfig('dispFits','reuse');clf;hold on
xaxis(-stimradius,stimradius);  
yaxis(-stimradius,stimradius);
grid on
for iVoxel = 1:length(results.x)
  % get a color
  c = getSmoothColor(iVoxel,length(results.x),'hsv');
  % plot the center
  plot(results.x(iVoxel),results.y(iVoxel),'o','MarkerFaceColor',c);
  % draw a circular RF
  x = []; y = [];
  for ang = 0:2*pi/359:2*pi;
    x(end+1) = results.x(iVoxel) + results.rfHalfWidth(iVoxel)*cos(ang);
    y(end+1) = results.y(iVoxel) + results.rfHalfWidth(iVoxel)*sin(ang);    
  end
  plot(x,y,'-','Color',c);
end
set(gca,'XTick',-stimradius:1:stimradius);
set(gca,'YTick',-stimradius:1:stimradius);

%%%%%%%%%%%%%%%%%
%    cleanUp    %
%%%%%%%%%%%%%%%%%
function results = cleanUp(results)

% display error string
if (results.status == -1) && ~isempty(results.errorString)
  disp(results.errorString);
end
% function to clean up after things have run
if isdir(results.dataDir)
  % remove everything in temporary directory
  system(sprintf('rm -rf %s',results.dataDir));
end

% remove the dataDir
results = rmfield(results,'dataDir');

% quit MLR and clean-up views
mrQuit;

