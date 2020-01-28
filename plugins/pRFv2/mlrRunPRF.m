% mlrRunPRF.m
%
%      usage: results = mlrRunPRF(homedir, datafile, stimfile, stimsize)
%         by: justin gardner
%       date: 01/25/20
%    purpose: Wrapper for Gari's validation call. This is meant to be run as a standalone function
%             to run the MLR pRF code. It is written so that it can be dockarized.
%
%             homedir: A string that is the directory where this function will temporarily
%                      create (and then delete) an MLR direcotry for processing
%
%             datafile: A string that contains the filename of the datafile. This should be a nifti
%                       of BOLD data. The nifti header should have correct info - in particular the
%                       framePeriod (time per each volume) will be read off the nifti header which is important for
%                       fitting the model as the hemodynamic impulse response needs to be calculated
%                       correctly.
%
%             stimfile: A string that contains the filename of a nifti file which contains
%                       the stimulus. This should be essentially a mask of zeros and ones
%                       that indicate timeframe by timeframe where the stimulus is. Note
%                       that this should have the same number of frames as the datafile. Each frame
%                       is an image with width and height determined by the stimsize argument.
%
%             stimsize: An array [width height] specifying the size of the stimulus image in degrees of visual angle
%                       If this is a scalar, then will assume a square stimulus image where
%                       width = height = stimsize.
%
%             You can set optional arguments (these are set by adding the argument as a string followed by
%             an argument for the value: e.g. ...,'dispFigs=1',1):
%
%             dispFigs: (default 0) display figures with fits
%             quickFit: (default 0) just do a quick prefit to test code
%             parallel: (default 1) Set to number of parallel workers to startup, 0 to not use parallel toolbox
%               set to 1 to bring up a dialog box that will ask to set workers
%             canonicalParams: (Default [5.4 5.2 10.8 7.35 0.35]) array of parameters for canonical difference of
%               gamma function [peak1 fwhm1 peak2 fwhm2 dip]. See rmHrfTwogammas in Vistasoft for more info.
%             rfType: (default 'gaussian') The name of the pRF model we are fitting. See the models directory
%               of pRFv2 to see models that have been codedpR
%
%             Return variable (results) has the following fields
%
%             status: -1 for failure, 1 for success
%             errorString: only set for failure - a string that gives info on what the error was
%             r2: r2 of pRF fit for each voxel
%             polarAngle, eccentricity: RF location for each voxel in polar coordinates
%             x, y: RF location for each voxel in cartesian coordinates
%             rfHalfWidth: Std of gaussian fit for each voxel
%
%       e.g.: 
%r = mlrRunPRF('~/Desktop/temp','~/Desktop/bold.nii.gz','~/Desktop/stim.nii.gz',[20 20],'dispFigs=1');
%
function results = mlrRunPRF(homedir, datafile, stimfile, stimsize, varargin)

% set up return value (0 means in progress, -1 failure, 1 success)
results.status = 0;

% check arguments
if (nargin < 4) || (nargout ~= 1)
  help mlrRunPRF
  results.status = -1;
  results.errorString = '(mlrRunPRF) Invalid number of arguments';
  return
end

% check and interpret arguments, load stim and data
[results inputs flags] = checkArguments(homedir, datafile, stimfile, stimsize, varargin, results);
if results.status == -1,return,end

% setup MLR session and get a view
[results v] = setupMLR(inputs, results);
if results.status == -1,return,end

% setup stimulus image
stim = setupStimulus(inputs);

% get the parameters for running the pRF analysis
[v params] = getPRFParams(v,stim,inputs,flags);

% and run the pRF analysis
v = pRF(v,params);

% extract the overlays and analysis into the results
results = getAnalysisResults(v,results);

% set status to good
results.status = 1;

% do clean up
results = cleanUp(results,inputs);

% display the hdr if set to do that
if flags.doDispHDR,dispHDR(v);end

% display fits if set to do that
if flags.doDispFits,dispFits(results,inputs);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getAnalysisResults    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = getAnalysisResults(v,results)
  
% what scan to get data from (should always be one)
scanNum = viewGet(v,'curScan');
  
% get the overlays
overlayNames = viewGet(v,'overlayNames');

% cycle over each one, extracting from the view
for iOverlay = 1:length(overlayNames)
  % get which overlayNum it is
  overlayNum = viewGet(v,'overlayNum',overlayNames{iOverlay});
  % then read the data into results
  results.(overlayNames{iOverlay}) = viewGet(v,'overlayData',scanNum,overlayNum);
end

% convert to cartesian
results.x = cos(results.polarAngle) .* results.eccentricity;
results.y = sin(results.polarAngle) .* results.eccentricity;

%%%%%%%%%%%%%%%%%%%%%%
%    getPRFParams    %
%%%%%%%%%%%%%%%%%%%%%%
function [v params] = getPRFParams(v,stim,inputs,flags);

% setup the pRF analysis
[v params] = pRF(v,[],'justGetParams=1','defaultParams=1');

% add stimulus to params
params.pRFFit.stim{1} = stim;

% setup parameters for doing a quick analysis (for debugging)
if flags.doQuickfit 
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
params.pRFFit.canonicalParams = inputs.canonicalParams;

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

% set model name
params.pRFFit.rfType = inputs.rfType;

%%%%%%%%%%%%%%%%%%%%%%%
%    setupStimulus    %
%%%%%%%%%%%%%%%%%%%%%%%
function stim = setupStimulus(inputs);
 
% squeeze out any singleton dimensions from stim
stim.im = squeeze(inputs.stimulusImage);

% set the stimulus file, needs to be a structure with field
% stim.im, stim.x, stim.y, stim.t where stim.im is and w x h x t
% volume of stim images, stim.x and stim.y are w x h matrices
% with the x and y position in degrees and stim.t is a t array
% with the times of each frame

% the following are not necessary for the pRF code, but
% setting them in the struct for possible future ref
stim.size = size(stim.im);
stim.width = inputs.stimulusWidth;

% compute the X and Y of every position in the image
stimX = -inputs.stimulusWidth/2:(inputs.stimulusWidth)/(stim.size(1)-1):inputs.stimulusWidth/2;
stimY = -inputs.stimulusHeight/2:(inputs.stimulusHeight)/(stim.size(1)-1):inputs.stimulusHeight/2;
[stim.x stim.y] = meshgrid(stimX,stimY);
% set the t
stim.t = 0:inputs.framePeriod:inputs.framePeriod*(stim.size(3)-1);

%%%%%%%%%%%%%%%%%%
%    setupMLR    %
%%%%%%%%%%%%%%%%%%
function [results v] = setupMLR(inputs, results);

% make sure that MLR is in a default quit state
mrQuit;

% now make that as a temporary directory
makeEmptyMLRDir(inputs.dataDir,'description=Temporary directory for mlrRunPRF to process pRF data','defaultParams=1');

% get a view from the temporary directory
cd(inputs.dataDir);
v = newView;

% save bold file as a Raw/TSeries
v = saveNewTSeries(v,inputs.boldData,[],inputs.boldHeader.hdr);

% check to see if we got a scan
if viewGet(v,'nScans') ~= 1
  % set error
  results.status = -1;
  % set error string
  results.errorString = sprintf('(mlrRunPRF:setupMLR) Found %i scans in %s. Scanfile may not have loaded properly',viewGet(v,'nScans'),inputs.dataDir);
  % do clean up and return
  results = cleanUp(results,inputs);
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    checkArguments    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [results inputs flags] = checkArguments(homedir, datafile, stimfile, stimsize, varargin, results);

% check arguments
getArgs(varargin,{'dispFigs=0','quickFit=0','parallel=1','canonicalParams=[5.4 5.2 10.8 7.35 0.35]','rfType=gaussian'},'verbose=1');

% deal with parallel workers for parfor loops
global mlrNoParallel;
if parallel == 0
  mlrNoParallel = 1;
elseif parallel >= 1
  mlrNoParallel = 0;
  mlrNumWorkers(parallel);
end
    
% some flags for different debugging displays
flags.doDispHDR = dispFigs;
flags.doDispFits = dispFigs;

% also for debugging, do a quickfit
flags.doQuickfit = quickFit;

% check for proper homedir
if ~isdir(homedir)
  mkdir(homedir);
end

% name for temporary MLR directory to process data
inputs.dataDir = fullfile(homedir,'mlrRunPRF');

% remove any existing driectory
if isdir(inputs.dataDir) & ~isempty(inputs.dataDir);
  disp(sprintf('(mlrRunPRF:checkArguments) Removing existing directory: %s',inputs.dataDir));
  system(sprintf('rm -rf %s',inputs.dataDir));
end

% check canonical params
if length(canonicalParams) ~= 5
  % set error
  results.status = -1;
  % set error string
  results.errorString = '(mlrRunPRF:checkArguments) Canonical parameters should be an array of length 5: [peak1 fwhm1 peak2 fwhm2 dip] see rmHrfTwogammas in Vistasoft for more info';
  % do clean up and return
  results = cleanUp(results,inputs);
  return
end
inputs.canonicalParams = canonicalParams;

% load the BOLD data and check
[inputs.boldData inputs.boldHeader] = mlrImageLoad(datafile);
if isempty(inputs.boldData)
  % set error
  results.status = -1;
  % error string
  results.errorString= sprintf('(mlrRunPRF:checkArguments) Could not read %s',datafile);
  % do clean up and return
  results = cleanUp(results,inputs);
  return
end

% set qform and sfrom to 1 and identity
inputs.boldHeader.hdr.qformcode = 1;
inputs.boldHeader.hdr.sformcode = 1;
inputs.boldHeader.hdr.qform = eye(4);
inputs.boldHeader.hdr.sform = eye(4);

% check frame period
if length(inputs.boldHeader.pixdim) >= 4
  inputs.framePeriod = inputs.boldHeader.pixdim(4);
  disp(sprintf('(mlrRunPRF:checkArguments) Frame period of BOLD scan is: %0.2fs',inputs.framePeriod));
else
  % set error
  results.status = -1;
  % set error string
  results.errorString('(mlrRunPRF:checkArguments) Could not determine frame period of BOLD scan because pixdim array does not have the 4th dimension set');
  % do clean up and return
  results = cleanUp(results,inputs);
  return
end  

% load the stimulus image
[inputs.stimulusImage inputs.stimulusHeader] = mlrImageLoad(stimfile);
if isempty(inputs.stimulusImage)
  % set error
  results.status = -1;
  % set error string
  results.errorString('(mlrRunPRF:checkArguments) Could not read %s',stimfile);
  % do clean up and return
  results = cleanUp(results,inputs);
  return
end

% check stimsize
if ~any(length(stimsize) == [1 2])
  % set error
  results.status = -1;
  % set error string
  results.errorString('(mlrRunPRF:checkArguments) stimsize should be a scalar or a vector of two elements');
  % do clean up and return
  results = cleanUp(results,inputs);
  return
elseif length(stimsize) == 1
  % set stimulus width and height
  inputs.stimulusWidth = stimsize;
  inputs.stimulusHeight = stimsize;
else
  inputs.stimulusWidth = stimsize(1);
  inputs.stimulusHeight = stimsize(1);
end

% keep direcgtory where we started
inputs.startDir = cd;

% set rfType (that is the model that we are fitting)
inputs.rfType = rfType;

%%%%%%%%%%%%%%%%%
%    dispHDR    %
%%%%%%%%%%%%%%%%%
function dispHDR(v)

% get analysis strucuture
a = viewGet(v,'Analysis');
d = a.d{1};

% make figure with canonical plot
mlrSmartfig('mlrRunPRF_hdr','reuse');clf;
plot(d.canonicalModel.time,d.canonicalModel.hrf);
xlabel('Time (s)');
ylabel('Response amplitude (% signal change)');
title('(mlrRunPRF) Canonical hemodynamic response');

  
%%%%%%%%%%%%%%%%%%
%    dispFits    %
%%%%%%%%%%%%%%%%%%
function dispFits(results,inputs)

mlrSmartfig('mlrRunPRF_dispFits','reuse');clf;hold on
xaxis(-inputs.stimulusWidth/2,inputs.stimulusWidth/2);  
yaxis(-inputs.stimulusHeight/2,inputs.stimulusHeight/2);
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
set(gca,'XTick',-inputs.stimulusWidth/2:1:inputs.stimulusWidth/2);
set(gca,'YTick',-inputs.stimulusHeight/2:1:inputs.stimulusHeight/2);

%%%%%%%%%%%%%%%%%
%    cleanUp    %
%%%%%%%%%%%%%%%%%
function results = cleanUp(results,inputs)

% display error string
if (results.status == -1) && ~isempty(results.errorString)
  disp(results.errorString);
end
% function to clean up after things have run
if isdir(inputs.dataDir)
  % remove everything in temporary directory
  system(sprintf('rm -rf %s',inputs.dataDir));
end

% quit MLR and clean-up views
mrQuit;

% switch back to staruting directory
if isfield(inputs,'startDir') & isdir(inputs.startDir)
  cd(inputs.startDir);
end
