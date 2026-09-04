% pRFGetStimImageFromStimfile.m
%
%        $Id:$
%      usage: stim = pRFGetStimImageFromStimfile(stimfile,<timePoints>)
%         by: justin gardner
%       date: 10/17/11
%    purpose: Pass in a stimfile (can be either a string filename, or a strucutre
%             with myscreen/task) created by mglRetinotopy (make sure this is a
%             stimfile created with a version of mglRetinotopy past 10/2011 which
%             has the proper variables stored to enable reconstruction). Will
%             create a volume of dimensions x,y,t with the stimulus image (load
%             in mlrVol to view). stim.x and stim.y are the X and Y coordinates
%             in degrees of every point. stim.t is the array of times at
%             which image is taken.
%
%             Optionally arguments:
%
%             timePoints: array for which the stim image should be computed.
%
%             The default returned stimulus is in the canonical world-coordinate
%             orientation used by pRF. A saved pRFStimImage must use this same
%             orientation. xFlip and yFlip are optional, explicit reversals of
%             that canonical image.
%
%             compareStimfiles: return true when the supplied stimfiles have
%             the same saved pRFStimImage or the same stimulus/task settings.
%             This check does not open an MGL display or construct images.
%
%             saveStimImageTargets: internal pRFFit option. When a verified
%             averaged scan is reconstructed from one component, save the
%             same canonical pRFStimImage to every supplied target stimfile.
%
%             Note for developers - this function needs to keep up-to-date with
%             any changes in the display loop of mglRetinotopy to interpret
%             the stimfiles correctly
%
function stim = pRFGetStimImageFromStimfile(stimfile,varargin)

% set default return arguments
stim = [];

% check arguments
if nargin < 1
  help pRFGetStimImageFromStimfile
  return
end

% parse arguments
timePoints = [];screenWidth = [];screenHeight = [];volTrigRatio = [];
xFlip = [];yFlip = [];timeShift = [];verbose = [];debugStimulus = [];compareStimfiles = [];
getArgs(varargin,{'timePoints=[]','screenWidth=[]','screenHeight=[]','volTrigRatio=[]','xFlip=0','yFlip=0','timeShift=0','verbose=0','debugStimulus=0','saveStimImage=0','recomputeStimImage=0','compareStimfiles=0','saveStimImageTargets=[]'});

if compareStimfiles
  stim = haveSameStimulus(stimfile,volTrigRatio,recomputeStimImage);
  return
end

% handle cell array
if iscell(stimfile) && ((length(stimfile)>1) || (length(stimfile{1})>1))
  for i = 1:length(stimfile)
    % get current volTrigRatio
    if isempty(volTrigRatio)
      thisVolTrigRatio = [];
    else
      thisVolTrigRatio = volTrigRatio{i};
    end
    stim{i} = pRFGetStimImageFromStimfile(stimfile{i},'timePoints',timePoints,'screenWidth',screenWidth,'screenHeight',screenHeight,'volTrigRatio',thisVolTrigRatio,'xFlip',xFlip,'yFlip',yFlip,'timeShift',timeShift,'verbose',verbose,'debugStimulus',debugStimulus,'saveStimImage',saveStimImage,'recomputeStimImage',recomputeStimImage);
    if isempty(stim{i}),stim = [];return;end
  end
  return
end

% check volTrigRatio
if iscell(volTrigRatio)
  if length(volTrigRatio) > 1
    disp(sprintf('(pRFGetStimImageFromStimfile) volTrigRatio should not be of length greater than one (length=%i) using only the first value of %i',length(volTrigRatio),volTrigRatio{1}));
  end
  volTrigRatio = volTrigRatio{1};
end

% get the stimfile
s = getStimfile(stimfile);
if isempty(s),return,end

% Force reconstruction from the task information, not from an embedded image.
% Removing the field here also lets checkStimfile recover taskNum correctly.
if recomputeStimImage && isstruct(s) && isfield(s,'pRFStimImage')
  s = rmfield(s,'pRFStimImage');
end

% check that we have a stimfile that is interpretable
% by this program
[tf s taskNum] = checkStimfile(s);
if ~tf,return,end

% check to see if a stimImage exists
if ~isfield(s,'pRFStimImage') || recomputeStimImage
  % make the traces
  s.myscreen = makeTraces(s.myscreen,verbose);

  % get task variables
  e = getTaskParameters(s.myscreen,s.task{taskNum});

  % get some traces of things of interest
  s.time = s.myscreen.time;
  s.maskPhase = s.myscreen.traces(s.task{taskNum}{1}.maskPhaseTrace,:);
  s.blank = s.myscreen.traces(s.task{taskNum}{1}.blankTrace,:);
  s.vol = s.myscreen.traces(1,:);
  s.trialVol = e.trialVolume;
  s.trialTicknum = e.trialTicknum;
  s.blank = e.randVars.blank;
  switch s.stimulus.stimulusType
    case {1,2}
      % Wedge and ring stimuli do not use bar-angle state.
    case {3,4,5}
      s.barAngle = e.parameter.barAngle;
      s.elementAngle = e.randVars.elementAngle;
    otherwise
      error('pRFGetStimImageFromStimfile:UnsupportedStimulusType', ...
        'Unsupported stimulus type: %s',mat2str(s.stimulus.stimulusType));
  end

  % if no timepoints, then get one for each volume
  if isempty(timePoints)
    timePoints = s.time(find(s.vol))-s.time(first(find(s.vol)));

    % add timepoints for sense acceleration
    if ~isempty(volTrigRatio)
      framePeriod = median(diff(timePoints));
      accTimePoints = [];
      for i = 1:length(timePoints);
      	accTimePoints(end+1) = timePoints(i);
      	for j = 1:(volTrigRatio-1)
      	  accTimePoints(end+1) = timePoints(i)+framePeriod*j/volTrigRatio;
      	end
      end
      timePoints = accTimePoints;
    end
    % add one frame extra since the stimulus usually changes after
    % the volume
    timePoints = timePoints + s.myscreen.frametime;
  end
  stim.t = timePoints;

  % check the pixel dimensions
  if ~isempty(screenWidth) && ~isempty(screenHeight)
  elseif isempty(screenWidth) && isempty(screenHeight)
    screenWidth = round(s.myscreen.screenWidth/10);
    screenHeight = round(s.myscreen.screenHeight/10);
  elseif isempty(screenWidth)
    % no screenWith, compute based on aspect ratio
    screenWidth = round(screenHeight*s.myscreen.screenWidth/s.myscreen.screenHeight);
  elseif isempty(screenHeight)
    screenHeight = round(screenWidth*s.myscreen.screenHeight/s.myscreen.screenWidth);
  end

  % open the screen
  % if mglGetParam('displayNumber') ~= -1,mglClose;end
  % mglSetParam('offscreenContext',1);
  mglOpen(0,screenWidth,screenHeight);
  stimulusDisplayCleanup = onCleanup(@closeStimulusDisplay);
  mglFrameGrab('init');

  mglVisualAngleCoordinates(s.myscreen.displayDistance,s.myscreen.displaySize);

  % create stim.x and stim.y
  imageWidth = s.myscreen.imageWidth;
  imageHeight = s.myscreen.imageHeight;
  [stim.x stim.y] = ndgrid(-imageWidth/2:imageWidth/(screenWidth-1):imageWidth/2,-imageHeight/2:imageHeight/(screenHeight-1):imageHeight/2);

  if verbose,disppercent(-inf,'(pRFGetStimImageFromStimfile) Computing stimulus images');end
  warnOnStimfileMissingInfo = true;
  for iImage = 1:length(stim.t)
    im = createMaskImage(s,stim.t(iImage),debugStimulus);
    % if no image, that probably means the stimfile ended early
    if isempty(im)
      % check to see if we are within (arbitrarily) 5% of the end
      % and use that as a cutoff for asking the user if something drastically
      % wrong has occurred.
      if warnOnStimfileMissingInfo
      	if (1-iImage/length(stim.t)) > 0.05
      	  if askuser('Your stimfile is missing information for volume %i of %i. This might be because you have linked the wrong stimfile or that the stim program ended before the scan or that there is some other problem with the stimfile. It would be a good idea to try to fix this cause this may now be generating the wrong stimulus. Continue anyway?',0,1)
      	    warnOnStimfileMissingInfo = false;
      	  else
      	    % user did not agree to continue, bail out
            stim=[];
      	    return
      	  end
      	end
      end
      % just put up the warning
      disp(sprintf('(pRFGetStimImageFromStimfile) !!! Missing stimulus info for volume %i of %i. Setting to blank image. !!!',iImage,length(stim.t)));
      im = zeros(screenWidth,screenHeight);
    end
    % Allocate from the first reconstructed frame so the movie retains the
    % exact numeric class that indexed growth would have produced.
    if iImage == 1
      stim.im = zeros(screenWidth,screenHeight,length(stim.t),'like',im);
    end
    stim.im(1:screenWidth,1:screenHeight,iImage) = im';
    if verbose,disppercent(iImage/length(stim.t));end
  end
  if verbose,disppercent(inf);end

  % Close the offscreen display now; onCleanup also covers every early return
  % and error after mglOpen.
  clear stimulusDisplayCleanup

  % mglFrameGrab raster order is opposite the ascending Cartesian y grid in
  % stim.y. Normalize that representation once during reconstruction so the
  % public return value and any saved pRFStimImage share one canonical
  % world-coordinate orientation.
  stim.im = stim.im(:,end:-1:1,:);
  
  % save stim back to stimFile if called for
  if saveStimImage
    if isempty(saveStimImageTargets)
      saveStimImageTargets = stimfile;
    end
    saveStimImageToStimfile(stim,saveStimImageTargets);
  end
else
  % stim image was stored, just reclaim it
  disp(sprintf('(pRFGetStimImageFromStimfile) Loaded stim image from stimfile.'));
  stim = s.pRFStimImage;
end

% apply transformations if called for
if xFlip
  disp(sprintf('(pRFGetStimImageFromStimfile) X flipping stimulus image'));
  stim.im = stim.im(end:-1:1,:,:);
end
if yFlip
  disp(sprintf('(pRFGetStimImageFromStimfile) Y flipping stimulus image'));
  stim.im = stim.im(:,end:-1:1,:);
end
if timeShift
  disp(sprintf('(pRFGetStimImageFromStimfile) Time shifting stimulus image by %i',timeShift));
  stim.im = circshift(stim.im,[0 0 timeShift]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    closeStimulusDisplay    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeStimulusDisplay

try
  mglFrameGrab('end');
catch
end
try
  mglClose;
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    haveSameStimulus          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = haveSameStimulus(stimfile,volTrigRatio,recomputeStimImage)

tf = false;
if ~iscell(stimfile) || numel(stimfile) < 2
  return
end

% Put the per-component volume/trigger setting into the comparison.
if isempty(volTrigRatio)
  componentVolTrigRatio = repmat({[]},size(stimfile));
elseif iscell(volTrigRatio)
  if numel(volTrigRatio) == 1
    componentVolTrigRatio = repmat(volTrigRatio,size(stimfile));
  elseif numel(volTrigRatio) == numel(stimfile)
    componentVolTrigRatio = reshape(volTrigRatio,size(stimfile));
  else
    return
  end
elseif isscalar(volTrigRatio)
  componentVolTrigRatio = repmat({volTrigRatio},size(stimfile));
elseif numel(volTrigRatio) == numel(stimfile)
  componentVolTrigRatio = reshape(num2cell(volTrigRatio),size(stimfile));
else
  return
end

reference = [];
for iStimfile = 1:numel(stimfile)
  s = getStimfile(stimfile{iStimfile});
  if isempty(s) || iscell(s)
    return
  end

  % A saved image is the simplest exact case. Compare the values that are
  % returned to pRFFit directly; no hashes or secondary cache are needed.
  if isfield(s,'pRFStimImage') && ~recomputeStimImage
    requiredFields = {'im','x','y','t'};
    if ~isstruct(s.pRFStimImage) || ~all(isfield(s.pRFStimImage,requiredFields))
      return
    end
    signature.kind = 'savedImage';
    for iField = 1:numel(requiredFields)
      fieldName = requiredFields{iField};
      signature.(fieldName) = s.pRFStimImage.(fieldName);
    end
  else
    if isfield(s,'pRFStimImage')
      s = rmfield(s,'pRFStimImage');
    end
    [validStimfile s taskNum] = checkStimfile(s);
    if ~validStimfile || isempty(taskNum) || ...
        ~isfield(s,'stimulus') || ~isfield(s.stimulus,'stimulusType')
      return
    end

    s.myscreen = makeTraces(s.myscreen,0);
    e = getTaskParameters(s.myscreen,s.task{taskNum});
    if isempty(e) || ~isfield(e,'parameter') || ...
        ~isfield(e,'randVars') || ~isfield(e,'trialVolume')
      return
    end

    % Remove only end-of-run drawing state. The remaining stimulus structure
    % contains the actual geometry and stimulus settings.
    stimulus = s.stimulus;
    runtimeFields = {'phaseNum','phaseNumRect','last','currentMask', ...
      'cycleTime','blank','elementRotMatrix','maskBarRotMatrix'};
    stimulus = rmfield(stimulus,intersect(runtimeFields,fieldnames(stimulus)));

    signature.kind = 'parameters';
    signature.taskFilename = canonicalPRFTaskFilename(s.task{taskNum}{1}.taskFilename);
    signature.stimulus = stimulus;
    signature.parameter = e.parameter;
    signature.randVars = e.randVars;
    signature.trialVolume = e.trialVolume;
    signature.volTrigRatio = componentVolTrigRatio{iStimfile};

    screenFields = {'displayDistance','displaySize','imageWidth','imageHeight', ...
      'screenWidth','screenHeight','frametime'};
    for iField = 1:numel(screenFields)
      fieldName = screenFields{iField};
      if isfield(s.myscreen,fieldName)
        signature.screen.(fieldName) = s.myscreen.(fieldName);
      else
        signature.screen.(fieldName) = [];
      end
    end
  end

  if iStimfile == 1
    reference = signature;
  elseif ~isequaln(reference,signature)
    return
  end
  clear signature
end
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    canonicalPRFTaskFilename    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function taskFilename = canonicalPRFTaskFilename(taskFilename)

if ischar(taskFilename)
  taskFilename = lower(taskFilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    isSupportedPRFTaskFilename  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isSupportedPRFTaskFilename(taskFilename)

supportedTaskFilenames = {'mglretinotopy.m','gruretinotopy.m', ...
  'offsetretinotopy.m','mgldoublebars.m','mglmetalretinotopy.m'};
tf = ischar(taskFilename) && any(strcmp(canonicalPRFTaskFilename(taskFilename),supportedTaskFilenames));

%%%%%%%%%%%%%%%%%%%%%%%%%
%    createMaskImage    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function maskImage = createMaskImage(s,t,debugStimulus)

% find the beginning of the experiment
firstTimepoint = find(s.vol);
firstTimepoint = firstTimepoint(1);

% find the timepoint that corresponds to this time
thisTimepoint = s.time(firstTimepoint)+t;
thisTimepoint = find(thisTimepoint <= s.time);
if isempty(thisTimepoint)
  disp(sprintf('(pRFGetStimImageFromStimfile) Timepoint %0.1fs does not exist in stimfile. This might have happened if you have linked the wrong stimfile with the scan - in which case, setStimfile to change stimfile linking',t));
  maskImage = [];
  return
end
% get timepoint
thisTimepoint = thisTimepoint(1);

% get current volume number
thisVol = cumsum(s.vol);
thisVol = thisVol(thisTimepoint);

% get curent trial number
thisTrial = find(s.trialVol <= thisVol);
if isempty(thisTrial)
  disp(sprintf('(pRFGetStimImageFromStimfile) Trial for volume %i not found. The first trial starts at volume %i. Consider discarding the scan associated with stimfile %s',thisVol,first(s.trialVol),getLastDir(s.myscreen.stimfile)));
  maskImage = [];
  return;
end
thisTrial = thisTrial(end);

% make sure that the timepoint is valid for the trial
thisTimepoint = max(thisTimepoint,s.trialTicknum(thisTrial));
if length(s.trialTicknum) > thisTrial
  thisTimepoint = min(thisTimepoint,s.trialTicknum(thisTrial+1));
end

% pull out stimulus variable
global stimulus;
stimulus = s.stimulus;

if isfield(s,'barAngle')
  % now make a rotation matrix for the background angle
  elementAngle = s.elementAngle(thisTrial);
  co = cos(pi*elementAngle/180);
  si = sin(pi*elementAngle/180);
  stimulus.elementRotMatrix = [co si;-si co];

  % now make a rotation matrix for the bar angle we want to present
  barAngle = s.barAngle(thisTrial);
  co = cos(pi*barAngle/180);
  si = sin(pi*barAngle/180);
  stimulus.maskBarRotMatrix = [co si;-si co];

  % see whether this is a blank
  if barAngle == -1,blank = true;else blank = false;end

  % clear screen
  mglClearScreen(1);
else
  % for rings and wedges see if it is a blank
  blank = s.blank(thisTrial);
  mglClearScreen(0.5);
  mglFillOval(0,0,[stimulus.maxRadius*2 stimulus.maxRadius*2],[1 1 1]);
end

% set the current mask
stimulus.currentMask = s.maskPhase(thisTimepoint);
if stimulus.currentMask == 0,blank = true;end

% update bars only if this is not a blank frame
if ~blank
  updateRetinotopyStimulus2(stimulus,s.myscreen);
else
  % clear screen
  mglClearScreen(0.5);
end

mglFlush;
% grab the screen
maskImage = mglFrameGrab;

% make into a black and white image
% Treat the gray background as zero and all non-background pixels as the
% stimulus mask. This preserves the established reconstruction behavior
% when mglClearScreen returns gray rather than black.
% maskImage((maskImage > 0.51) | (maskImage < 0.49)) = 1;
% maskImage((maskImage < 0.51) & (maskImage > 0.49)) = 0;
maskImage((maskImage > 0.51) | (maskImage < 0.49)) = 0;
maskImage((maskImage < 0.51) & (maskImage > 0.49)) = 1;
maskImage = maskImage(:,:,1);

if blank
  maskImage = zeros(size(maskImage));
end

% Optional debug display. Keeping this off avoids a figure redraw for every
% reconstructed stimulus frame.
if debugStimulus
  disp(sprintf('Trial %i maskPhase: %i',thisTrial,stimulus.currentMask))
  mlrSmartfig('pRF','reuse');clf;imagesc(maskImage);
  title(sprintf('Trial %i maskPhase: %i',thisTrial,stimulus.currentMask))
  drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw retinotopy stimulus to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateRetinotopyStimulus2(stimulus,myscreen)

switch stimulus.stimulusType
  % if doing standard bars, then draw sliding elements
  case 3

    % draw background, 1 is sliding wedges
    if stimulus.backgroundType == 1
      % update the phase of the sliding wedges
      stimulus.phaseNumRect = 1+mod(stimulus.phaseNumRect,stimulus.nRect);

      % draw the whole stimulus pattern, rotate to the element angle
      x = stimulus.xRect{stimulus.phaseNumRect};
      y = stimulus.yRect{stimulus.phaseNumRect};
      coords(1:2,:) = stimulus.elementRotMatrix*[x(1,:);y(1,:)];
      coords(3:4,:) = stimulus.elementRotMatrix*[x(2,:);y(2,:)];
      coords(5:6,:) = stimulus.elementRotMatrix*[x(3,:);y(3,:)];
      coords(7:8,:) = stimulus.elementRotMatrix*[x(4,:);y(4,:)];
      %mglQuad(coords(1:2:8,:)+stimulus.xOffset,coords(2:2:8,:)+stimulus.yOffset,stimulus.cRect{stimulus.phaseNumRect},1);
    elseif stimulus.backgroundType == 2
      westheimer('do=update','stimulus',stimulus.westheimer,'frameNum',myscreen.tick);
    end

    % compute the center of the bar
    barCenter = repmat(stimulus.barCenter(stimulus.currentMask,:),size(stimulus.maskBarLeft,1),1);
    % compute the left and right masks (covering up everything except the bar)
    % by shifting by the barCenter and rotating the coordinates for the angle we want
    maskBarLeft = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarLeft)';
    maskBarRight = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarRight)';

    % draw the bar masks
    % Use black masks against the gray reconstruction background; the binary
    % conversion below treats non-background pixels as stimulus locations.
    % mglPolygon(maskBarLeft(1,:)+stimulus.xOffset,maskBarLeft(2,:)+stimulus.yOffset,myscreen.backgroundColor);
    % mglPolygon(maskBarRight(1,:)+stimulus.xOffset,maskBarRight(2,:)+stimulus.yOffset,myscreen.backgroundColor);
    mglPolygon(maskBarLeft(1,:)+stimulus.xOffset,maskBarLeft(2,:)+stimulus.yOffset,[0 0 0]);
    mglPolygon(maskBarRight(1,:)+stimulus.xOffset,maskBarRight(2,:)+stimulus.yOffset,[0 0 0]);
  case 4
    % doing dots task, so draw dots, first clear screen
    mglClearScreen;
    % draw the dots
    xOffset = stimulus.xOffset+stimulus.barCenter(stimulus.currentMask,1);
    yOffset = stimulus.yOffset+stimulus.barCenter(stimulus.currentMask,2);
    drawDots(stimulus.dots.middle,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.upper,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.lower,stimulus.maskBarRotMatrix,xOffset,yOffset);
    % update the dots
    stimulus.dots.upper = updateDots(stimulus.dots.upper);
    stimulus.dots.middle = updateDots(stimulus.dots.middle);
    stimulus.dots.lower = updateDots(stimulus.dots.lower);
    % draw the fixation cross
    global fixStimulus;
    mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);
    mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,stimulus.fixColor,fixStimulus.pos);
  case 5
    % doing dots task, so draw dots, first clear screen
    mglClearScreen;
    % draw the dots
    xOffset = stimulus.xOffset+stimulus.barCenter(stimulus.currentMask,1);
    yOffset = stimulus.yOffset+stimulus.barCenter(stimulus.currentMask,2);
    drawDots(stimulus.dots.middle,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.upper,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.lower,stimulus.maskBarRotMatrix,xOffset,yOffset);
    % update the dots
    stimulus.dots.upper = updateDots(stimulus.dots.upper);
    stimulus.dots.middle = updateDots(stimulus.dots.middle);
    stimulus.dots.lower = updateDots(stimulus.dots.lower);
  case {1,2}
  % update the phase of the sliding wedges
  stimulus.phaseNum = 1+mod(stimulus.phaseNum,stimulus.n);
  % draw the whole stimulus pattern
  %mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.c{stimulus.phaseNum},1);
  
  % mask out to get a wedge
  if stimulus.stimulusType == 1
    % For the mask is made of two overlapping semicircles.
    topX = stimulus.maskWedgeX{stimulus.currentMask}{1};
    topY = stimulus.maskWedgeY{stimulus.currentMask}{1};
    mglPolygon(topX+stimulus.xOffset,topY+stimulus.yOffset,0.5);

    bottomX = stimulus.maskWedgeX{stimulus.currentMask}{2};
    bottomY = stimulus.maskWedgeY{stimulus.currentMask}{2};
    mglPolygon(bottomX+stimulus.xOffset,bottomY+stimulus.yOffset,0.5);
    % or mask out to get a ring
  else
    mglPolygon(stimulus.maskInnerX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskInnerY{stimulus.currentMask}+stimulus.yOffset,0.5);
    mglQuad(stimulus.maskOuterX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskOuterY{stimulus.currentMask}+stimulus.yOffset,stimulus.maskOuterC{stimulus.currentMask});
  end
  otherwise
    error('pRFGetStimImageFromStimfile:UnsupportedStimulusType', ...
      'Unsupported stimulus type: %s',mat2str(stimulus.stimulusType));
end





%%%%%%%%%%%%%%%%%%%%%
%    getStimfile    %
%%%%%%%%%%%%%%%%%%%%%
function s = getStimfile(stimfile)

s = [];

% deal with a cell array of stimfiles (like in an average)
if iscell(stimfile)
  for i = 1:length(stimfile)
    s{i} = getStimfile(stimfile{i});
    if isempty(s{i}),return;end
  end
  return
end

% load stimfile
if isstr(stimfile)
  stimfile = setext(stimfile,'mat');
  if ~isfile(stimfile)
    disp(sprintf('(pRFGetStimImageFromStimfile) Could not open stimfile: %s',stimfile));
    return
  end
  s = load(stimfile);
elseif isstruct(stimfile)
  % see if this is a myscreen
  if isfield(stimfile,'imageWidth')
    % check for task field
    if isfield(stimfile,'task')
      s.task = stimfile.task;
      stimfile = rmfield(stimfile,'task');
    end
    % check for stimulus field
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
      stimfile = rmfield(stimfile,'stimulus');
    end
    % set myscreen field
    s.myscreen = stimfile;
    % else a variable with myscreen, task and stimulus or pRFStimImage
  elseif isfield(stimfile,'myscreen') || isfield(stimfile,'pRFStimImage')
    % copy fields over
    if isfield(stimfile,'myscreen')
      s.myscreen = stimfile.myscreen;
    end
    if isfield(stimfile,'task')
      s.task = stimfile.task;
    end
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
    end
    if isfield(stimfile,'pRFStimImage')
      s.pRFStimImage = stimfile.pRFStimImage;
    end
  end
end

% if you have a pRFStimImage then don't bother with the rest of the fields
if ~isfield(s,'pRFStimImage')
  % check fields
  checkFields = {'myscreen','task','stimulus'};
  for i = 1:length(checkFields)
    if ~isfield(s,checkFields{i})
      stimfileName = '';
      if isfield(s,'myscreen') && isfield(s.myscreen,'stimfile')
      	stimfileName = getLastDir(s.myscreen.stimfile);
      end
      disp(sprintf('(pRFGetStimImageFromStimfile) !!! Missing variable: %s in stimfile %s !!!',checkFields{i},stimfileName));
      s = [];
      return
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%
function [tf s taskNum] = checkStimfile(s)

tf = true;
s = cellArray(s);
taskNum = [];

stimulusType = [];
barAngle = [];
direction = [];

for i = 1:length(s)
  thiss = s{i};
  if isempty(thiss)
    disp(sprintf('(pRFGetStimImageFromStimfile) Missing stimfile'));
    tf = false;
    return
  end
  % if this has a pRFStimImage then we are ok
  if isfield(thiss,'pRFStimImage')
    continue;
  end
  dispstr = sprintf('%s: vols=%i',thiss.myscreen.stimfile,thiss.myscreen.volnum);
  if ~isfield(thiss,'stimulus') || ~isfield(thiss.stimulus,'stimulusType') || ...
      ~isnumeric(thiss.stimulus.stimulusType) || ...
      ~isscalar(thiss.stimulus.stimulusType) || ...
      ~ismember(thiss.stimulus.stimulusType,1:5)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp('(pRFGetStimImageFromStimfile:checkStimfile) Stimulus type must be one of the supported numeric values 1 through 5.');
    tf = false;
    return
  end
  % first check if this is a retinotpy stimfile - it should
  % have a task which is mglRetinotopy
  taskNum = [];
  for iTask = 1:2
    % Keep the bounds/type checks before dereferencing the task. Some valid
    % stimfiles contain only one task. Task filenames are matched without
    % case so mglDoublebars.m and mglDoubleBars.m are the same supported
    % program on case-sensitive filesystems.
    if (length(thiss.task) >= iTask) && iscell(thiss.task{iTask}) && ...
        ~isempty(thiss.task{iTask}) && isstruct(thiss.task{iTask}{1}) && ...
        isfield(thiss.task{iTask}{1},'taskFilename')
      fname = thiss.task{iTask}{1}.taskFilename;
      if isSupportedPRFTaskFilename(fname)
        taskNum = iTask;
      end
    end
  end
  if isempty(taskNum)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by mglRetinotopy'));
    tf = false;
    return
  end

  % check for proper saved fields
  missing = '';
  if ~isfield(thiss.task{taskNum}{1},'randVars') missing = 'randVars';end
  if ~isfield(thiss.task{taskNum}{1},'parameter') missing = 'parameter';end
  if ~any(strcmp('maskPhase',thiss.myscreen.traceNames)) missing = 'maskPhase';end
  if ~any(strcmp('blank',thiss.myscreen.traceNames)) && strcmp(canonicalPRFTaskFilename(thiss.task{taskNum}{1}.taskFilename),'mglretinotopy.m') missing = 'blank';end
  if ~isempty(missing)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the field %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one.',missing));
    tf = false;
    return
  end

  % check for necessary variables
  e = getTaskParameters(thiss.myscreen,thiss.task{taskNum}{1});

  % now check for each variable that we need
  varnames = {'blank'};
  if strcmp(canonicalPRFTaskFilename(thiss.task{taskNum}{1}.taskFilename),'mgldoublebars.m') % mglDoubleBars doesn't use blanks
    varnames = {};
  end
  for iVar = 1:length(varnames)
    varval = getVarFromParameters(varnames{iVar},e);
    if isempty(varval)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the variable %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one',varnames{iVar}));
      tf = false;
      return
    end
  end

  % check for matching stimfiles
  thisStimulusType = thiss.stimulus.stimulusType;
  if ~isempty(stimulusType) && ~isequal(stimulusType,thisStimulusType)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! Have you averaged together scans with different stimulus conditions?',dispstr));
  end
  stimulusType = thisStimulusType;
  if ~strcmp(canonicalPRFTaskFilename(thiss.task{taskNum}{1}.taskFilename),'mgldoublebars.m') && ismember(thisStimulusType,[3 4 5])
    varval = getVarFromParameters('barAngle',e);
    if ~isempty(barAngle) && ~isequal(varval,barAngle)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The barAngles are different! Have you averaged together scans with different stimulus conditions?',dispstr));
    end
    barAngle = varval;
  elseif strcmp(canonicalPRFTaskFilename(thiss.task{taskNum}{1}.taskFilename),'mgldoublebars.m')
    varval = getVarFromParameters('conditionNum', e);
    barAngle = thiss.stimulus.conditions(varval, 3:4);
  else
    if ~isempty(direction) && ~isequal(thiss.stimulus.direction,direction)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The directions are different! Have you averaged together scans with different stimulus conditions?',dispstr));
    end
    direction = thiss.stimulus.direction;
  end
end

s = s{end};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    saveStimImageToStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveStimImageToStimfile(stim,stimfile)

% make sure stimfile is a cell array
stimfile = cellArray(stimfile);

% first reload the stimfile
for iStimfile = 1:length(stimfile)
  if isfield(stimfile{iStimfile},'filename')
    s = load(stimfile{iStimfile}.filename);
    if isempty(s)
      disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Could not load stimfile %s. Unable to save stim image back to stimfile',stimfile{iStimfile}.filename));
    else
      % append the stim image and save back
      s.pRFStimImage = stim;
      save(stimfile{iStimfile}.filename,'-struct','s');
      disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Saved pRFStimImage to %s.',stimfile{iStimfile}.filename));
    end
  else
    disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Missing filename in stimfile structure, could not save pRFStimImage back to stimfile'));
  end
end
