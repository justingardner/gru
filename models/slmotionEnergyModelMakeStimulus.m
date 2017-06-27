% slmotionEnergyModelMakeStimulus.m
%
%        $Id$
%      usage: slmotionEnergyModelMakeStimulus
%         by: justin gardner
%             modified by steeve laquitaine
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Generates a stimulus for the motionEnergyMode, make a screen using mglEditScreenParams
%             called offscreen, with the desired resolution (typically, something smallish, like 100x100 pixels
%             setting the distance to 57 cm and the dimensions to 10 x 10 (will give you approximately 10 x 10 deg)
%             If you set the screenName to be 'offscreen', then the following calls will use those screen settings
%
%             [s msc] = slmotionEnergyModelMakeStimulus('screenName=offscreen');
%
%             % display stimulus
%             mlrVol(s{1}.s);
%
%             % compute motionEnergy filter response
%             slmotionEnergyModelMakeStimulus(s{1}.s,'myscreen',msc);
%
%             % run for multiple stimulus types
%             [s msc] = slmotionEnergyModelMakeStimulus('screenName=offscreen','coherence=[0:0.5:1]','direction=0:180:359','n=2');
%
function [stimulus myscreen] = slmotionEnergyModelMakeStimulus(varargin)

% parse args
getArgs(varargin,{'screenName',[],'coherence=1','direction=0','speed=2.8',...
    'density=16.7','stimulusLength=1','n=1'});

% use off screen context - which will display to a memory
% buffer so that we can just mglFrameGrab to get the images
mglSetParam('useCGL',0);
mglSetParam('offscreenContext',1);

% initalize the screen
myscreen = initScreen(screenName);

% task just has dots
task{1}.waitForBacktick = 0;
task{1}.seglen = [stimulusLength+0.5 0.5];
task{1}.numBlocks = n;
task{1}.parameter.dir = direction;
task{1}.parameter.coherence = coherence;

% number of frames to compute for
myscreen.nFrames = stimulusLength*myscreen.framesPerSecond;
myscreen.iFrame = 1;

% initialize our task
[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@updateScreenCallback);

% init the stimulus
global stimulus;
stimulus.dots.speed = speed;
stimulus.dots.density = density;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while phaseNum <= 1
  % update the dots
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
clear global stimulus;
mglSetParam('offscreenContext',0);

stimulus = myscreen.stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
if (task.thistrial.thisseg == 1)
  stimulus.dots.coherence = task.thistrial.coherence;
  myscreen.iFrame = 1;
  myscreen.stimulus{task.trialnum}.coherence = task.thistrial.coherence;
  myscreen.stimulus{task.trialnum}.dir = task.thistrial.dir;
  myscreen.stimulus{task.trialnum}.n = task.blocknum;
  % precompute data arrray for computing stimulus
  myscreen.stimulus{task.trialnum}.s = uint8(zeros(mglGetParam('screenWidth'),mglGetParam('screenHeight'),myscreen.nFrames));
  disp(sprintf('(motionEnergyModelMakeStimulus) Making stimulus with coherence: %0.2f direction: %0.1f',task.thistrial.coherence,task.thistrial.dir));
else
  stimulus.dots.coherence = 0;
end
stimulus.dots.dir = task.thistrial.dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
stimulus = updateDots(stimulus,myscreen);

% grab the frame, and average over color channels
if myscreen.iFrame <= myscreen.nFrames
  thisFrame = mglFrameGrab;
  myscreen.stimulus{task.trialnum}.s(:,:,myscreen.iFrame) = uint8(255*mean(thisFrame,3));

  % update frame counter
  myscreen.iFrame = myscreen.iFrame + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = myscreen.imageWidth;,end
% if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 2.5;,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 3;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 16.7;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 1;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = 2.8;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end

% actually a square patch of dots that get stenciled
% so calculate width and height
stimulus.dots.width = stimulus.dots.rmax*2;
stimulus.dots.height = stimulus.dots.rmax*2;

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;

% get initial position
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% get the dots step
stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;

% pick a random set of dots
stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
%stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;

% draw the dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
mglStencilSelect(0);




