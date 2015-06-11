% initStaircase.m
%
%        $Id:$ 
%      usage: doStaircase(command,varargin)
%         by: justin gardner
%       date: 07/27/10
%    purpose: wrapper function for various staircase methods, so that they can be used
%             interchageably in stimulus programs. (i.e. can be used to initialize
%             and run upDownStaircase and Quest staircases). The sequence of operations
%             is to first initialize the staircase that you want:
%
%             % e.g. for an 1 up 3 down upDown Staircase
%             s = doStaircase('init','upDown','nup=1','ndown=3','initialThreshold=0.5','initialStepsize=0.05','nTrials=40')
%             % or for a quest staircase
%             s = doStaircase('init','quest','tGuess',log10(0.5),'tGuessSd=2','nTrials=40');
%             % or for method of constant stimuli
%             s = doStaircase('init','fixed','fixedVals=0.3:0.1:0.7');
%
%             % then run in a loop getting and testing the current recommended test value
%             while ~doStaircase('stop',s)
%
%               [testValue s] = doStaircase('testValue',s);
%
%               % to update the staircase with the subject's response (1 or 0 for correct or incorrect)
%               s = doStaircase('update',s,r);
%
%             end
%
%             % then compute the final threshold
%             threshold = doStaircase('threshold',s);
%
%             You can get help for how to run by doing e.g.:
%             doStaircase('init','upDown','help');
%             doStaircase('init','quest','help');
%             doStaircase('init','fixed','help');
%             doStaircase('threshold','help');
%
%             If you want to update the staircase with a test-value that was not one recommended
%             by doing doStaircase('testValue'), then you add the test-value as the third argument on update
%             s = doStaircase('update',s,r,testValue);
% 
%             For a test run (with a simulated observer) testType can be quest/pest/levitt/sdt/ratings
%             doStaircase('test','testType=quest')
%             doStaircase('test','testType=pest')
%             doStaircase('test','testType=levitt')
%             doStaircase('test','testType=sdt','dprime=1','criterion=0')
%             doStaircase('test','testType=ratings','dprime=1','ratings=5')
%             These all take optional arguments
%             'nTrials=100' (number of trials)
%             'nStaircases=20' (number of staircase to run)
%
%
%             If you want to start a new staircase at the threshold determined by the last
%             staircase (with the same parameters), you can do:
%             s(end+1) = doStaircase('init',s(end));
%
%             You can also run a detection experiment, by using the sdt method. Responses in this case
%             are still correct or incorrect (i.e. *not* subject saying response present/absent)
%             doStaircase('init','sdt','help');
%             Note that in this case, doStaircase does not actually staircase the contrast values
%             This is something that could potentially be implemented, but is not currently
%
%             You can also run a ratings experiment. Responses in this case are *ratings*, that
%             is, if you have ratings set to 5, then the responses should be a value between
%             1 and 5.
%             doStaircase('init','ratings','ratings=5');
%
%             To combine the data from two staircases of the same type (only implemented for sdt and for
%             upDown). 
%             sCombined = doStaircase('combine',s1,s2);
%             You can also combine a cell or struct array of staircases
%             sCombined = doStaircase('combine',s)
%
function [retval retval2] = doStaircase(command,varargin)

retval = [];
retval2 = [];

% check arguments
if nargin == 0
  help doStaircase
  return
end

% check command
switch lower(command)
 case {'init','i','initstaircase'}
    retval = initStaircase(varargin);
 case {'update','u'}
    retval = updateStaircase(varargin);
 case {'threshold','t','getthreshold'}
    retval = getThreshold(varargin);
 case {'testvalue','v','gettestvalue'}
    [retval retval2] = getTestValue(varargin);
 case {'gethistory','h','history','hist'}
    retval = getHistory(varargin);
 case {'getpsycho','p','psycho','psychometric','psychometricfunction'}
    retval = getPsycho(varargin);
 case {'stop','s'}
    retval = getStop(varargin);
 case {'combine','c','combinestaircase','combinestaircases'}
    retval = combineStaircases(varargin);
 case {'test'}
    testStaircase(varargin);
 otherwise
  disp(sprintf('(doStaircase) Unknown command: %s',command));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function s = initStaircase(args)

s = [];
if length(args) < 1
  disp(sprintf('(doStaircase:initStaircase) Must specify type of staircase to init'));
  return
end

% split out type argument
type = args{1};
args = {args{2:end}};

% get some arguments for this function
help = [];
getArgs(args,{'help=0'},'suppressUnknownArgMessage=1');

% display help if necessary
if help
  disp(sprintf('initStaircase takes a type argument which can be either ''quest'' or ''upDown'''));
  disp(sprintf('To start a subsequent staircase based on the threshold found on this one:'));
  disp(sprintf('  s(end+1) = doStaircase(''init'',s(end));'));
  disp(sprintf('dispFig=0\t\t\tDisplay figure to show subject responses'));
  disp(sprintf('subplotRows=[]\t\t\tSet how many subplot rows to have'))
  disp(sprintf('subplotCols=[]\t\t\tSet how many subplot cols to have'))
  disp(sprintf('subplotNum=[]\t\t\tSet which subplot to dispay fig in'))
end

% check to see if we were passed in a staircase structure in which
% we restart the staircase based on that structure
if isStaircase(type)
  oldS = type;
  args = getRestartArgs(oldS,1);
  type = oldS(end).type;
else
  oldS = [];
end

% check arguments here
initArgs = args;dispFig = [];nTrials=[];subplotNum = [];subplotRows = [];subplotCols = [];
[argNames argValues args] = getArgs(args,{'dispFig=0','nTrials=[]','subplotNum=1','subplotRows=1','subplotCols=1'});

% check the type and farm off to specialied init functions
switch lower(type)
 case 'quest'
  s = initQuestStaircase(args);
 case 'updown'
  s = initUpDownStaircase(args);
 case 'fixed'
  s = initFixed(args);
 case 'sdt'
  s = initSdt(args);
 case 'ratings'
  s = initRatings(args);
 otherwise
  disp(sprintf('(doStaircase) Unknown staircase type: %s',type));
  return
end

if isempty(s),return,end

% setup default values
s.dispFig = dispFig;
s.subplotNum = subplotNum;
s.subplotRows = subplotRows;
s.subplotCols = subplotCols;
[tf s] = isStaircase(s);
s.initArgs = initArgs;

if ~isempty(nTrials)
  s.stopCriterionType = 'nTrials';
  s.stopCriterion = nTrials;
end

% setup figure
if dispFig
  smartfig('doStaircase','reuse');
  subplot(s(end).subplotRows,s(end).subplotCols,s(end).subplotNum);
  cla;
  title(sprintf('%s staircase',s.dispType));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    initUpDownStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = initUpDownStaircase(args)

s = [];

% get initial arguments
nup=[];
ndown=[];
initialThreshold=[];
initialStepsize=[];
minStepsize=[];
maxStepsize=[];
minThreshold=[];
maxThreshold=[];
verbose=[];
stepRule=[];
help=[];
getArgs(args,{'nup=1','ndown=3','initialThreshold=[]','initialStepsize=[]','minStepsize=[]','maxStepsize=[]','minThreshold=[]','maxThreshold=[]','verbose=1','stepRule=0','help=0'});

% display help
if help
  disp(sprintf('nup=%i\t\t\t\tHow many incorrect responses before we increase the current test value',nup));
  disp(sprintf('ndown=%i\t\t\t\tHow many correct responses before we decrease the current test value',ndown));
  disp(sprintf('initialThreshold\t\tThe start threshold to use'));
  disp(sprintf('initialStepsize\t\t\tthe starting stepsize by which test values are incremented by'));
  disp(sprintf('minStepsize\t\t\tthe minimum allowable stepsize - use in conjunction with levitt or pest rules'));
  disp(sprintf('maxStepsize\t\t\tthe maximum allowable stepsize - use in conjunction with the pest rule'));
  disp(sprintf('minThreshold\t\t\tThe minimum allowable threshold'));
  disp(sprintf('maxThreshold\t\t\tThe maximum allowable threshold'));
  disp(sprintf('verbose=%i\t\t\tDisplay text feedback',verbose));
  disp(sprintf('stepRule\t\t\tThe rule by which to change the stepsize. Should be either ''none'',''levitt'' or ''pest'' - see upDownStaircase for more info'));
  return
end
  
% check for mandatory arguments
mandatoryArguments = {'initialThreshold','initialStepsize'};
for i = 1:length(mandatoryArguments)
  if isempty(eval(mandatoryArguments{i}))
    disp(sprintf('(doStaircase:initUpDownStaircase) **** Must specify an %s ****',mandatoryArguments{i}));
    return
  end
end

% convert stepRule into string
switch lower(stepRule)
 case {0,'none'}
  stepRule = 'None';
 case {1,'levitt'}
  stepRule = 'Levitt';
 case {2,'pest'}
  stepRule = 'Pest';
end

% for pest make sure we have minThreshold
if strcmp(lower(stepRule),'pest')
  if isempty(minStepsize)
    if verbose
      disp(sprintf('(doStaircase) For PEST staircases it is best to specify a minimum stepsize'));
    end
  else
    initialStepsize(2) = minStepsize;
  end
  if isempty(maxStepsize)
    if verbose
      disp(sprintf('(doStaircase) For PEST staircases it is best to specify a maximum stepsize'));
    end
  else
    initialStepsize(3) = maxStepsize;
  end
end

% display what we are doing
if verbose
  disp(sprintf('(doStaircase) Initializing upDown staircase: nUp=%i nDown=%i initialThreshold=%f initialStepsize=%s stepRule=%s',nup,ndown,initialThreshold,num2str(initialStepsize),stepRule));
end

% ok, setup staircase
s.s = upDownStaircase(nup,ndown,initialThreshold,initialStepsize,stepRule);

% set min and max threshold
if ~isempty(minThreshold)
  s.s.minThreshold = minThreshold;
end

if ~isempty(maxThreshold)
  s.s.maxThreshold = maxThreshold;
end

s.type = 'upDown';
if ~strcmp(stepRule,'None')
  s.dispType = sprintf('%i-up %i-down %s',nup,ndown,stepRule);
else
  s.dispType = sprintf('%i-up %i-down',nup,ndown);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    initQuestStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = initQuestStaircase(args)

s = [];

% get initial arguments
initialThreshold=[];tGuess=[];
initialThresholdSd=[];tGuessSd=[];
pThreshold=[];
beta=[];
delta=[];
gamma=[];
grain=[];
range=[];
verbose=[];
stepRule=[];
help=[];
getArgs(args,{'initialThreshold=[]','initialThresholdSd=[]','tGuess=[]','tGuessSd=[]','pThreshold',1/sqrt(2),'beta=3.5','delta=0.01','gamma=0.5','verbose=1','logValues=[]','help=[]'});

if exist('QuestCreate') ~= 2
  disp(sprintf('(doStaircase:initQuestStaircase) You need to install the Quest package. See: http://www.psych.nyu.edu/pelli/software.html'));
  return
end

if help
  disp(sprintf('initialThreshold\t\tInitial threshold guess - this should be a linear (not log 10 based) value. You can speciy either initialThreshold or tGuess, but not both'));
  disp(sprintf('tGuess\t\t\t\tInitial threshold guess - this should be a log 10 based value.'));
  disp(sprintf('initialThresholdSd\t\tThe standard deviation of the prior. Should be a linear (not log 10 based value). You can speciy eith initialThresholdSd or tGuessSd but not both'));
  disp(sprintf('tGuessSd\t\t\tThe standard deviation of the prior - this should be a log 10 based value'));
  disp(sprintf('pThreshold=%f\t\tThe percent correct value at which threshold should be calculated',pThreshold));
  disp(sprintf('beta=%f\t\t\tSlope of the weibull function',beta));
  disp(sprintf('delta=%f\t\t\tPercent of trials where subject guesses',delta));
  disp(sprintf('gamma=%f\t\t\tPercent correct when subject is guessing - e.g. 0.5 for 2AFC',gamma));
  disp(sprintf('verbose=%i\t\t\tDisplay verbose information'));
  disp(sprintf('logValues\t\t\tIf set to 1, then all input and output will be log 10 based. Otherwise linear'));
  disp(sprintf('NOTE: Quest assumes that the psychometric function is weibull in linear coordinates *not log coordinates*. '));
  disp(sprintf('Do help on QuestCreate for more info on parameters'));
  return
end

% check for bad input values
if ~isempty(initialThreshold) && ~isempty(tGuess)
  disp(sprintf('(doStaircase:initQuestStaircase) You should only specify one of initialThreshold (linear threshold value) or tGuess (log base 10 threshold value)'));
  return
end

if ~isempty(initialThresholdSd) && ~isempty(tGuessSd)
  disp(sprintf('(doStaricase:initQuestStaircase) You should only specify one of initialThresholdSd (linear threshold prior std) or tGuessSd (log base 10 threshold prior std)'));
  return
end

% get thresholds
if ~isempty(initialThreshold),tGuess=log(initialThreshold);end
if ~isempty(initialThresholdSd),tGuessSd=abs(log(initialThresholdSd));end


% check for mandatory arguments
mandatoryArguments = {'tGuess','tGuessSd'};
for i = 1:length(mandatoryArguments)
  if isempty(eval(mandatoryArguments{i}))
    disp(sprintf('(doStaircase:initQuestStaircase) Must specify a %s',mandatoryArguments{i}));
    return
  end
end


% check to see whether we should return values as log or not
if isempty(logValues)
  if ~isempty(initialThreshold)
    logValues = 0;
    if verbose
      disp(sprintf('(doStaircase:initQuestStaircase) Since you have specified an initialThreshold as a linear value, then doStaircase will return testValues as linear values. If you want log 10 values instead, you can set the input variable logValues=1'));
    end
  else
    disp(sprintf('(doStaircase:initQuestStaircase) Since you have specified a tGuess as a log 10 value, then doStaircase will return testValues as log values. If you want linear values instead, you can set the input variable logValues=1'));
    logValues = 1;
  end
end

% display thresholds
if verbose
  disp(sprintf('(doStaircase:initQuestStaircase) tGuess=%f tGuessSd=%f (initialThreshold prior = %f<-%f->%f)',tGuess,tGuessSd,10^(tGuess-tGuessSd),10^(tGuess),10^(tGuess+tGuessSd)));
end

% display what we are doing
if verbose
  disp(sprintf('(doStaircase:initQuestStaircase) Initializing QUEST staircase: tGuess=%f tGuessSd=%f pThreshold=%f beta=%f delta=%f gamma=%f',tGuess,tGuessSd,pThreshold,beta,delta,gamma));
end

% ok, setup staircase
s.s = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
s.s.normalizePdf = 1;

s.type = 'Quest';
s.dispType = 'Quest';
s.logValues = logValues;

%%%%%%%%%%%%%%%%%%%
%    initFixed    %
%%%%%%%%%%%%%%%%%%%
function s = initFixed(args)

s = [];

fixedVals = [];verbose = [];help = [];
% get initial arguments
getArgs(args,{'fixedVals=[]','verbose=1','help=0'});

% display help
if help
  disp(sprintf('fixedVals\t\t\t\tSet of fixed values to test'));
  return
end
  
% check for mandatory arguments
mandatoryArguments = {'fixedVals'};
for i = 1:length(mandatoryArguments)
  if isempty(eval(mandatoryArguments{i}))
    disp(sprintf('(doStaircase:initFixed) **** Must specify %s ****',mandatoryArguments{i}));
    return
  end
end

s.s.fixedVals = fixedVals;
s.s.blockLen = length(fixedVals);
s.s.blockOrder = [];
s.s.blockTrialNum = s.s.blockLen+1;

s.type = 'Fixed';
s.dispType = 'Method of constant stimuli';
s.logValues = 0;

%%%%%%%%%%%%%%%%%%%%%
%    initRatings    %
%%%%%%%%%%%%%%%%%%%%%
function s = initRatings(args)

help = [];ratings = [];
[argNames argValues args] = getArgs(args,{'ratings=5','help=0'});

if help
  disp(sprintf('ratings=%i\tHave subject return ratings (if using ratings then set to the number of ratings the subject is asked to report, e.g. 5)',ratings));
  return
end

% initialize just like an sdt
s = initSdt(args);

% set the ratings type
s.type = 'Ratings';

% and keep the ratings 
s.s.ratings = ratings;


%%%%%%%%%%%%%%%%%
%    initSdt    %
%%%%%%%%%%%%%%%%%
function s = initSdt(args)

s = [];

strength = [];p = [];verbose = [];help = [];
% get initial arguments
getArgs(args,{'strength=[]','p=[]','verbose=1','help=0'});

% display help
if help
  disp(sprintf('strength\tSignal strength to test'));
  disp(sprintf('p=%f\t\tStimulus probability',p));
  return
end
  
% check for mandatory arguments
mandatoryArguments = {'strength'};
for i = 1:length(mandatoryArguments)
  if isempty(eval(mandatoryArguments{i}))
    disp(sprintf('(doStaircase:initSdt) **** Must specify %s ****',mandatoryArguments{i}));
    return
  end
end

% now if p was not specified make it so that there are an even number of absent and present
if isempty(p)
  p(1:length(strength)) = 0.5/(length(strength));
else
  %check p
  if sum(p) > 1
    disp(sprintf('(doStaircase:initSdt) **** Sum of p must be less than 1 ****'));
    keyboard
  elseif length(p) < length(strength)
    disp(sprintf('(doStaircase:initSdt) Must specify a p value for each signal strength'));
  end
end

s.s.strength = strength;
s.s.p = p;

s.type = 'Sdt';
s.dispType = 'Signal detection';
s.logValues = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%    updateStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function s = updateStaircase(args)


% check that we are passed in a staircase
s = [];
if (length(args) < 1) || ~isStaircase(args{1})
  disp(sprintf('(doStaircase:updateStaircase) ***** Must pass in staircase in to update ******'));
  keyboard
  return
end
if (length(args) < 2) 
  disp(sprintf('(doStaircase:updateStaircase) **** Must pass in a subject response (1 or 0 for true or false) to update ****'));
  keyboard
  return
end

% split out arguments
s = args{1};
r = args{2};
strength = [];
if length(args) > 2
  strength = args{3};
end

if isempty(r)
  disp(sprintf('(doStaircase) Response is empty. Should be 0 or 1'))
  return
end

% get what test value was recommended
if isempty(strength)
  if isempty(s(end).lastTestValue)
    disp(sprintf('(doStaircase:updateStaircase) Staircase does not have a lastTestValue - use should call getTestValue first before updating the staircase'));
    return
  end
  strength = s(end).lastTestValue;
end

switch(s(end).type)
 case 'Quest'
  [s(end) strength] = updateQuestStaircase(s(end),r,strength);
 case 'upDown'
  [s(end) strength] = updateUpDownStaircase(s(end),r,strength);
 case 'Fixed'
  [s(end) strength] = updateFixed(s(end),r,strength);
 case {'Sdt','Ratings'}
  [s(end) strength] = updateSdt(s(end),r,strength);
 otherwise
  disp(sprintf('(doStaircase:updateStaircase) Unknown staircase type %s',s(end).type));
end

% keep our own copy of the strength and response
s(end).testValues(end+1) = strength;
s(end).response(end+1) = r;
s(end).trialNum = s(end).trialNum+1;

% display the history
if (s(end).dispFig == 1)
  history = doStaircase('getHistory',s);
  smartfig('doStaircase','reuse');
  subplot(s(end).subplotRows,s(end).subplotCols,s(end).subplotNum);
  cla;
  % draw line
  plot(history.testValues,'k-');
  hold on
  xlabel('Trial #');
  ylabel('Test value');
  % for a normal non-ratings experiment
  if isempty(history.ratings)
    % plot the correct trials in green
    correctTrialNums = find(history.response==1);
    if ~isempty(correctTrialNums)
      plot(correctTrialNums,history.testValues(correctTrialNums),'go','MarkerFaceColor','g');
    end
    % plot the incorrect trials in red
    incorrectTrialNums = find(history.response==0);
    if ~isempty(incorrectTrialNums)
      plot(incorrectTrialNums,history.testValues(incorrectTrialNums),'ro','MarkerFaceColor','r');
    end
  else
    % get ratings info
    ratings = unique(history.ratings);
    maxRating = s.s.ratings;
    % size of markers to show rating
    minMarkerSize = 3;
    maxMarkerSize = 20;
    for rating = ratings
      % get the marker size to display at
      markerSize = minMarkerSize+(maxMarkerSize-minMarkerSize)*(rating-1)/(maxRating-1);
      % plot the correct trials in green
      correctTrialNums = find((history.response==1) & (history.ratings == rating));
      if ~isempty(correctTrialNums)
	plot(correctTrialNums,history.testValues(correctTrialNums),'go','MarkerFaceColor','g','MarkerSize',markerSize);
      end
      % plot the incorrect trials in red
      incorrectTrialNums = find((history.response==0) & (history.ratings == rating));
      if ~isempty(incorrectTrialNums)
	plot(incorrectTrialNums,history.testValues(incorrectTrialNums),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
      end
      
    end
  end
  drawnow
elseif s(end).dispFig == 2
  disp(sprintf('(doStaircase:update) %i:Strength: %f Response: %i',s(end).trialNum,strength,r));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    updateQuestStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s strength] = updateQuestStaircase(s,r,strength)

% convert to log, if necessary, and update quest
if (s.logValues==0)
  s.s = QuestUpdate(s.s,log10(strength),r);
else
  s.s = QuestUpdate(s.s,strength,r);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    updateUpDownStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s strength] = updateUpDownStaircase(s,r,strength)

s.s = upDownStaircase(s.s,r);

%%%%%%%%%%%%%%%%%%%%%
%    updateFixed    %
%%%%%%%%%%%%%%%%%%%%%
function [s strength] = updateFixed(s,r,strength)

s.s.blockTrialNum = s.s.blockTrialNum+1;


%%%%%%%%%%%%%%%%%%%
%    updateSdt    %
%%%%%%%%%%%%%%%%%%%
function [s strength] = updateSdt(s,r,strength)

%%%%%%%%%%%%%%%%%%%%%%
%    getTestValue    %
%%%%%%%%%%%%%%%%%%%%%%
function [testValue s] = getTestValue(args)

% check that we are passed in a staircase
s = [];testValue = [];
if (length(args) < 1) || ~isStaircase(args{1})
  disp(sprintf('(doStaircase:getTestValue) ******** Must pass in staircase in to update **********'));
  keyboard
  return
end

% check to make sure that we are being called with two output variables
if evalin('caller','nargout') ~= 2
  disp(sprintf('(doStaircase:getTestValue) Must be called with two outputs: [testValue s] = doStaircase(''testValue'',s);'));
  return
end

% split out arguments
s = args{1};
args = {args{2:end}};

switch(s(end).type)
 case 'Quest'
  testValue = getQuestTestValue(s(end),args);
 case 'upDown'
  testValue = getUpDownTestValue(s(end),args);
 case 'Fixed'
  [testValue s(end)] = getFixedTestValue(s(end),args);
 case {'Sdt','Ratings'}
  [testValue s(end)] = getSdtTestValue(s(end),args);
 otherwise
  disp(sprintf('(doStaircase:updateStaircase) Unknown staircase type %s',s(end).type));
end

s(end).lastTestValue = testValue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getQuestTestValue    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testValue = getQuestTestValue(s,args)

testValue = [];
type = 'Quantile';
help = 0;
getArgs(args,{'type','help'},'suppressUnknownArgMessage=1');

if help
  disp(sprintf('type=%s\t\tWhat method to get test value. Can be ''Quantile'',''Mean'' or ''Mode''',type));
  return
end

switch lower(type)
 case 'quantile'
  testValue = QuestQuantile(s.s);
 case 'mean'
  testValue = QuestMean(s.s);
 case 'mode'
  testValue = QuestMode(s.s);
 otherwise
  disp(sprintf('(doStaircase:getQuestTestValue) Unknown type: %s',type));
end

% convert to linear, if logValues is not set
if (s.logValues==0)
  testValue = 10^(testValue);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getUpDownTestValue    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testValue = getUpDownTestValue(s,args)

% just take the current value
testValue = s.s.threshold;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getFixedTestValue    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [testValue s] = getFixedTestValue(s,args)

% see if we have to restart block
if s.s.blockTrialNum > s.s.blockLen
  s.s.blockOrder(end+1,:) = randperm(s.s.blockLen);
  s.s.blockTrialNum = 1;
end

% get value
testValue = s.s.fixedVals(s.s.blockOrder(end,s.s.blockTrialNum));

%%%%%%%%%%%%%%%%%%%%%%%%%
%    getSdtTestValue    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [testValue s] = getSdtTestValue(s,args)

% the possible values of signal strength
strength = [s.s.strength 0];
% now choose a random value from the above array
% according to the probability desired
testValue = strength(first(find(rand<[cumsum(s.s.p) 1])));

%%%%%%%%%%%%%%%%%%%%%%
%    getThreshold    %
%%%%%%%%%%%%%%%%%%%%%%
function threshold = getThreshold(args)

threshold = [];

% check that we are passed in a staircase
s = [];
if (length(args) < 1)
  disp(sprintf('(doStaircase:getThreshold) Must pass in staircase in to update'));
  return
end

% split out arguments
s = args{1};
if isStaircase(s)
  args = {args{2:end}};
else
  s = [];
end
type = [];
[argNames argValues args] = getArgs(args,'type=[]');
help=[];
getArgs(args,'help=0','suppressUnknownArgMessage=1');
if help || ~isStaircase(s)
  disp(sprintf('type\t\tIf set to ''weibull'' does a weibull fit on all of the data instead of default'));
end

% default threshold type
if isempty(type) && isStaircase(s)
  type = s(end).type;
end

if isempty(type),return,end

threshold = [];averageThresholds = false;
switch(lower(type))
 case 'quest'
  % compute each individual threshold
  for i = 1:length(s)
    thisThreshold(i) = getQuestThreshold(s(i),args);
  end
  averageThresholds = true;  
 case 'updown'
  for i = 1:length(s)
    try
    thisThreshold(i) = getUpDownThreshold(s(i),args);
    catch
      disp('(doStaircase) Warning: DAN CODE... avoids failure here');
      thisThreshold(i) = NaN;
    end
  end
  averageThresholds = true;  
 case {'weibull','fixed'}
  threshold = getWeibullThreshold(s,args);
 case {'ratings'}
  threshold = getRatingsThreshold(s,args);
 case {'sdt'}
  threshold = getSdtThreshold(s,args);
 otherwise
  disp(sprintf('(doStaircase:updateStaircase) Unknown staircase type %s',s.type));
end

% for thresholds that are calculated per staircase (like a quest staircase)
% we keep all the individual data, and then compute a mean threshold across
% all staircases
if averageThresholds
  % get number of trials in each staircase
  for i = 1:length(s)
    h(i) = getHistory({s(i)}); 
  end
  n = [h.n];totalN = sum(n);
  % now concatenate all the fields of the thresholds together
  fields = fieldnames(thisThreshold);
  for i = 1:length(fields)
    % concatenate scalar values into an array and other values into a cell array
    if isscalar(thisThreshold(end).(fields{i}));
      threshold.(fields{i}) = [thisThreshold.(fields{i})];
    else
      threshold.(fields{i}) = {thisThreshold.(fields{i})};
    end
  end
  % keep each threshold
  threshold.eachThreshold = threshold.threshold;
  % calculate mean, weighted by number of trials
  threshold.threshold = sum((n/totalN).*threshold.eachThreshold);
  threshold.thresholdSTE = std(threshold.eachThreshold)/sqrt(length(threshold.eachThreshold));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getQuestThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getQuestThreshold(s,args)

method=[];
getArgs(args,{'method=default'});

threshold.threshold = QuestMean(s.s);
threshold.sd = QuestSd(s.s);

if s.logValues==0
  threshold.threshold = 10^(QuestMean(s.s));
  threshold.sd = 10^(QuestSd(s.s));
else
  threshold.threshold = QuestMean(s.s);
  threshold.sd = QuestSd(s.s);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getUpDownThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getUpDownThreshold(s,args)

threshold = [];
method=[];
getArgs(args,{'method=default','dispFig=0','maxIter=1000',...
  'dogoodnessoffit=0','dobootstrap=0'});

switch method
  case 'default'
   s.s = upDownStaircase(s.s,'dispFig',dispFig,'maxIter',maxIter,...
     'dogoodnessoffit',dogoodnessoffit,'dobootstrap',dobootstrap);
   threshold = s.s.computedThresholds;
   threshold.threshold = threshold.weibull;
 otherwise
  disp(sprintf('(doStaircase:getUpDownThreshold) Unknown threshold computation type %s',method));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getWeibullThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getWeibullThreshold(s,args)

p = [];help = [];gamma = [];dispFig=[];dispPsycho=[];doFit=[];
getArgs(args,{'p=[]','gamma=0.5','help=0','dispFig=0','dispPsycho=0','verbose=1','doFit=1'});

if help
  disp(sprintf('gamma=%f\t\tFits the weibull with the guessing rate (i.e. percent correct when no signal is present)',gamma));
  disp(sprintf('p\t\tComputes the threshold for the percent correct specified. Default returns the Threshold parameter of the weibull fit'));
  disp(sprintf('dispFig=%i\t\tDisplays figure with fit plot',dispFig));
  disp(sprintf('dispPsycho=%i\t\tDisplays just the psychometric function in the current axes',dispPsycho));
  disp(sprintf('doFit=%i\t\tIf set to 0, this will not fit the data, useful in combination with dispPsycho to just display psychometric function without a fit',dispPsycho));
  disp(sprintf('verbose=%i\t\tDisplay verbose info',verbose));
  disp(sprintf('NOTE that you can pass this routine an array or cell array of staircase structures and then it will compute the fit for all of the data at once'));
  threshold = [];
  return
end

testValues = [];
response = [];
staircaseEnd = [];
for iStaircase = 1:length(s)
  if iscell(s)
    h = doStaircase('getHistory',s{iStaircase});
  else
    h = doStaircase('getHistory',s(iStaircase));
  end
  testValues = [testValues h.testValues];
  response = [response h.response];
  staircaseEnd(end+1) = length(response);
end

% fit the weibull for linear coordinates


if doFit
  if (verbose) disppercent(-inf,sprintf('(doStaircase:getWeibullThreshold) Fitting weibull function to %i staircases (nTrials=%i)',length(s),length(testValues)));end  
  threshold.fit = fitweibull(testValues,response,[],gamma);
  if (verbose) disppercent(inf); end
else
  % put in bogus values for fit parameters
  threshold.fit.fitparams = nan(1,3);
end

% get the threshold value at the approriate p
if isempty(p)
  threshold.threshold = threshold.fit.fitparams(1);
  threshold.p = [];
else
  threshold.threshold = weibullInv(p,[threshold.fit.fitparams gamma]);
  threshold.p = p;
end

if (dispFig || dispPsycho)  && ~isempty(testValues)
  % get x range
  xMin = min(testValues);
  xMax = max(testValues);
  x = xMin:(xMax-xMin)/100:xMax;
  % get fit params
  fitParams = [threshold.fit.fitparams gamma];

  if strcmp(s(end).type,'Fixed')
    fixedVals = unique(testValues);
    for i = 1:length(fixedVals)
      %% 
      thisResponse = response(testValues==fixedVals(i));
      pCorrect(i) = sum(thisResponse)/length(thisResponse);
      pCorrectSte(i) = pCorrect(i)*(1-pCorrect(i))/sqrt(length(thisResponse));
    end
    binCenter = fixedVals;
    xDispMin = min(fixedVals);
    xDispMax = max(fixedVals);
  else
    % bin the responses to calculate % correct
    % first get a range inside the fit function of values to bin over
    nBins = min(6,length(unique(testValues)));

    % calculate by taking equal number of points
    sortedTestValues = sort(testValues);
    binLocs = min(round(1:length(sortedTestValues)/(nBins+1):length(sortedTestValues)),length(sortedTestValues));
    binLocs(end) = length(sortedTestValues);
    bins = sortedTestValues(binLocs);

    % no get the trials in each bin
    [counts indexes] = histc(testValues,bins);
    % and compute pCorrect and ste for each bin
    pCorrect = [];pCorrectSte = [];
    for i = 1:length(counts)-1
      thisIndexes = find(indexes==i);
      if length(thisIndexes) > 0
	pCorrect(i) = sum(response(thisIndexes))/length(thisIndexes);
	pCorrectSte(i) = pCorrect(i)*(1-pCorrect(i))/sqrt(length(thisIndexes)); 
      else
	pCorrect(i) = nan;
	pCorrectSte(i) = nan;
      end
      % get bin centers for plotting
      binCenter(i) = median(testValues(indexes==i));
    end
    % calculate bins according to placement on psychometric function
    minP = gamma+0.1;
    maxP = 1-fitParams(3)-0.1;
    pBins = minP:(maxP-minP)/(nBins-1):maxP;
    pSpacedBins = weibullInv(pBins,fitParams);
  
    % get range to display
    xDispMin = min(pSpacedBins)-1.5*(max(pSpacedBins)-min(pSpacedBins));
    xDispMax = max(pSpacedBins)+1.5*(max(pSpacedBins)-min(pSpacedBins));
  end

  threshold.fit.xdata = binCenter;
  threshold.fit.ydata = pCorrect;
  threshold.fit.ydataSte = pCorrectSte;
  
  if dispFig
    %now display the staircases
    f = smartfig('doStaircaseFit','reuse');
    clf(f)
    staircaseStart = 1;
    maxStaircaseLen = max(diff([0 staircaseEnd]));
    allStaircases = nan(length(staircaseEnd),maxStaircaseLen);
    % display each staircase
    for i = 1:length(staircaseEnd)
      % plot side by side
      subplot(2,3,1:3)
      plot(staircaseStart:staircaseEnd(i),testValues(staircaseStart:staircaseEnd(i)),'k-');
      hold on
      % plot on top of each other
      subplot(2,3,4)
      plot(testValues(staircaseStart:staircaseEnd(i)),'k-');hold on
      correctVals = response(staircaseStart:staircaseEnd(i));
      plot(find(correctVals),testValues(staircaseStart+find(correctVals)-1),'g.');
      plot(find(~correctVals),testValues(staircaseStart+find(~correctVals)-1),'r.');
      allStaircases(i,1:staircaseEnd(i)-staircaseStart+1) = testValues(staircaseStart:staircaseEnd(i));
      staircaseStart = staircaseEnd(i)+1;
    end

    % display correct and incorrect valuesnn
    subplot(2,3,1:3);
    correctVals = response;
    plot(find(correctVals),testValues(find(correctVals)),'g.');
    incorrectVals = ~response;
    plot(find(incorrectVals),testValues(find(incorrectVals)),'r.');
    xlabel('Trial #');
    ylabel('Test value');
    title(sprintf('%s staircases (n=%i) ',h.dispType,length(staircaseEnd)));
    vline(staircaseEnd);
    hline(threshold.threshold);
    subplot(2,3,4);
    plot(nanmean(allStaircases),'k-','LineWidth',2);
    xlabel('Trial #');
    ylabel('Test value');
    title('All staircases + average');
    hline(threshold.threshold);

    subplot(2,3,5);
  end

  % display psychometric function
  hline(threshold.threshold);
  cla;
  plot(x,weibull(x,fitParams),'k-','LineWidth',2);
  hold on
  myerrorbar(binCenter,pCorrect,'yError',pCorrectSte,'MarkerSize=4','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Symbol','o');
  plot(x,weibull(x,fitParams),'k-','LineWidth',2);
  
  xlabel('Signal strength');
  ylabel('Percent correct');
  legend('Weibull Fit','Binned data','Location','NorthWest');
  title(sprintf('Threshold (p=%f): %f (n=%i)\nT: %f Beta: %f delta: %f gamma: %f\n',p,threshold.threshold,length(testValues),fitParams(1),fitParams(2),fitParams(3),fitParams(4)));
  
  vline(threshold.threshold);
  if ~isempty(p)
    hline(p);
  end
  xaxis(xDispMin,xDispMax);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getRatingsThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = getRatingsThreshold(s,args)

% compute the sdt info for each ratings
for i = 1:s.s.ratings-1
  t.byRating(i) = getSdtThreshold(s,{'ratingsThreshold',i+1/2,'verbose=0'});
end

% calculate dprime as the mean dprime for now.
t.dprime = ([t.byRating.dprime]);
t.dprime = t.dprime(~isinf(t.dprime) & ~isnan(t.dprime));
t.dprime = mean(t.dprime);

% plot the ROC curve
plot([t.byRating.falseAlarmRate],[t.byRating.hitRate],'ko');hold on
xaxis(0,1);
xlabel('False Alarms');
yaxis(0,1);
ylabel('Hits');
dline;
axis square;

disp(sprintf('(doStaircase:getRatingsThreshold) mean dprime: %s (all dprimes %s)',mynum2str(t.dprime),mynum2str([t.byRating.dprime])));
  
%%%%%%%%%%%%%%%%%%%%%%%%%
%    getSdtThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function t = getSdtThreshold(s,args)

% this part is used so that getRatingsThreshold 
% can repeatedly call getSdtThreshold to compute
% hits/fa etc for each rating level
ratingsThreshold = [];verbose = [];
getArgs(args,{'ratingsThreshold=[]','verbose=1'});

% compute d'
h = doStaircase('history',s,'ratingsThreshold',ratingsThreshold);

% get unique signals
t.uniqueSignals = setdiff(unique(h.testValues),0);

for i = 1:length(t.uniqueSignals)
  % get signal present and signal absent trials
  signalPresent = h.testValues == t.uniqueSignals(i) ;
  signalAbsent = h.testValues == 0;
  t.nSignalPresent(i) = sum(signalPresent);
  t.nSignalAbsent(i) = sum(signalAbsent);
  
  % categorize
  t.hits(i) = sum(signalPresent & h.response);t.hitRate(i) = t.hits(i)/t.nSignalPresent(i);
  t.misses(i) = sum(signalPresent & ~h.response);t.missRate(i) = t.misses(i)/t.nSignalPresent(i);
  t.correctRejects(i) = sum(signalAbsent & h.response);t.correctRejectRate(i) = t.correctRejects(i)/t.nSignalAbsent(i);
  t.falseAlarms(i) = sum(signalAbsent & ~h.response);t.falseAlarmRate(i) = t.falseAlarms(i)/t.nSignalAbsent(i);

  % compute dprime and criterion
  t.zHits(i) = -sqrt(2)*erfcinv(2*t.hitRate(i));
  t.zFalseAlarms(i) = -sqrt(2)*erfcinv(2*t.falseAlarmRate(i));
  t.dprime(i) = t.zHits(i)-t.zFalseAlarms(i);
  t.criterion(i) = -0.5*(t.zHits(i)+t.zFalseAlarms(i));

  if verbose
    disp(sprintf('(doStaircase:getSdtThreshold) Signal strength %f: hits %0.1f%% (%i/%i) misses %0.1f%% (%i/%i) correct rejects %0.1f%% (%i/%i) false alarms %0.1f%% (%i/%i)',t.uniqueSignals(i),100*t.hitRate(i),t.hits(i),t.nSignalPresent(i),100*t.missRate(i),t.misses(i),t.nSignalPresent(i),100*t.correctRejectRate(i),t.correctRejects(i),t.nSignalAbsent(i),100*t.falseAlarmRate(i),t.falseAlarms(i),t.nSignalAbsent(i)));
    disp(sprintf('                              dprime=%f criterion=%f zHits=%f zFalseAlarms=%f',t.dprime(i),t.criterion(i),t.zHits(i),t.zFalseAlarms(i)));
  end
end

%%%%%%%%%%%%%%%%%%%
%    getPsycho    %
%%%%%%%%%%%%%%%%%%%
function p = getPsycho(args)

p = [];

% display help if necessary
if isstr(args{1})
  disp(sprintf('Returns a structure with the pyschometric function. Useful for when you have used a fixed (method of constant stimuli)'));
  return
end

if (length(args) < 1) || ~isStaircase(args{1})
  disp(sprintf('(doStaircase:getPsycho) Must pass in staircase'));
  return
end

s = args{1};
%args = {args{2:end}};
%getArgs(args,{'linearValues'});

% make sure this is fixed
if ~strcmp(s.type,'Fixed')
  disp(sprintf('(doStaircase:getPsycho) Psychometric function only available for type fixed (Method of constant stimuli). This staircase is ''%s''.',s.type));
  return
end

% get the history
h = doStaircase('getHistory',s);

% get the test values
p.stimStrength = unique(h.testValues);

% for each test value, compute percent correct
for i = 1:length(p.stimStrength)
  whichTrials = find(h.testValues==p.stimStrength(i));
  p.n(i) = length(whichTrials);
  p.pCorrect(i) = sum(h.response(whichTrials))/p.n(i);
end

%%%%%%%%%%%%%%%%%%%%
%    getHistory    %
%%%%%%%%%%%%%%%%%%%%
function h = getHistory(args)


h = [];

% display help if necessary
if isstr(args{1})
  disp(sprintf('Returns a structure with the staircase history. Includes fields'));
  disp(sprintf('testValues: the values chosen to test at'));
  disp(sprintf('response: the responses (correct or incorrect) of the subject'));
  disp(sprintf('n: Number of trials run'));
  disp(sprintf('dispType: A string that says what type of staircase the strucutre is'));
  return
end

if (length(args) < 1) || ~isStaircase(args{1})
  disp(sprintf('(doStaircase:getHistory) Must pass in staircase'));
  return
end

s = args{1};
args = {args{2:end}};
linearValues = [];ratingsThreshold = [];
getArgs(args,{'linearValues=1','ratingsThreshold=[]'});

h.testValues = s(end).testValues;
h.response = s(end).response;
h.n = length(s(end).testValues);
h.dispType = s(end).dispType;

if s(end).logValues && linearValues
  h.testValues = 10.^h.testValues;
end

% see if this is a ratings experiment in which the response values 
% are treated as a rating and compared against a ratingsThreshold
if strcmp(lower(s(end).type),'ratings')
  % default to setting correct incorrect based on the middle of the ratings scale
  if isempty(ratingsThreshold), ratingsThreshold = s(end).s.ratings/2;end
  % now fill in response
  h.ratings = h.response;
  % correct trials are either hits or correct rejects
  h.response = ((h.response > ratingsThreshold) & (h.testValues>0)) | ((h.response <= ratingsThreshold) & (h.testValues==0));
else
  h.ratings = [];
end

  
%%%%%%%%%%%%%%%%%
%    getStop    %
%%%%%%%%%%%%%%%%%
function retval = getStop(args)

retval = [];

% display help if necessary
if isstr(args{1})
  disp(sprintf('Returns true/false of whether it is recommended to stop the staircase'));
  disp(sprintf('This is usually just done by checking to see if you have run more than'));
  disp(sprintf('the nTrials set in initStaircase, but in the future could use more sophisticated stop criterion'));
  return
end

s = [];
if length(args) < 1
  disp(sprintf('(doStaircase:getStop) Must pass staircase in'));
  return
end

% split out type argument
s = args{1};
if ~isStaircase(s(end))
  disp(sprintf('(doStaircase:getStop) Passed in argument is not a staircase'));
  return
end

% check stop criterion
if strcmp(s(end).stopCriterionType,'nTrials')
  if s(end).trialNum >= s(end).stopCriterion
    retval = 1;
  else
    retval = 0;
  end
end
  
%%%%%%%%%%%%%%%%%%%%%
%    isStaircase    %
%%%%%%%%%%%%%%%%%%%%%
function [tf s] = isStaircase(s)

tf = 0;
requiredFields = {'s','type'};
optionalFields = {'logValues',0;
		  'lastTestValue',[];
		  'testValues',[];
		  'response',[];
		  'dispType',[];
		  'initArgs',[];
		  'trialNum',0;
		  'stopCriterionType','nTrials';
		  'stopCriterion',50;
		  'dispFig',0};

% only one argument output then all fields are required
if nargout == 1
  requiredFields = {optionalFields{:,1} requiredFields{:}};
  optionalFields = {};
end

if isempty(s),return,end

% check required fields
for i = 1:length(requiredFields)
  if ~isfield(s(end),requiredFields{i});
    return
  end
end

% setup/check optional fields
tf = 1;
for i = 1:size(optionalFields,1)
  if ~isfield(s,optionalFields{i,1})
    s.(optionalFields{i,1}) = optionalFields{i,2};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    getRestartArgs    %
%%%%%%%%%%%%%%%%%%%%%%%%
function initArgs = getRestartArgs(s,verbose)

if ~isequal(s(end).type,'Fixed')
  % return arguments to restart with the initial threshold set to the
  % calculated threshold
  threshold = getThreshold({s});

  % remove initialThreshold and 
  [argNames argVals initArgs] = getArgs(s(end).initArgs,{'initialThreshold=0','tGuess=0'},'suppressUnknownArgMessage=1');

  % now set the initial threshold in the args and return
  initArgs = {initArgs{:} 'initialThreshold' threshold.threshold};

  % display restart threshold
  if verbose
    disp(sprintf('(doStaircase) Setting restart parameters with initialThreshold=%f based on computed threshold',threshold.threshold));
  end
else
  initArgs = s(end).initArgs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    combineStaircases    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = combineStaircases(args)

if length(args) == 1
  % passed in a cell array
  if iscell(args{1})
    % check if each is a stair
    if ~isStaircase(args{1}{1})
      disp(sprintf('(doStaircase:combineStaircases) !!! Must pass in a cell array of staircases !!!!'));
      return
    end
    % passed test, so combine each one
    s = args{1}{1};
    for iStair = 2:numel(args{1})
      s = doStaircase('combine',s,args{1}{iStair});
    end
    return
  end
  % if a struct, do same thing, but struct referncing
  for iStair = 1:length(args{1})
    if ~isStaircase(args{1}(iStair))
      disp(sprintf('(doStaircase:combineStaircases) !!! Must pass in an array of staircases !!!!'));
      return
    end
  end
  % found an array of staircases
  s = args{1}(1);
  for iStair = 2:length(args{1})
    s = doStaircase('combine',s,args{1}(iStair));
  end
  return
end


s = [];
if (length(args) < 2) || ~isStaircase(args{1}) || ~isStaircase(args{2})
  disp(sprintf('(doStaircase:combineStaircases) ***** Must pass in the two staircases you want to combine ******'));
  return
end
s1 = args{1};
s2 = args{2};

% now make sure they are the same kind of staircase
if ~strcmp(s1.type,s2.type)
  disp(sprintf('(doStaircase:combineStaircases) Staircases have mismatched type: %s vs %s',s1.type,s2.type));
  return
end

if strcmp(lower(s1.type),'sdt')
  % compare settings
  if ~isequal(s1.s,s2.s)
    disp(sprintf('(doStaircase:combineStaircases) Staircases have mismatched settings'));
    disp(sprintf('s1.s: '));
    s1.s
    disp(sprintf('s2.s: '));
    s2.s
    if ~askuser('Continue anyway'),return,end
  end

  % ok, now we are ready to combine the two
  s1.testValues = [s1.testValues s2.testValues];
  s1.response = [s1.response s2.response];
  s1.trialNum = s1.trialNum + s2.trialNum;
  s1.lastTestValue = s2.lastTestValue;

  s = s1;
elseif strcmp(lower(s1.type),'updown')
  % check to make sure fields match
  checkFields = {'type','upn','downn','stepsizeRule'};
  for iCheck = 1:length(checkFields)
    if ~isequal(s1.s.(checkFields{iCheck}),s2.s.(checkFields{iCheck}))
      disp(sprintf('(doStaircase:combineStaircases) Staircase have mismatched settings for: %s',checkFields{iCheck}));
      if ~askuser('Continue anyway'),return,end
    end
  end
  % conact over the upDownStaircase fields
  concatFields = {'response','group','reversal','strength'};
  for iConcat = 1:length(concatFields)
    s1.s.(concatFields{iConcat}) = [s1.s.(concatFields{iConcat}) s2.s.(concatFields{iConcat})];
  end
  % copy fields
  copyFields = {'direction','threshold','stepsize'};
  for iCopy = 1:length(copyFields)
    s1.s.(copyFields{iCopy}) = s2.s.(copyFields{iCopy});
  end
  % concat + add num trials
  concatAddFields = {'reversals'};
  for iConcatAdd = 1:length(concatAddFields)
    % check for existence of field in both structures
    if isfield(s1.s,concatAddFields{iConcatAdd}) && isfield(s2.s,concatAddFields{iConcatAdd}) 
      s1.s.(concatAddFields{iConcatAdd}) = [s1.s.(concatAddFields{iConcatAdd}) s1.s.n+s2.s.(concatAddFields{iConcatAdd})];
    % if only in second structure, just bring that one over then
    elseif isfield(s2.s,concatAddFields{iConcatAdd}) 
      s1.s.(concatAddFields{iConcatAdd}) = s1.s.n+s2.s.(concatAddFields{iConcatAdd});
    end
  end
  % sum fields
  sumFields = {'n','reversaln'};
  for iSum = 1:length(sumFields)
    s1.s.(sumFields{iSum}) = s1.s.(sumFields{iSum})+s2.s.(sumFields{iSum});
  end

  % combined fields for doStaircase structure
  s1.testValues = [s1.testValues s2.testValues];
  s1.response = [s1.response s2.response];
  s1.trialNum = s1.trialNum + s2.trialNum;
  s1.lastTestValue = s2.lastTestValue;
  s = s1;
else
  disp(sprintf('(doStaircase:combineStaircases) Combining staircases of type %s not implemented yet',s1.type));
  return
end  

%%%%%%%%%%%%%%%%%%
%    dispCell    %
%%%%%%%%%%%%%%%%%%
function dispCell(x,sigfigs)

if ieNotDefined('sigfigs'), sigfigs = 2;end
for i = 1:length(x)
  if isstr(x{i})
    mydisp(x{i})
  elseif isnumeric(x{i})
    mydisp(sprintf('%s',mynum2str(x{i},'sigfigs',sigfigs,'compact=1')));
  else
    disp(x)
  end
  mydisp(' ');
end


%%%%%%%%%%%%%%%%%%%%%%%
%    testStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function testStaircase(args)

testType = [];
[argNames argValues args] = getArgs(args,{'testType=quest','nTrials=100','nStaircases=20'});
disp(sprintf('(doStaircase) Testing %s with doStaircase for %i staircases of %i trials',testType,nStaircases,nTrials));

% test sdt
if strcmp(testType,'sdt')
  testSdt(args,nTrials,nStaircases);
  return
elseif strcmp(testType,'ratings')
  testRatings(args,nTrials,nStaircases);
  return
end
  
% threshold of observer and parameters of their
% weibull function.
threshold = 0.4;
observerParams = [log10(threshold) 3.5 0.01 0.5];

% set the initial threshold
initialThreshold = threshold*rand;
initialStepsize = 0.1;
allTestValues = [];

% init figure
smartfig('doStaircase','reuse');clf;

for iStaircase = 1:nStaircases

  % display simulated observer psychometric function
  x = 0.1:0.01:1;
  smartfig('doStaircase','reuse');
  subplot(2,3,4);cla;
  semilogx(x,weibull(log10(x),observerParams,10));
  vline(threshold);
  xlabel('Contrast (log axis)');
  ylabel('Simulated observer percent correct');
  title(sprintf('Simulated observer psychometric function\nThreshold: %f',threshold));
  subplot(2,3,1:3);cla;

  % initialize the staircase
  if iStaircase == 1
    switch lower(testType)
     case {'quest','q'}
      % initialize a quest staircase
      s = doStaircase('init','quest','dispFig=1','initialThreshold',initialThreshold,'tGuessSd',2,'pThreshold=0.75','nTrials',nTrials);
     case {'levitt','l'}
      % initialize an up Down staircase with levitt rule
      s = doStaircase('init','upDown','dispFig=1','initialThreshold',initialThreshold,'initialStepsize',initialStepsize,'stepRule','Levitt','nTrials',nTrials);
     case {'pest','p'}
      % initialize an up Down staircase with pest rule
      s = doStaircase('init','upDown','dispFig=1','initialThreshold',initialThreshold,'initialStepsize',initialStepsize,'stepRule','Pest','minStepsize=0.01','maxStepsize=0.2','nTrials',nTrials);
     case {'constant','c','methodofconstantstimuli','fixed'}
      % initialize an up Down staircase with pest rule
      s = doStaircase('init','fixed','fixedVals',0:0.1:0.7,'nTrials',nTrials,'dispFig=1');
     otherwise
      disp(sprintf('(doStaircase:testStaircase) Unknown testType: %s',testType));
    end
  else
    % restart based on current computed threshold from last staircase
    s(end+1) = doStaircase('init',s(end));
  end

  % simulate trials
  while ~doStaircase('stop',s)
    % get the current test value
    [testValue s]= doStaircase('testValue',s);

    % simulate observer response as a weibull
    randDraw = rand;
    pObserverCorrect = weibull(log10(testValue),observerParams,10);
    response = randDraw < pObserverCorrect;

    % update the staircase
    s = doStaircase('update',s,response);
  end

  % get threshold
  t(iStaircase) = doStaircase('threshold',s(iStaircase));

  % get history
  h = doStaircase('getHistory',s);

  % plot bar graph of how many values were used along psychometric function
  subplot(2,3,5);cla
  allTestValues = [allTestValues h.testValues];
  [y x] = hist(allTestValues,min(10*iStaircase,100));
  bar(x,y/max(y));
  hold on
  plot(x,2*(weibull(log10(x),observerParams,10)-0.5));
  xlabel('Contrast - linear axis');
  ylabel('Percent of trials');
  title('Histogram of test values');

  % display threshold histogram
  if iStaircase > 1
    subplot(2,3,6);cla
    hist([t.threshold]);
    xlabel('Computed threshold');
    ylabel('Number of staircases');
    title(sprintf('Weibull threshold per staircase\nMean %f Std %f',mean([t.threshold]),std([t.threshold])));
  end
end

% fit to all of data
doStaircase('t',s,'type=weibull','dispFig=1','p',1/sqrt(2));

%%%%%%%%%%%%%%%%%
%    testSdt    %
%%%%%%%%%%%%%%%%%
function testSdt(args,nTrials,nRuns)

dprime = [];criterion = [];
getArgs(args,{'dprime=1','criterion=0'});

% convert criterion so that 0.5 means in between the two distributions
criterion = (criterion+1)/2;

% given the dprime, set the observerStd assuming signal strength of 1
observerStd = 1/dprime;

smartfig('doStaircase','reuse');clf;

for iRun = 1:nRuns
  subplot(2,2,1:2);cla
  % init the sdt
  s = doStaircase('init','sdt','dispFig=1','strength=1','p=0.5','nTrials',nTrials);

  % simulate trials
  while ~doStaircase('stop',s)
    % get the current test value
    [testValue s]= doStaircase('testValue',s);
    
    % simulate observer response
    observerResponse = randgauss(testValue,observerStd,1);
    if observerResponse>criterion
      % subject responds target present
      if testValue == 0,response = false;else response = true;end
    else
      % subject responds target absent
      if testValue == 0,response = true;else response = false;end
    end

    % update the staircase
    s = doStaircase('update',s,response);
  end

  % compute dprime
  t = doStaircase('threshold',s);
  
  % get dprime and criterion
  simulatedDprime(iRun) = t.dprime;
  simulatedCriterion(iRun) = t.criterion;
  
  % plot as distribution
  subplot(2,2,3);cla;
  thisDist = simulatedDprime(~isinf(simulatedDprime));
  if length(thisDist)>1,myhist(thisDist);end
  xlabel('dprime');
  ylabel('n runs');

  subplot(2,2,4);cla;
  thisDist = simulatedCriterion(~isinf(simulatedCriterion));
  if length(thisDist)>1,myhist(thisDist);end
  xlabel('criterion');
  ylabel('n runs');
end


%%%%%%%%%%%%%%%%%%%%%
%    testRatings    %
%%%%%%%%%%%%%%%%%%%%%
function testRatings(args,nTrials,nRuns)

dprime = [];criterion = [];ratings = [];
getArgs(args,{'dprime=1','ratings=5'});

% convert criterion so that 0.5 means in between the two distributions
criterion = (criterion+1)/2;

% given the dprime, set the observerStd assuming signal strength of 1
observerStd = 1/dprime;

% compute what values the simulated observer will use
% to catergoize into what rating
ratingMag = [-inf -1:3/(ratings-2):2 inf];

% set up figure
smartfig('doStaircase','reuse');clf;

for iRun = 1:nRuns
  subplot(2,2,1:2);cla;
  % init the sdt
  s = doStaircase('init','ratings','ratings',ratings,'dispFig=1','strength=1','p=0.5','nTrials',nTrials);

  % simulate trials
  while ~doStaircase('stop',s)
    % get the current test value
    [testValue s]= doStaircase('testValue',s);
    
    % simulate observer response
    observerResponse = randgauss(testValue,observerStd,1);
    response = max(find(observerResponse>ratingMag));
    
    % update the staircase
    s = doStaircase('update',s,response);
  end

  % compute dprime
  subplot(2,2,3);
  t = doStaircase('threshold',s);
  
  % get dprime and criterion
  simulatedDprime(iRun) = t.dprime;
  
  % plot as distribution
  subplot(2,2,4);cla;
  thisDist = denan(simulatedDprime(~isinf(simulatedDprime)));
  if length(thisDist)>1,myhist(thisDist);end
  xlabel('dprime');
  ylabel('n runs');
end





