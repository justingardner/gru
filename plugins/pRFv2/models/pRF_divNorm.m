% pRF_divNorm.m
%
%        $Id:$ 
%      usage: pRF_divNorm(varargin)
%         by: austin kuo
%       date: 04/02/25
%    purpose: divisive normalization model for pRFFit (based on Aqil, Knapen, and Dumoulin, 2021)
%             
%             - Simply add the model type to pRFGUI and create a 
%               file like this one for your model with the appropriate 
%               methods filled in.
% 
%             See pRF_gaussian for general documentation

function output = pRF_divNorm(varargin)

if nargin < 2
  disp(sprintf('Not enough arguments'));
  disp(sprintf('Number of arguments: %d', nargin));
  celldisp(varargin)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pRF_gaussian('getModelResponse', fitParams, rfModel, hrf, p, i) %%%%
%%%                  Called from getModelResidual                   %%%%
if strcmp(varargin{1}, 'getModelResponse')

  fitParams = varargin{2};
  rfModel1 = varargin{3};
  hrf = varargin{4};
  p = varargin{5};
  i = varargin{6};

  rfModel2 = exp(-(((fitParams.stimX - p.x).^2)/(2*(p.stdSurround^2))+((fitParams.stimY-p.y).^2)/(2*(p.stdSurround^2)))); 
  nFrames = fitParams.concatInfo.runTransition(i,2);
  rPlus = convolveModelWithStimulus(rfModel1,fitParams.stim{i},nFrames);
  rMinus = convolveModelWithStimulus(rfModel2, fitParams.stim{i}, nFrames);

  % and convolve in time.
  pPlus = convolveModelResponseWithHRF(rPlus,hrf);
  pMinus = convolveModelResponseWithHRF(rMinus,hrf);
  
  % Model response
  thisModelResponse = (p.gainCenter*pPlus + p.actConst) ./ (p.gainSurround*pMinus + p.divConst) - p.actConst/p.divConst; 

  % drop junk frames here
  thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);

  % return the calculated model response
  output = thisModelResponse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pRF_gaussian('setParams', fitParams) %%%%
%%%        Called from setParams         %%%%

elseif strcmp(varargin{1}, 'setParams')

  fitParams = varargin{2};

  fitParams.paramNames = {'x','y','centerWidth', 'surroundWidth', 'centerAmplitude', 'surroundAmplitude', 'actConst', 'divConst'};
  fitParams.paramDescriptions = {'RF x position','RF y position','RF width of pRF center gaussian (std)', 'Width of surround gaussian (std)', 'Beta weight amplitude of center Gaussian', 'Beta weight amplitude of surround Gaussian', 'Activation constant', 'Divisive constant'};
  fitParams.paramIncDec = [1 1 1 1 1 1 1 1];
  % fitParams.paramMin = [-inf -inf 0 0 -inf -inf -inf 0];
  % fitParams.paramMax = [inf inf inf inf inf inf inf inf];
  % set min/max and init
  fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0 0 0 0];
  fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf inf inf inf inf inf];
  fitParams.initParams = [0 0 4 4 1 1 10 100];

  % return fitParams with modified values
  output = fitParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  pRF_gaussian(command, fitParams, params)  %%%%%
%%%%         Called from getFitParams           %%%%% 

elseif strcmp(varargin{1}, 'getFitParams')

  fitParams = varargin{2};
  params = varargin{3};

  p.rfType = fitParams.rfType;

  % Define your parameters here
  p.x = params(1);
  p.y = params(2);
  p.std = params(3); % std is stdCenter
  p.stdSurround = params(4);
  p.gainCenter = params(5);
  p.gainSurround = params(6);
  p.actConst = params(7);
  p.divConst = params(8);
  % use a fixed single gaussian
  p.canonical.type = 'gamma';
  p.canonical.lengthInSeconds = 25;
  p.canonical.timelag = fitParams.timelag;
  p.canonical.tau = fitParams.tau;
  p.canonical.exponent = fitParams.exponent;
  p.canonical.offset = 0;
  p.canonical.diffOfGamma = fitParams.diffOfGamma;
  p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
  p.canonical.timelag2 = fitParams.timelag2;
  p.canonical.tau2 = fitParams.tau2;
  p.canonical.exponent2 = fitParams.exponent2;
  p.canonical.offset2 = 0;

  output = struct(p);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
function modelResponse = convolveModelWithStimulus(rfModel,stim,nFrames)

% get number of frames
nStimFrames = size(stim.im,3);

% preallocate memory
modelResponse = zeros(1,nStimFrames);

for frameNum = 1:nStimFrames
  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
  modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

