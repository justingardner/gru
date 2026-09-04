% pRF_divNorm.m
%
%        $Id:$ 
%      usage: pRF_divNorm(varargin)
%         by: Austin Kuo with Codex
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
%%% pRF_divNorm('getModelResponse', fitParams, rfModel, hrf, p, i) %%%%
%%%                  Called from getModelResidual                   %%%%
if strcmp(varargin{1}, 'getModelResponse')

  fitParams = varargin{2};
  rfModel1 = varargin{3};
  hrf = varargin{4};
  p = varargin{5};
  i = varargin{6};

  validateDivNormParameters(p);

  rfModel2 = pRFMakeGaussian(fitParams,p.x,p.y,p.stdSurround);
  rPlus = pRFProjectStimulus(rfModel1,fitParams.stim{i});
  rMinus = pRFProjectStimulus(rfModel2,fitParams.stim{i});

  % Apply divisive normalization to the neural response, then convolve the
  % completed response with the HRF. Scale all gains by the divisive constant
  % first. This is algebraically equivalent to (a*E+b)/(c*S+d)-b/d, returns
  % exact zero on blanks, and avoids overflow when a,b,c,d share a large scale.
  centerGain = p.gainCenter/p.divConst;
  surroundGain = p.gainSurround/p.divConst;
  activation = p.actConst/p.divConst;
  scaledParameters = [centerGain surroundGain activation];
  if any(~isfinite(scaledParameters))
    error('pRF_divNorm:ParameterRatioOverflow', ...
      'Divisive-normalization parameter ratios must remain finite.');
  end
  normalizationDrive = surroundGain*rMinus + 1;
  if any(~isfinite(normalizationDrive)) || any(normalizationDrive <= 0)
    error('pRF_divNorm:InvalidNormalizationDrive', ...
      'The divisive-normalization denominator must remain finite and positive.');
  end
  neuralResponse = ...
    (centerGain*rPlus - activation*surroundGain*rMinus) ./ normalizationDrive;
  thisModelResponse = convolveModelResponseWithHRF(neuralResponse,hrf);

  % drop junk frames here
  thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);

  % return the calculated model response
  output = thisModelResponse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pRF_divNorm('setParams', fitParams) %%%%
%%%        Called from setParams         %%%%

elseif strcmp(varargin{1}, 'setParams')

  fitParams = varargin{2};

  fitParams.paramNames = {'x','y','centerWidth', 'surroundWidth', 'centerAmplitude', 'surroundAmplitude', 'actConst', 'divConst'};
  fitParams.paramDescriptions = {'RF x position','RF y position','RF width of pRF center gaussian (std)', 'Width of surround gaussian (std)', 'Beta weight amplitude of center Gaussian', 'Beta weight amplitude of surround Gaussian', 'Activation constant', 'Divisive constant'};
  fitParams.paramIncDec = [1 1 1 1 1 1 1 1];
  fitParams.paramMin = [-inf -inf 0.1 0.1 0 0 0 0.1];
  fitParams.paramMax = [inf inf inf inf inf inf inf inf];
  % set min/max and init
  fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0.1 0.1 0 0 0 0.1];
  fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf inf inf inf inf inf];
  fitParams.initParams = [0 0 4 4 1 1 10 100];

  % return fitParams with modified values
  output = fitParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  pRF_divNorm(command, fitParams, params)  %%%%%
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
  validateDivNormParameters(p);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   validateDivNormParameters %%
function validateDivNormParameters(p)

values = [p.x p.y p.std p.stdSurround p.gainCenter p.gainSurround p.actConst p.divConst];
if ~isreal(values) || any(~isfinite(values))
  error('pRF_divNorm:InvalidParameters', ...
    'Divisive-normalization parameters must be finite, real scalar values.');
end
if p.std <= 0 || p.stdSurround <= 0
  error('pRF_divNorm:InvalidWidths', ...
    'Divisive-normalization center and surround widths must be greater than zero.');
end
if p.gainCenter < 0 || p.gainSurround < 0 || p.actConst < 0
  error('pRF_divNorm:InvalidAmplitudes', ...
    'Divisive-normalization gains and the activation constant must be nonnegative.');
end
if p.divConst <= 0
  error('pRF_divNorm:InvalidDivisiveConstant', ...
    'The divisive constant must be greater than zero.');
end
