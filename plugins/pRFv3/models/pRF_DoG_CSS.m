% pRF_DoG_CSS.m
%
%        $Id:$ 
%      usage: pRF_DoG_CSS(varargin)
%         by: akshay jagadeesh
%       date: 09/02/16
%    purpose: Model of Difference of Gaussians with Compressive Static Non linearity 
%
%             - Allows you to create a new model type 
%             - Simply add the model type to pRFGUI and create a 
%               file like this one for your model with the appropriate 
%               methods filled in.
%
%     This model template is set up for the standard Gaussian model.

function output = pRF_DoG_CSS(varargin)

if nargin < 2
  disp(sprintf('Not enough arguments'));
  disp(sprintf('Number of arguments: %d', nargin));
  celldisp(varargin)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pRF_DoG_CSS('getModelResponse', fitParams, rfModel, hrf, p, i) %%%%
%%%                  Called from getModelResidual                   %%%%
if strcmp(varargin{1}, 'getModelResponse')

  fitParams = varargin{2};
  rfModel1 = varargin{3};
  hrf = varargin{4};
  p = varargin{5};
  i = varargin{6};

  validateDoGCSSParameters(p);

  rfModel2 = pRFMakeGaussian(fitParams,p.x,p.y,p.std2);
  rPlus = pRFProjectStimulus(rfModel1,fitParams.stim{i});
  rMinus = pRFProjectStimulus(rfModel2,fitParams.stim{i});

  % Form the neural center-surround response before the hemodynamic
  % convolution. A signed power retains suppressive surround responses while
  % keeping fractional exponents real-valued.
  neuralResponse = rPlus*p.B1 + rMinus*p.B2;
  neuralResponse = sign(neuralResponse).*power(abs(neuralResponse),p.exp);
  thisModelResponse = convolveModelResponseWithHRF(neuralResponse,hrf);

  % drop junk frames here
  thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);

  % return the calculated model response
  output = thisModelResponse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pRF_DoG_CSS('setParams', fitParams) %%%%
%%%        Called from setParams         %%%%

elseif strcmp(varargin{1}, 'setParams')

  fitParams = varargin{2};

  fitParams.paramNames = {'x','y','rfWidth', 'surroundWidth', 'centerAmplitude', 'surroundAmplitude', 'exponent'};
  fitParams.paramDescriptions = {'RF x position','RF y position','RF width (std of gaussian)', 'RF width of surround pool', 'Beta weight amplitude of positive Gaussian', 'Beta weight amplitude for negative Gaussian', 'Exponent of the signed static nonlinearity'};
  fitParams.paramIncDec = [1 1 1 1 1 1 1];
  fitParams.paramMin = [-inf -inf 0.1 0.1 0 -inf 0.1];
  fitParams.paramMax = [inf inf inf inf inf 0 1];
  % set min/max and init
  fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0.1 0.1 0 -inf 0.1];
  fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf inf inf 0 1];
  fitParams.initParams = [0 0 4 4 1 0 1];

  % return fitParams with modified values
  output = fitParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  pRF_DoG_CSS(command, fitParams, params)  %%%%%
%%%%         Called from getFitParams           %%%%% 

elseif strcmp(varargin{1}, 'getFitParams')

  fitParams = varargin{2};
  params = varargin{3};

  p.rfType = fitParams.rfType;

  % Define your parameters here
  p.x = params(1);
  p.y = params(2);
  p.std = params(3);
  p.std2 = params(4);
  p.B1 = params(5);
  p.B2 = params(6);
  p.exp = params(7);
  validateDoGCSSParameters(p);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   validateDoGCSSParameters   %%
function validateDoGCSSParameters(p)

values = [p.x p.y p.std p.std2 p.B1 p.B2 p.exp];
if ~isreal(values) || any(~isfinite(values))
  error('pRF_DoG_CSS:InvalidParameters', ...
    'DoG-CSS parameters must be finite, real scalar values.');
end
if p.std <= 0 || p.std2 <= 0
  error('pRF_DoG_CSS:InvalidWidths', ...
    'DoG-CSS center and surround widths must be greater than zero.');
end
if p.B1 < 0 || p.B2 > 0
  error('pRF_DoG_CSS:InvalidAmplitudes', ...
    'DoG-CSS requires a nonnegative center amplitude and a nonpositive surround amplitude.');
end
if p.exp <= 0
  error('pRF_DoG_CSS:InvalidExponent', ...
    'The DoG-CSS exponent must be greater than zero.');
end
