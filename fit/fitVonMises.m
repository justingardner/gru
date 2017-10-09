% fitVonMises.m
%
%      usage: fit = fitVonMises(x,y,<dispfit=0>)
%         by: justin gardner
%       date: 03/28/17
%    purpose: x is in degrees, y is response amplitude
%             returns vonMises fit parameters
%
%       e.g.: 
%             x = [-90:22.5:89];
%             y = [0.088 0.091 0.124 0.195 0.279 0.190 0.112 0.084];
%             fitVonMises(x,y,'dispfit=1');
%
%             to set init parameters, note that mu is in degrees (not radians);
%             fitVonMises(x,y,'dispfit=1','mu=30','kappa=5','amp=1','offset=0');
%
%             if y is empty, will compute the VonMises with the specified init params. e.g.
%             x = 0:360;
%             y = fitVonMises(x,[],'mu=0','kappa=5','amp=1','offset=0');
%
%             You can convert kappa parameters to half-width-at-half-height in degrees
%             by setting x to empty:
%
%             halfWidth = fitVonMises([],[],'kappa=3');
% 
%             or the inverse (where half-width-at-half-height is in deg)
%
%             kappa = fitVonMises([],[],'halfWidthAtHalfHeight=30');
%
%             can also do a parametric bootstrap fit by randomly sampling
%             y values from the gaussian distribution with mean y and standard
%             deviation specified by yste. This will be done parametricBootrap number
%             of times and the parameters returned in bootParams:
%
%             x = [-90:22.5:89];
%             y = [0.088 0.091 0.124 0.195 0.279 0.190 0.112 0.084];
%             yste = [0.016 0.020 0.019 0.011 0.031 0.015 0.016 0.013];
%
%             fit = fitVonMises(x,y,'yste',yste,'parametricBootstrap=1000','dispFit=0')
%
function retval = fitVonMises(x,y,varargin)

% check arguments
if nargin < 2
  help fitVonMises
  return
end

% default retval
retval = [];

% check arguments (this was really slow, so providing a fixed way to call in
% which arguments just get passed in)
if nargin > 2
  if isstr(varargin{1})
    getArgs(varargin,{'dispFit=2','mu=0','kappa=5','amp=1','offset=0','halfWidthAtHalfHeight=[]','wrapAt=360','parametricBootstrap=0','yste=[]'});
  else
    kappa = 5;amp = 1; offset = 0; halfWidthAtHalfHeight = [];wrapAt = 360;dispFit = 2;yste=[];parametricBootstrap = 0;
    mu = varargin{1};
    if nargin > 3,kappa = varargin{2};end
    if nargin > 4,wrapAt = varargin{3};end
  end
end

% set how much info to display on fits
if dispFit > 0
  displayType = 'final';
else
  displayType = 'off';
end

% convert halfWidthAtHalfHeight to kappa
if ~isempty(halfWidthAtHalfHeight)
  % convert to kappa
  retval = halfWidthAtHalfHeight2kappa(halfWidthAtHalfHeight);
 
  % if x is not empty, then it means we are getting y-values so just set kappa
  if ~isempty(x)
    kappa = retval(1);
  else
    % otherwise return with the kappa values
    return
  end
end

% convert to radians
x = d2r(x(:));mu = d2r(mu);
y = y(:);

% initial parameters
initParams = setVonMisesParams(mu,kappa,amp,offset);

% check for empty x
if isempty(x)
  % convert to half-width-at-half-height
  retval = kappa2halfWidthAtHalfHeight(kappa);
  return
end
  
% check for empty y, then just compute vonMises
if isempty(y)
  % compute von mises, notice extra parameter
  % is so that the function can wrap not at 360 deg
  % but at smaller (e.g. 180 - useful for orientations).
  retval = myVonMises(initParams,x,360/wrapAt);
  return
end

% check lengths
if length(x) ~= length(y)
  disp(sprintf('(fitVonMises) Length of x (%i) must match y (%i)',length(x),length(y)));
return
end

% set options for fminsearch
options = optimset('MaxFunEvals',inf,'Display',displayType);

% do fit using nelder-mead
optimParams = fminsearch(@vonMisesError,initParams,options,y,x,0,dispFit);

% fit again with levenberg-marquardt (to get jacobian)
options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',inf,'Display',displayType);
[fitVals,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@vonMisesError,optimParams,[],[],options,y,x,1,dispFit);

% reduced chi squared is a factor that decreseas the value of the
% parameter variance estimates according to how many degrees of
% freedom there are.
reducedChiSquared = (residual(:)'*residual(:))/(length(y)-length(initParams));
% this gets the covariance matrix as the inverse of the hessian matrix
covar = reducedChiSquared * inv(jacobian'*jacobian);

% save params in return structure
retval.optimParams = optimParams;
[retval.params.mu retval.params.kappa retval.params.amp retval.params.offset] = getVonMisesParams(optimParams);
retval.params.halfWidthAtHalfHeight = kappa2halfWidthAtHalfHeight(retval.params.kappa);
retval.covar = covar;
retval.residual = residual;
retval.squaredError = sum(residual.^2);
% save function fit
retval.yFit = myVonMises(optimParams,x);
x = r2d(x);
retval.xSmooth = min(x):max(x);
retval.yFitSmooth = myVonMises(optimParams,d2r(retval.xSmooth));
retval.x = x;
retval.y = y;

% display fit
if dispFit
  f = mlrSmartfig('vonMisesFit','reuse=1');
  clf;
  plot(x,y,'k.');
  hold on
  plot(retval.xSmooth,retval.yFitSmooth,'r-');
  title(sprintf('%s: %0.4f',num2str(retval.optimParams),sqrt(sum(residual.^2))));
  drawnow
end

% do parameteric bootstrap, y sampling data from the gaussian implied by 
% the yste
if parametricBootstrap
  if isempty(yste)
    disp(sprintf('(fitVonMises) yste must be specified for parameteric bootstrap'));
    keyboard
  end
  
  retval.bootParams = [];
  disppercent(-inf,sprintf('(fitVonMises) Doing bootstrap'));
  for iBoot = 1:parametricBootstrap
    % sample data from gaussian distribution with 
    % mean and var specified by y and yste
    yBoot = mvnrnd(y,diag(yste.^2));
    % fit
    retval.bootParams(iBoot,:) = fminsearch(@vonMisesError,optimParams,options,yBoot,d2r(x),0,0);
    retval.yBoot(iBoot,:) = yBoot;
    % get boot params half width
    [~,kappa] = getVonMisesParams(retval.bootParams(iBoot,:));
    retval.yBootHalfWidth(iBoot) = kappa2halfWidthAtHalfHeight(kappa);
    % update disppercent
    disppercent(iBoot/parametricBootstrap);
  end
  disppercent(inf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    get error between von mises and data    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = vonMisesError(params,y,x,returnResidual,dispFit)

% compute vonMises 
yFit = myVonMises(params,x);

% get fit error
residual = y(:) - yFit(:);
squaredError = sum(residual.^2);

% return residual or squared error
if returnResidual
  retval = residual;
else
  retval = squaredError;
end

% display the fit
if dispFit
  f = mlrSmartfig('vonMisesFit','reuse=1');
  clf;
  plot(r2d(x),y,'k.');
  hold on
  plot(r2d(x),yFit,'r-');
  yaxis(0,0.4);
  title(sprintf('%s: %0.4f',num2str(params),squaredError));
  drawnow
end

%%%%%%%%%%%%%%%%%%%%
%    myVonMises    %
%%%%%%%%%%%%%%%%%%%%
function y = myVonMises(params,x,wrapScaleFactor)

% the wrapScaleFactor is used to make the x-range go
% to a different value - for example to 180 degrees
% for orientations. To achieve that we scale 
if nargin == 2
  wrapScaleFactor = 1;
end
  
% get the parameters
[mu kappa amp offset] = getVonMisesParams(params);

% if we are going to rescale the x-axis, but we want
% to make sure that kappa means the same thing
% in terms of halfWidthAtHalfHeight
if ~isequal(wrapScaleFactor,1)
  halfWidthAtHalfHeight = kappa2halfWidthAtHalfHeight(kappa);
  % rescale kappa
  kappa = halfWidthAtHalfHeight2kappa(halfWidthAtHalfHeight*wrapScaleFactor);
  % now rescale
  x = x*wrapScaleFactor;
  mu = mu*wrapScaleFactor;
end

% von mises normalized so that amplitude is the difference between the lowest and highest points 
y = offset + (amp-offset) * (exp(kappa * cos(x - mu)) - exp(-kappa))/(exp(kappa)-exp(-kappa));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    convert von mises parameters into values    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu kappa amp offset] = getVonMisesParams(params)

mu = params(1);
kappa = params(2);
amp = params(3);
offset = params(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    convert von mises parameters into array    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = setVonMisesParams(mu,k,amp,offset)

params = [mu k amp offset];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    kappa2halfWidthAtHalfHeight    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function halfWidthAtHalfHeight = kappa2halfWidthAtHalfHeight(kappa)

if isinf(kappa)
  halfWidthAtHalfHeight = 0;
else
  halfWidthAtHalfHeight = r2d(acos(log(((exp(kappa)-exp(-kappa))/2)+exp(-kappa))./kappa));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    halfWidthAtHalfHeight2kappa    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kappa = halfWidthAtHalfHeight2kappa(halfWidthAtHalfHeight)

for iWidth = 1:length(halfWidthAtHalfHeight)
  if halfWidthAtHalfHeight(iWidth) == 0
    kappa(iWidth) = inf;
  else
    % get cos of the desired angle
    cosTheta = cos(d2r(halfWidthAtHalfHeight(iWidth)));
    % then solve the equation to find what kappa makes it that that theta
    % value gives a value of 0.5 (init at kappa of 1) - don't think there
    % is easy closed form expression so use fzero to search for kappa
    kappa(iWidth) = fzero(@(kappa) ((exp(kappa*cosTheta) - exp(-kappa))/(exp(kappa)-exp(-kappa)) - 0.5),1);
  end
end
