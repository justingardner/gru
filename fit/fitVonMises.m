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
function retval = fitVonMises(x,y,varargin)

% check arguments
if nargin < 2
  help fitVonMises
  return
end

% default retval
retval = [];

% check arguments
getArgs(varargin,{'dispfit=1','mu=0','kappa=5','amp=1','offset=0','halfWidthAtHalfHeight=[]'});

% convert halfWidthAtHalfHeight to kappa
if ~isempty(halfWidthAtHalfHeight)
  % get cos of the desired angle
  cosTheta = cos(d2r(halfWidthAtHalfHeight));
  % then solve the equation to find what kappa makes it that that theta
  % value gives a value of 0.5 (init at kappa of 1) - don't think there
  % is easy closed form expression so use fzero to search for kappa
  retval = fzero(@(kappa) ((exp(kappa*cosTheta) - exp(-kappa))/(exp(kappa)-exp(-kappa)) - 0.5),1)
  return
end

% convert to radians
x = d2r(x(:));mu = d2r(mu);
y = y(:);

% initial parameters
initParams = setVonMisesParams(mu,kappa,amp,offset);

% check for empty x
if isempty(x)
  retval = r2d(acos(log(((exp(kappa)-exp(-kappa))/2)+exp(-kappa))/kappa));
  return
end
  
% check for empty y, then just compute vonMises
if isempty(y)
  retval = myVonMises(initParams,x);
  return
end

% check lengths
if length(x) ~= length(y)
  disp(sprintf('(fitVonMises) Length of x (%i) must match y (%i)',length(x),length(y)));
  return
end

% set options for fminsearch
options = optimset('MaxFunEvals',inf);

% do fit using nelder-mead
optimParams = fminsearch(@vonMisesError,initParams,options,y,x,0,dispfit);

% fit again with levenberg-marquardt (to get jacobian)
options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',inf);
[fitVals,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@vonMisesError,optimParams,[],[],options,y,x,1,dispfit);

% reduced chi squared is a factor that decreseas the value of the
% parameter variance estimates according to how many degrees of
% freedom there are.
reducedChiSquared = (residual(:)'*residual(:))/(length(y)-length(initParams));
% this gets the covariance matrix as the inverse of the hessian matrix
covar = reducedChiSquared * inv(jacobian'*jacobian);

% save params in return structure
retval.optimParams = optimParams;
[retval.params.mu retval.params.kappa retval.params.amp retval.params.offset] = getVonMisesParams(optimParams);
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
if dispfit
  f = mlrSmartfig('vonMisesFit','reuse=1');
  clf;
  plot(x,y,'k.');
  hold on
  plot(retval.xSmooth,retval.yFitSmooth,'r-');
  title(sprintf('%s: %0.4f',num2str(retval.optimParams),sqrt(sum(residual.^2))));
  drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    get error between von mises and data    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = vonMisesError(params,y,x,returnResidual,dispfit)

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
if dispfit
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
function y = myVonMises(params,x)

% get the parameters
[mu kappa amp offset] = getVonMisesParams(params);

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




