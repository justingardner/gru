% fitCumulativeGaussian
%
%      usage: fitCumulativeGaussian(x,y)
%         by: justin gardner
%       date: 10/11/2019
%    purpose: fit a cumulative gaussian to psychophysical data
%             y should be correct proportions (values between 0 and 1)
%    options: 'minParams=[-inf 0 0]' The minimum values parameters can go to
%             'maxParams=[inf inf inf] The maximum values parameters can go to
%             'initParams=[]' Starting values of parameters - defaults to median and std of x values
%             'maxIter=inf' How many iterations to search for best parameters
%             'fitType=basic' Fits just mean and standard deviation. Set to 'lapse' to also fit a lapse rate
%             'guessRate=0' Sets the guess rate, if you have a 2AFC you should set to 0.5
%
%       e.g.: 
%             x = [1 5 10 20 30 40];
%             y = [0.02 0.1 0.15 0.3 0.7 0.9];
%             fit = fitCumulativeGaussian(x,y,'fitType=lapse')
%             clf; plot(x,y,'ko'); hold on; plot(fit.fitX,fit.fitY,'r-');
%             xlabel('stimulus value');ylabel('Correct (proportion)');
%             title(sprintf('Mean %f Std %f lambda %f',fit.mean,fit.std,fit.lambda));
%              
%
%
function bestfit = fitCumulativeGaussian(x,y,varargin)

% check arguments
if nargin < 2
  help fitCumulativeGaussian;
end

% parse input arguments
getArgs(varargin,{'minParams',[-inf 0 0],'maxParams',[inf inf inf],'initParams',[],'maxIter=inf','fitType=basic','guessRate=0'});

% make sure we have a column vector
x = x(:)';y = y(:)';

% set the initial parameters
% params are [mean std lambda] 
if isempty(initParams)
  initParams = [median(x) std(x) 0];
end

% number of parameters are different with different types of fit
switch (fitType)
 case 'basic'
  % just mean / std
  nParams = 2;
 case 'lapse'
  % has lapse rate
  nParams = 3;
 otherwise
  disp(sprintf('(fitCumulativeGaussian) Unkown fitType: %s',fitType));
  return
end

% make parameters the right length
initParams = initParams(1:nParams);
minParams = minParams(1:nParams);
maxParams = maxParams(1:nParams);
  
% set optimization parametrs
optimParams = optimset('MaxIter',maxIter);

% some globals to keep track of what lsqnonlin does
global numIters;numIters = 0;

% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@cumulativeGaussianErr,initParams,minParams,maxParams,optimParams,x,y,fitType,guessRate);

% Taken from Numerical Recipies, 
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
jacobian = jacobian'*jacobian;
reducedChiSquared = (residual*residual')/(length(y)-length(initParams));
covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit = extractParams(fitParams,fitType);
bestfit.params = fitParams;
bestfit.covar = covar;
bestfit.output = output;

% get the fit
[bestfit.err bestfit.fit] = cumulativeGaussianErr(bestfit.params,x,y,fitType,guessRate);

% compute r2 of fit
bestfit.r2 = 1-var(bestfit.err)/var(y);

% compute a smoother fit (i.e. nFitPoints along x axis)
nFitPoints = 1000;
bestfit.fitX = min(x):(max(x)-min(x))/(nFitPoints-1):max(x);

% note here that the y variable is just a dummy value since we don't
% care about the error 
[~,bestfit.fitY] = cumulativeGaussianErr(bestfit.params,bestfit.fitX,zeros(1,nFitPoints),fitType,guessRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cumulativeGaussianErr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fit] = cumulativeGaussianErr(fitParams,x,y,fitType,guessRate)

% get the parmas
p = extractParams(fitParams, fitType);

% calculate the gaussian
fit = guessRate+p.lambda+(normcdf(x,p.mean,p.std)*((1-guessRate)-2*p.lambda));

% update number of iterations
global numIters;
numIters = numIters+1;

err = y-fit;

%%%%%%%%%%%%%%%%%%%%%%%
%    extractParams    %
%%%%%%%%%%%%%%%%%%%%%%%
function p = extractParams(fitParams, fitType)

% extrat the parameters
if strcmp(fitType,'basic')
  p.mean = fitParams(1);
  p.std = fitParams(2);
  p.lambda = 0;
elseif strcmp(fitType,'lapse')
  p.mean = fitParams(1);
  p.std = fitParams(2);
  p.lambda = fitParams(3);
else
  disp(sprintf('(fitCumulativeGaussian:extractParams) Unknown fitType: %s',fitType));
  keyboard
end