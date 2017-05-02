% fitsigmoid
%
%      usage: fitsigmoid(contrasts,responses,minfit,maxfit)
%         by: justin gardner
%       date: 10/07/03
%    purpose: fit a sigmoid function to data
%
function bestfit = fitsigmoid(contrasts,responses,minfit,maxfit)

% check arguments
if (nargin == 2)
  minfit = [0   0   -inf -inf];
  maxfit = [inf 1 inf inf];
elseif (nargin ~= 4)
  help fitsigmoid
  return
end

% make sure we have a column vector
if (size(contrasts,1) == 1)
  contrasts = contrasts';
end
if (size(responses,1) == 1)
  responses = responses';
end

maxiter = inf;
bestresnorm = inf;

% these will be set in the fits (not allowed to fix)
% so there value is passed to the objective function
% as these globals
global n;
% set fixed parameters
n = 2;

% set initparams
%initparams = [1 median(contrasts) 0];
%initparams = [1 rand(1) 0];
initparams = [1 median(contrasts) 1.5 0.2];
%minfit = [0 0 0 -inf];
%maxfit = [inf inf inf inf];

global numIters;numIters = 0;
%initparams(2) = 1000;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitparams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@sigmoiderr,initparams,minfit,maxfit,[],contrasts,responses);

% FIXFIXFIXFIX
% Taken from Numerical Recipies, 
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
jacobian = jacobian'*jacobian;
reducedChiSquared = (residual'*residual)/(length(responses)-length(initparams));
covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit.params = fitparams;
bestfit.covar = covar;
bestfit.output = output;

numIters = nan;
[bestfit.err bestfit.fit] = sigmoiderr(bestfit.params,contrasts,responses);

bestfit.fitx = 0:.0001:1;
[dummy bestfit.fity] = sigmoiderr(bestfit.params,bestfit.fitx,bestfit.fitx);
bestfit.r2 = 1-var(bestfit.err)/var(responses);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fittnig gamma to estimated hdrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fitfun] = sigmoiderr(fitparams,contrasts,responses)

% if we are passed three parameters then the
% rest come from globals
if (length(fitparams) == 3)
  global n;
  amp = fitparams(1);
  c50 = fitparams(2);
  offset = fitparams(3);
else
  amp = fitparams(1);
  c50 = fitparams(2);
  n = fitparams(3);
  offset = fitparams(4);
end

% calculate function
fitfun = amp*(contrasts.^n)./(contrasts.^n+c50.^n)+offset;

global numIters;
numIters = numIters+1;
% display the fit as we go along
if 0%~isnan(numIters)
  [sortcontrasts sortindex] = sort(contrasts);
  clf;semilogx(100*contrasts,responses,'ko');hold on
  semilogx(100*sortcontrasts,fitfun(sortindex),'r-');
  vline(c50*100);
  if (length(fitparams) == 4)
    title(sprintf('amp=%0.2f, c50=%0.2f, n=%0.2f, offset=%0.2f numIters=%i',fitparams(1),100*fitparams(2),fitparams(3),fitparams(4),numIters));
  else
    title(sprintf('amp=%0.2f, c50=%0.2f, offset=%0.2f numIters=%i',fitparams(1),100*fitparams(2),fitparams(3),numIters));
  end
    drawnow
end
% calculate error
err = responses-fitfun;