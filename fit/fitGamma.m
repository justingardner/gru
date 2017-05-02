% fitgamma.m
%
%      usage: fitgamma(hdr,<time>,<minfit>,<maxfit>,<dispfit>)
%         by: justin gardner
%       date: 08/21/03
%    purpose: fit a gamma function to hemodynamic response
%             the parameters returned are
%             [amp tau timelag offset exponent]
%             offset is usually forced to zero
%             exponent is either 5 or 6
%             minfit is an array of length three specifying
%             the minimum values for amp tau and timelag respectively
%             maxfit is the maximum
%
%        e.g.:
%tr = 0.8;
%time = tr:tr:25*tr;
%hrf = mygamma(time-3,5,1);hrf = hrf/max(hrf);
%hrfnoise = hrf + rand(1,length(hrf))/2;
%bestfit = fitgamma(hrfnoise,time,[],[],1);
%
function bestfit = fitgamma(hdr,time,minfit,maxfit,dispfit)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help fitgamma
  return
end

if ~exist('time','var') || isempty(time)
  time = 0:(length(hdr)-1);
end
if ~exist('minfit','var')
  minfit = [];
end
if ~exist('maxfit','var')
  maxfit = [inf inf 5 inf inf];
end
if ~exist('dispfit','var')
  dispfit = 0;
end
% whether to display errors or not
if (dispfit),displsqnonlin = 'final';,else,displsqnonlin = 'off';,end

% make sure we have a column vector
if (size(hdr,1) == 1)
  hdr = hdr';
end

global numcalls;numcalls = 0;
maxiter = inf;
bestresnorm = inf;

% these will be set in the fits (not allowed to fix)
% so there value is passed to the objective function
% as these globals
global exponent;
global offset;

for i = 5:6
  % set fixed parameters
  exponent = i;
  offset = 0;
  % set initparams
  initparams = [max(hdr) 2 0];
  % fit function using lsqnonlin in LevenbergMarquardt mode.
  [fitparams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gammaerr,initparams,minfit,maxfit,optimset('LevenbergMarquardt','on','MaxIter',maxiter,'Display',displsqnonlin),hdr,time,dispfit);

  % Taken from Numerical Recipies, 
  % the leastsq function seems to return the transposed gradient
  % instead of the jacobian...
  jacobian = jacobian'*jacobian;
  reducedChiSquared = (residual*residual')/(length(hdr)-length(initparams));
  covar = sqrt(reducedChiSquared * inv(jacobian));

  % if we have the best fit then keep it.
  if (resnorm < bestresnorm)
    bestfit.params = [fitparams(1) fitparams(2) fitparams(3) offset exponent];
    bestfit.covar = covar;
    bestfit.output = output;
    bestresnorm = resnorm;
  end
end

[bestfit.err bestfit.fit] = gammaerr(bestfit.params,hdr,time,dispfit);

% now make a fit that has more time points in it
bestfit.smoothX = interp1(1:length(hdr),time,1:0.1:length(hdr));
hdrinterp = interp1(1:length(hdr),hdr,1:0.1:length(hdr));
[bogus bestfit.smoothFit] = gammaerr(bestfit.params,hdrinterp',bestfit.smoothX,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fittnig gamma to estimated hdrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fitfun] = gammaerr(fitparams,ehdr,time,dispfit)

% number of times this routine has been called
global numcalls;
numcalls = numcalls + 1;

% if we are passed three parameters then the
% rest come from globals
global exponent;
global offset;
if (length(fitparams) == 3)
  amp = fitparams(1);
  tau = fitparams(2);
  timelag = fitparams(3);
% else all parameters passed in
elseif (length(fitparams) == 4)
  amp = fitparams(1);
  tau = fitparams(2);
  timelag = fitparams(3);
  offset = fitparams(4);
else
  amp = fitparams(1);
  tau = fitparams(2);
  timelag = fitparams(3);
  offset = fitparams(4);
  exponent = round(fitparams(5));
end

% calculate function
gammafun = mygamma(time-timelag,exponent,tau);
gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
fitfun = (amp*gammafun+offset);
err = ehdr' - fitfun;

if (dispfit && (10*(floor((numcalls-1)/10))==(numcalls-1)))
clf;
  plot(time,ehdr(:),'k.');
  hold on
  plot(time,fitfun,'r-');

  hline(0);
  yaxis(min(-1,-2.5*mean(std(ehdr'))),max(1.5,2.5*mean(std(ehdr'))));
  title(sprintf('amp=%0.2f tau=%0.2f lag=%0.2f funcalls=%i',fitparams(1),fitparams(2),fitparams(3),numcalls));
  xlabel('time (sec)');
  ylabel('BOLD response (%)');
  drawnow
end