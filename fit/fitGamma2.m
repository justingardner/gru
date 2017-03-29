% fitgamma.m
%
%      usage: fitgamma2(hdr,<time>,<minfit>,<maxfit>,<dispfit>)
%         by: justin gardner
%       date: 08/21/03
%    purpose: fit a gamma function to hemodynamic response
%
function bestfit = fitgamma2(hdr,time,minfit,maxfit,dispfit)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help fitgamma2;
  return
end

if ~exist('time','var') || isempty(time)
  time = 0:(length(hdr)-1);
end
if ~exist('minfit','var')
  minfit = [];
end
if ~exist('maxfit','var')
  maxfit = [];
end
if ~exist('dispfit','var')
  dispfit = 0;
end

if (isempty(minfit))
  minfit = [-10 -inf    -inf -5 -inf  -inf];
end
if (isempty(maxfit))
  maxfit = [10   inf   inf  5  inf  inf];
end

% whether to display errors or not
if (dispfit),displsqnonlin = 'final';,else,displsqnonlin = 'off';,end

% make sure we have a column vector
if (size(hdr,1) ~= length(time))
  hdr = hdr';
end

% how many times function is called
global numcalls;numcalls = 0;
maxiter = inf;
bestresnorm = inf;

% exp1 and exp2 are for the difference of gamma model
% the offset for these is from above
global offset;offset = 0;
global exp1;
global exp2;
exp1 = 5;
exp2 = 5;

bestresnorm = inf;
if (size(hdr,2) ~= 1)
  [maxhdr maxindex] = max(mean(hdr'));
  [minhdr minindex] = min(mean(hdr'));
else
  [maxhdr maxindex] = max(hdr);
  [minhdr minindex] = min(hdr);
end
if (maxhdr > abs(minhdr))
  initparams = [maxhdr 1 2 minhdr 1 time(maxindex(1))];
else
  initparams = [minhdr 1 2 maxhdr 1 time(minindex(1))];
end
if (~dispfit),warning off,end
[fitparams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gammaerr2,initparams,minfit,maxfit,optimset('LevenbergMarquardt','on','MaxIter',maxiter,'Display',displsqnonlin),hdr,time,dispfit);
if (~dispfit),warning backtrace,end

% reduced chi squared is a factor that decreseas the value of the parameter 
% variance estimates according to how many degrees of freedom there are.
reducedChiSquared = (residual*residual')/(prod(size(hdr))-length(initparams));
% this gets the covariance matrix as the inverse of the hessian matrix
covar = reducedChiSquared * inv(jacobian'*jacobian);

if (resnorm < bestresnorm)
  bestfit.params = fitparams;
  bestfit.covar = covar;
  bestfit.minfit = minfit;
  bestfit.maxfit = maxfit;
  bestfit.output = output;
  bestresnorm = resnorm;
end

[bestfit.err bestfit.fit] = gammaerr2(bestfit.params,hdr,time,dispfit);
if (max(bestfit.fit) > abs(min(bestfit.fit)))
  [bestfit.max bestfit.timetomax] = max(bestfit.fit);
  bestfit.timetomax = time(bestfit.timetomax);
else
  [bestfit.max bestfit.timetomax] = min(bestfit.fit);
  bestfit.timetomax = time(bestfit.timetomax);
end

% now make a fit that has more time points in it
bestfit.smoothX = interp1(1:length(hdr),time,1:0.1:length(hdr));
hdrinterp = interp1(1:length(hdr),hdr,1:0.1:length(hdr));
[bogus bestfit.smoothFit] = gammaerr2(bestfit.params,hdrinterp',bestfit.smoothX,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns value of difference of gamma function for x values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = diffgammafun(x,params)

global exp1;
global exp2;
global offset;
y = gammafun(x,[params(1:3) offset exp1])+gammafun(x,[params(4:6) offset exp2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns value of gamma function for x values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = gammafun(x,params)

% if we are passed three parameters then the
% rest come from globals
global exponent;
global offset;
if (length(params) == 3)
  amp = params(1);
  tau = params(2);
  timelag = params(3);
% else all parameters passed in
else
  amp = params(1);
  tau = params(2);
  timelag = params(3);
  offset = params(4);
  exponent = params(5);
end

gammafun = mygamma(x-timelag,exponent,tau);
% make sure gamma isn't all zeros since this will cause
% an error message (if this happens don't bother scaling.
% otherwise set maximum value to 1.
if (max(gammafun)-min(gammafun) ~= 0)
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
y = (amp*gammafun+offset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fittnig gamma to estimated hdrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fitfun] = gammaerr2(fitparams,ehdr,time,dispfit)

% number of times this routine has been called
global numcalls;
numcalls = numcalls + 1;

fitfun = diffgammafun(time,fitparams);
err = [];
for j = 1:size(ehdr,2)
  err = [err (ehdr(:,j)' - fitfun)];
end

if (dispfit && (10*(floor((numcalls-1)/10))==(numcalls-1)))
  global exp1;
  global exp2;
  global offset;
  clf;
  for j = 1:size(ehdr,2)
    plot(time,ehdr(:,j),'k.');
    hold on
  end
  plot(time,fitfun,'r-');

  plot(time,gammafun(time,[fitparams(1:3) offset exp1]),'g');
  plot(time,gammafun(time,[fitparams(4:6) offset exp2]),'c');
  hline(0);
  yaxis(min(-1,-2.5*mean(std(ehdr'))),max(1.5,2.5*mean(std(ehdr'))));
  title(sprintf('amp=%0.2f tau=%0.2f lag=%0.2f amp=%0.2f tau=%0.2f lag=%0.2f funcalls %i',fitparams(1),fitparams(2),fitparams(3),fitparams(4),fitparams(5),fitparams(6),numcalls));
  xlabel('time (sec)');
  ylabel('BOLD response (%)');
  drawnow
end