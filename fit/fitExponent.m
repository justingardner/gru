% fitexponent.m
%
%      usage: fitexponent(x,y,dispfit)
%         by: justin gardner
%       date: 09/30/04
%    purpose: fit exponential curve of the following form to y
%             
%             y = a*e(-x/tau)+offset
%
%       e.g.: fitparams = fitexponent(time,data);
%             to display fits set dispfit to 1
%             fitparams = fitexponent(time,data,1);
function retval = fitexponent(x,y,dispfit)

if (nargin == 2)
  dispfit = 0;
elseif (nargin ~= 3)
  help fitexponent
  return
end

% make into row vectors
if (size(x,1) ~= 1),x = x';end
if (size(y,1) ~= 1),y = y';end

% set some of the parameters for the fitting
maxiter = inf;
minfit = [];
maxfit = [];
if (dispfit),displsqnonlin = 'final';,else,displsqnonlin = 'off';,end
initparams = [max(y)-min(y) -median(x) min(x)];

% get the fit
[fitparams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@experr,initparams,minfit,maxfit,optimset('Display','iter','LevenbergMarquardt','on','MaxIter',maxiter,'Display',displsqnonlin),x,y);

% reduced chi squared is a factor that decreseas the value of the
% parameter variance estimates according to how many degrees of
% freedom there are.
reducedChiSquared = (residual*residual')/(length(y)-length(initparams));
% this gets the covariance matrix as the inverse of the hessian matrix
covar = reducedChiSquared * inv(jacobian'*jacobian);

retval.fitparams = fitparams;
retval.covar = covar;
retval.a = fitparams(1);
retval.avar = covar(1,1);
retval.apervar = 100*retval.avar/retval.a;
retval.tau = fitparams(2);
retval.tauvar = covar(2,2);
retval.taupervar = 100*retval.tauvar/retval.tau;
retval.offset = fitparams(3);
retval.offsetvar = covar(3,3);
retval.offsetpervar = 100*retval.offsetvar/retval.offset;
retval.x = x;
retval.y = y;
retval.fit = fitparams(1)*exp(-x/fitparams(2)) + fitparams(3);
retval.r2 = 1-var(residual)/var(y);
if (dispfit)
  plot(x,y,'ko');
  hold on
  plot(x,retval.fit,'k-');
  xlabel('x');
  ylabel(sprintf('%0.2f * exp(-x/%0.2f) + %0.2f',fitparams(1),fitparams(2),fitparams(3)));
  title(sprintf('amp=%0.2f (%0.2f%%) tau=%0.2f (%0.2f%%) offset=%0.2f (%0.2f%%)',retval.a,retval.apervar,retval.tau,retval.taupervar,retval.offset,retval.offsetpervar));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fitting exponent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err] = experr(params,x,y)

model = params(1)*exp(-x/params(2)) + params(3);
err = y-model;

