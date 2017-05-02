% fitNakaRushton.m
%
%        $Id:$ 
%      usage: fit = fitNakaRushton(c,r)
%         by: justin gardner
%       date: 12/24/13
%    purpose: fit a naka-rushton to data (taken from selectionModel)
%
function fit = fitNakaRushton(c,r,varargin)

% check arguments
if nargin < 2
  help fitNakaRushton
  return
end

% parse arguments
getArgs(varargin,{'dispFit=0','evalFit=[]'});

% find contrast that evokes closest to half-maximal response
rMid = ((max(r)-min(r))/2) + min(r);
[dummy,rMidIndex] = min(abs(r-rMid));
initC50 = c(rMidIndex(1));

% parmaeters
             %Rmax          c50     n     offsets x 5
initParams = [max(r)        initC50 2  min(r)];
minParams =  [0             0       1  -inf];
maxParams =  [inf           1       5  inf];

% set model type
m.fixedN = 0;
m.dispFit = dispFit;

% optimization parameters
maxiter = inf;
optimParams = optimset('MaxIter',maxiter,'Display','off');

% now go fit
[params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@nakaRushtonResidual,initParams,minParams,maxParams,optimParams,c,r,m);

% parse params and return
fit = parseParams(params,m);

% if evalFit is set then eval at every contrast specified in that variable
if ~isempty(evalFit)
  fit.cFit = evalFit;
  fit.rFit = nakaRushton(evalFit,fit);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    nakaRushtonResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = nakaRushtonResidual(params,c,r,m)

% decode parameters
p = parseParams(params,m);
% calculate naka-rushton
fitR = nakaRushton(c,p); 

% display fit if called for
if m.dispFit
  f = smartfig('selectionModel_nakaRushtonResidual','reuse');  
  clf;
  semilogx(c,r,'ko');
  hold on
  semilogx(c,fitR,'k-')
  titleStr = sprintf('Rmax: %0.3f c50: %0.2f n: %0.3f\n',p.Rmax,p.c50,p.n);
  title(sprintf('%s offset: %f',titleStr,p.offset));
  drawnow
end

residual = r(:)-fitR(:);
residual = residual(:);

%%%%%%%%%%%%%%%%%%%%%
%    nakaRushton    %
%%%%%%%%%%%%%%%%%%%%%
function response = nakaRushton(c,p)

response = p.Rmax * ((c.^p.n) ./ ((c.^p.n) + p.c50.^p.n)) + p.offset;

%%%%%%%%%%%%%%%%%%%%%
%    parseParams    %
%%%%%%%%%%%%%%%%%%%%%
function p = parseParams(params,m)

if m.fixedN
  p.Rmax = params(1);
  p.c50 = params(2);
  p.n = m.n;
  p.offset = params(3);
else
  p.Rmax = params(1);
  p.c50 = params(2);
  p.n = params(3);
  p.offset = params(4);
end  


