% fitweibull.m
%
%      usage: fitweibull(signal,ncorrect,<'ntotal=[]'>,<'gamma=0.5'>,<'dispfig=0'>,<'dogoodnessoffit=0'>,<'dobootstrap=0'>,<'maxIter=10000'>,<'bootmode='binned'>)
%         by: justin gardner
%       date: 06/08/05
%    purpose: fit weibull distribution to psychophysical data
%             signal is signal strength
%             ncorrect is number of trials correct for each signal strength
%             ntotal is total trials for each signal strength
%             gamma (rate of correct when guessing - can be set to nan if you want to have it fit) 
%             is set to 0.5 (for 2AFC) if not given
%             maxIter is the maximum number of iterations to run the nelder-mead search routine for
%             See:
%             bootmode is the type of bootstrap. 'binned' means we resample
%             data within each level (preserves the overall sample
%             size at each signal level), 'full' means we resample single
%             trials from the entire dataset with replacement. 
%             Wichmann & Hill (2001) Percept Psychophys 63:1314-1329
%       e.g.: fitweibull(signal,ncorrect,'ntotal',ntotal);
%             to get a goodness of fit estimate with monte-carlo
%             simulation or bootstrap confidence intervals and with figure
%             fitweibull(signal,ncorrect,'ntotal',ntotal,'gamma',0.5,'dispfig',1,'dogoodnessoffit',1000,'dobootstrap',1000)
%
%             if ntotal = [], then ncorrect is interperted as the
%             responses of the subject (1=correct, 0=incorrect,-1=no response)
%             and signal contains the signal strength that was
%             tested. then the program computes percent correct
%  
%       e.g.: % create some fake data and test
% 
%             % now get some random signal strengths
%             n = 100;
%             s = rand(1,n);
%
%             % get the probability correct for each signal strength
%             p = weibull(s,[0.5 3.5 0.01 0.5]);
%
%             % generate data
%             correctIncorrect = p > rand(1,n);
%
%             % now fit the data
%             fit = fitweibull(s,correctIncorrect,'gamma',0.5,'dispfig',1);
%
function bestfit = fitweibull(varargin)

% parse arguments. This accepts the old style of calling which was the arguments:
% (signal,ncorrect,ntotal,gamma,dispfig,dogoodnessoffit,dobootstrap) in that order.
% Now you can use the getArgs system for setting arguments (e.g. fitweibull(signal,ncorrect,'dispfig',1) and
% so forth.
[signal ncorrect ntotal gamma dispfig dogoodnessoffit dobootstrap maxIter bootmode] = parseArgs(varargin);
if isempty(signal),return,end

% calculate the ncorrect and ntotal if we are passed in [] as
% ntotal
if (isempty(ntotal))
  uniqsignal = unique(sort(signal));
  response = ncorrect;ncorrect = [];
  % count the number of correct and total for each signal strength
  for i = 1:length(uniqsignal)
    thissig = find(signal == uniqsignal(i));
    ncorrect(i) = sum(response(thissig)==1);
    ntotal(i) = sum(response(thissig)~=-1);
  end
  signal = uniqsignal;
  % now display what we found if asked for
  if (dispfig)
    for i = 1:length(uniqsignal)
      mydisp(sprintf('%0.2f\t',signal(i)));
    end
    mydisp(sprintf('\n'));
    for i = 1:length(uniqsignal)
      mydisp(sprintf('%i/%i\t',ncorrect(i),ntotal(i)));
    end
    mydisp(sprintf('\n'));
    for i = 1:length(uniqsignal)
      mydisp(sprintf('%0.3f\t',ncorrect(i)/ntotal(i)));
    end
    mydisp(sprintf('\n'));
  end
end

% initial parameter guess
initparams = [median(signal) 3 0];

if isnan(gamma),initparams(end+1) = 0.5;end
  
% set optimization parameters
maxiter = maxIter;
optimparam = optimset('MaxIter',maxiter,'MaxFunEvals',inf);

% range of valid lambda (this is the rate of lapses)
% it controls the maximum height of the psychometric
% curve. 
maxlambda = 0.06;

% use nelder-mead to find maximum-likelihood
[bestfit.fitparams bestfit.fval bestfit.exitflag] = fminsearch(@weibulllike,initparams,optimparam,signal,ncorrect,ntotal,gamma,maxlambda);

% if gamma was not a fit parameter
if ~isnan(gamma)
  bestfit.fitparams(end+1) = gamma;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goodness of fit estimation with monte-carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dogoodnessoffit)
  % get the deviance of this parameter estimate
  D = weibulldeviance(bestfit.fitparams,signal,ncorrect,ntotal);
  
  % now simulate B draws from the model and get a distribution
  % of deviances (monte-carlo)
  if (dogoodnessoffit == 1)
    B = 10000;
  else
    B = dogoodnessoffit;
  end
  disppercent(-inf,'Calculating monte-carlo of deviances');
  for i = 1:B
    disppercent(i/B);
    % get the probabilities for the signal strength given the model
    pmodel = weibull(signal,bestfit.fitparams);
    for j=1:length(signal)
      nmodelcorrect(j) = sum(rand(1,ntotal(j))<=pmodel(j));
    end
    % get the deviance
    Ddist(i) = weibulldeviance(bestfit.fitparams,signal,nmodelcorrect,ntotal);
  end
  disppercent(inf);

  % get where we are in the probability distribution
  bestfit.D = D;
  bestfit.p = sum(D>sort(Ddist))/B;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do bootstrap of parameter estimates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dobootstrap)
  if (dobootstrap == 1)
    B = 1000;
  else
    B = dobootstrap;
  end
  % get probability correct
  pcorrect = ncorrect./ntotal;
  % calculate parameter estimates with B number of resamples
  disppercent(-inf,'Calculating bootstrap');
  bootmax = 1000;
  optimparam = optimset(optimparam,'Display','off','MaxIter',bootmax,'MaxFunEvals',bootmax);
  switch bootmode
    case 'binned'
      parfor i = 1:B
        %disppercent(i/B);
        % resample the distribution of correct answers
        nbootcorrect = [];
        for j=1:length(ntotal)
          nbootcorrect(j) = sum(rand(1,ntotal(j))<pcorrect(j));
        end
        % use nelder-mead to find maximum-likelihood of this resample
        bootparams(i,:) = fminsearch(@weibulllike,initparams,optimparam,signal,nbootcorrect,ntotal,gamma,maxlambda);
      end
    case 'full'
      % full bootstrap resample
      % need to map all the data out to single observations for resampling
      fullcorrect = [];
      fullsignal = [];
      for n = 1:numel(ncorrect)
        thisn = ntotal(n);
        thisc = zeros(1,ntotal(n));
        thisc(1:ncorrect(n)) = 1;
        fullcorrect = [fullcorrect thisc];
        fullsignal = [fullsignal repmat(signal(n),[1 thisn])];
      end
      nfull = numel(fullcorrect);
      parfor i = 1:B
        % sample with replacement
        randsamp = ceil(rand(1,nfull)*nfull);
        bootcorrect = fullcorrect(randsamp);
        bootsignal = fullsignal(randsamp);
        boottotal = ones(1,nfull);
        bootparams(i,:) = fminsearch(@weibulllike,initparams,optimparam,bootsignal,bootcorrect,boottotal,gamma,maxlambda);
      end
    otherwise
      error('unknown bootmode: %s',bootmode)
  end
  disppercent(inf);
  % sort the bootstrap parameters
  for i = 1:size(bootparams,2)
    bootparams(:,i) = sort(bootparams(:,i));
  end
  bestfit.bootparams = bootparams;
  % keep upper and lower bounds of bootstrap estimates
  bestfit.upper = prctile(bootparams,97.5,1);
  bestfit.lower = prctile(bootparams,2.5,1);
  if ~isnan(gamma)
    bestfit.upper(end+1) = gamma;
    bestfit.lower(end+1) = gamma;
  end

  % and calculate 95% confidence interval as the mean between the two
  bestfit.confint = ((bestfit.upper-bestfit.fitparams)+(bestfit.fitparams-bestfit.lower))/2;
end

bestfit.x = min(signal):(max(signal)-min(signal))/100:max(signal);
bestfit.y = weibull(bestfit.x,bestfit.fitparams);
bestfit.signal = signal;
bestfit.pcorrect = ncorrect./ntotal;
bestfit.pcorrectste = sqrt((bestfit.pcorrect.*(1-bestfit.pcorrect))./ntotal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display some stuff
if (dispfig)
  smartfig('fitweibull');
  nrows = 1+(dogoodnessoffit>0);
  % first plot fit
  subplot(nrows,1,1);
  % display fit
  x = min(signal):(max(signal)-min(signal))/100:max(signal);
  fit = weibull(x,bestfit.fitparams);
  hold on
  plot(x,fit,'k-');
  
  % plot data
  myerrorbar(bestfit.signal,bestfit.pcorrect,'yError',bestfit.pcorrectste,'Symbol','o','MarkerFaceColor','k','MarkerSize',8);
  xlabel('signal strength');
  ylabel('proportion correct');
  title(sprintf('T=%0.3f beta=%0.3f gamma=%0.3f lambda = %0.3f',bestfit.fitparams(1),bestfit.fitparams(2),gamma,bestfit.fitparams(3)));
  yaxis(min(0.4,min(bestfit.pcorrect-bestfit.pcorrectste)-0.1),max(1.1,max(bestfit.pcorrect+bestfit.pcorrectste)+0.1));
  vline(bestfit.fitparams(1),'b-');
  hline(gamma);
  hline(1-bestfit.fitparams(3));
  if (dogoodnessoffit)
    % plot distribution of deviances
    subplot(nrows,1,2);
    hold on
    % get the chi-squared approximation
    k = length(ntotal);
    x = min(Ddist):(max(Ddist)-min(Ddist))/100:max(Ddist);
    plot(x,B*chi2pdf(x,k)/sum(chi2pdf(x,k)),'LineWidth',2); 
    legend(sprintf('Chi-squared distribution (deg-of-freedom = %i)',k));
    hist(Ddist,100);
    vline(D,'r-');
    xlabel('deviance');
    ylabel('n');
    if (bestfit.p > 0.975)
      title(sprintf('Bad fit?: D = %0.4f p = %0.4f overdispersion',D,bestfit.p));
    elseif (bestfit.p < 0.025)
      title(sprintf('Fit is too good?: D = %0.4f p = %0.4f underdispersion',D,bestfit.p));
    else
      title(sprintf('Good fit: D = %0.4f p = %0.4f',D,bestfit.p));
    end
  end
  if (dobootstrap)
    smartfig('fitweibull_bootstrap');
    parametername = {'Threshold','Beta','Lambda'};
    xscale = {[min(signal):(max(signal)-min(signal))/1000:max(signal)],[0:0.1:3.5],[0:0.0001:maxlambda]};
    for i = 1:3
      subplot(3,1,i);
      thisdist = bootparams((bootparams(:,i) < max(xscale{i})) & (bootparams(:,i)>min(xscale{i})),i);
      hist(thisdist,xscale{i});
      hold on
      vline(bestfit.fitparams(i),'r-');
      vline(bestfit.upper(i),'g-');
      vline(bestfit.lower(i),'g-');
      xlabel(parametername{i});
      ylabel('n');
      title(sprintf('%s %0.4f <- %0.4f -> %0.4f (%0.2f missing)',parametername{i},bestfit.lower(1),bestfit.fitparams(1),bestfit.upper(1),1-length(thisdist)/size(bootparams,1)));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns likelihood of data given weibull of
% passed in parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function like = weibulllike(fitparams,signal,ncorrect,ntotal,gamma,maxlambda)

% add gamma to parameters if it is not set to be fit
if ~isnan(gamma),fitparams(end+1) = gamma;end

like = 0;
% get the likelihood
for i = 1:length(signal)
  % the probability of correct given the paramters for this
  % signal strength is:
  p = weibull(signal(i),fitparams);
  % the likelihood that we would observe ncorrect trials out
  % of ntotal trials with a probability of p, is assumed
  % to be a bernoulli process and so is binomailly distributed
  % i.e.
  % bernoulli = nCk * p^k * (1-p)^(n-k)
  % taking the log:
  % log bernoulli = log(n/k) + k*log(p) + (n-k)*log(1-p)
  warning off
  like = like + log(nchoosek(ntotal(i),ncorrect(i))) + ncorrect(i)*log(p) + (ntotal(i)-ncorrect(i))*log(1-p);
  warning on backtrace
end

% get the prior for the lambda, it can't be less
% than 0 or greater than maxlambda
if (fitparams(3) < 0) | (fitparams(3) > maxlambda)
  logw = -inf;
else
  logw = 0;
end

% and add it as a prior to the likelihood,
% also invert sign, since we want to maximize the likelihood
% and we are doing a min search
like = -(like+logw);

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns deviance of model given data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = weibulldeviance(params,signal,ncorrect,ntotal)

% calculate the deviance of this parameter estimate
% deviance is the log-likelihood ratio of model with no error
% to this model chosen
D = 0;
for i = 1:length(signal)
  % get the probability correct for the model
  pi = weibull(signal(i),params);
  pcorrect(i) = ncorrect(i)/ntotal(i);
  % and compute the likelihood ratio of the two
  if (ncorrect(i)>0)
    D = D + 2*(ncorrect(i)*log(pcorrect(i)/pi));
  end
  if (ncorrect(i)<ntotal(i))
    D = D + 2*((ntotal(i)-ncorrect(i))*log((1-pcorrect(i))/(1-pi)));
  end
end

%%%%%%%%%%%%%%%%%%%
%    parseArgs    %
%%%%%%%%%%%%%%%%%%%
function [signal ncorrect ntotal gamma dispfig dogoodnessoffit dobootstrap maxIter bootmode] = parseArgs(args)

% set defaults
signal = [];
ncorrect = [];
ntotal = [];
dispfig = 0;
dogoodnessoffit = 0;
dobootstrap = 0;
gamma = 0.5;
maxIter = 10000;
bootmode = 'binned';

% check arguments
if length(args) < 2 help fitweibull; end

% get two needed arguments
signal = args{1};
ncorrect = args{2};

% extra arguments
if length(args) > 2 
  % old style
  if ~isstr(args{3})
    % get number of arguments
    numArgs = length(args);
    % check arguments
    if (numArgs >= 3) ntotal = args{3};end
    if (numArgs >= 4) gamma = args{4};end
    if (numArgs >= 5) dispfig = args{5};end
    if (numArgs >= 6) dogoodnessoffit = args{6};end
    if (numArgs >= 7) dobootstrap = args{7};end
    if (numArgs >= 8) bootmode = args{8};end
    if (numArgs > 8) help fitweibull; end
    return
    % new style
  else
    getArgs({args{3:end}},{'ntotal',ntotal,'gamma',gamma,'dispfig',dispfig,'dogoodnessoffit',dogoodnessoffit,'dobootstrap',dobootstrap,'maxIter',maxIter,'bootmode',bootmode});
    return
  end
end

