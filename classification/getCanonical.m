% getCanonical.m
%
%        $Id:$ 
%      usage: canonical = getCanonical(v,roi,stimvol,type,varargin)
%         by: justin gardner
%       date: 06/01/11
%    purpose: Creates a "canonical" response that can be used for doing GLM.
%             Pass in view, roi and stimvol. Type should be set to one of:
%
%             type = 'all'  : uses all voxels to create a canonical response, returns
%                           the deconvolved response computed by averaging all
%                           time series together and treating every trial as the same type.
%             type = 'allfit1' : same as all, but returns a smooth gamma fit
%             type = 'allfit2' : same as all, but returns a smooth difference of gamma fit
%             'r2cutoff=[]'    : if this is set, first computes the r2 captured by timelocked
%                                activity. Then only uses voxels that survive the cutoff
%             'normType=max'   : How to normalize the curve. Set to max (default) for 
%                           dividing by the max value. For events all of the same length
%                           this will scale so that beta weights equal the max amp response
%                           set to 'area' if you want to normalize by the area under the curve
%                           this will make for beta weights that scale with length of stim
%
%             'alwaysPos=1' : Makes the amplitude of the canonical always positive (that way
%                           your beta weights will be negative when you have a negative response)
%             'rejectr2=0' : If the r2 of the fit to the canonical is less than rejectr2 will
%                           retrun as [].
%                                 
%
function [canonicalResponse r2 d] = getCanonical(v,roi,stimvol,type,varargin)

% check arguments
if nargin < 4
  help getCanonical
  return
end

% parse arguments
hdrlen = [];
getArgs(varargin,{'hdrlen=20','r2cutoff=[]','normType=max','pval=[]','alwaysPos=1','rejectr2=0'});

% default return argument
canonicalResponse = [];
r2 = [];
d = [];

% get stimvol that makes every event the same (i.e. we are going
% to look for the average response to all events)
stimvolCanonical{1} = cell2mat(stimvol);

% get info about time series
v = viewSet(v,'curGroup',roi.groupNum);
v = viewSet(v,'curScan',roi.scanNum);
concatInfo = viewGet(v,'concatInfo',roi.scanNum,roi.groupNum);
framePeriod = viewGet(v,'framePeriod',roi.scanNum,roi.groupNum);

% check r2 cutoff (if set then take mean of only those voxels that exceed cutoff)
% otherwise choose all voxels to mean over
if ~isempty(r2cutoff) || (nargout >= 2) || ~isempty(pval)
  % get fit
  d = fitTimecourse(roi.tSeries,stimvolCanonical,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=none','displayFit=0','hdrlen',hdrlen,'verbose=0');
  r2 = d.deconv.r2;
  if ~isempty(r2cutoff)
    % check voxels that survive the r2 cutoff
    survivors = find(d.deconv.r2 >= r2cutoff);
    disp(sprintf('(getCanonical) Found %i/%i (%0.2f%%) of voxels surviving r2 cutoff of %0.3f',length(survivors),d.deconv.dim(3),100*length(survivors)/d.deconv.dim(3),r2cutoff));
    if isempty(survivors),return,end
    tSeries = mean(roi.tSeries(survivors,:),1);
  else
    % get mean time series
    tSeries = mean(roi.tSeries,1);
  end
else
  % get mean time series
  tSeries = mean(roi.tSeries,1);
end

%see if we are supposed to use pval cutoff
if ~isempty(pval)
  % compute pval cutoff (randomize enough for 100,000 total values)
  cutoff = getr2cutoff(roi.tSeries,stimvolCanonical,framePeriod,hdrlen,'concatInfo',concatInfo,'n',round(10000/roi.n),'computeOverAllVoxels',1,'pval',pval);
  survivors = find(r2 >= cutoff.pval);
  disp(sprintf('(getCanonical) Found %i/%i (%0.2f%%) of voxels surviving pval cutoff of %0.3f',length(survivors),d.deconv.dim(3),100*length(survivors)/d.deconv.dim(3),pval));
  if isempty(survivors),return,end
  tSeries = mean(roi.tSeries(survivors,:),1);
end

switch (type)
  case {'allfit2'}
   % use fitTimecourse to get the fit deconvolution for all stimulus types
   d = fitTimecourse(tSeries,stimvolCanonical,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=fit2','hdrlen',hdrlen,'minTimelag2',5,'maxTimelag',4);
   % resample the output for each volume aquisition
   canonicalResponse = interp1(d.time,d.ehdr,(0:framePeriod:hdrlen)+framePeriod/2);
 case {'allfit1'}
   % use fitTimecourse to get the fit deconvolution for all stimulus types
   d = fitTimecourse(tSeries,stimvolCanonical,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=fit1','hdrlen',hdrlen,'minTimelag',1.5);
   % resample the output for each volume aquisition
   canonicalResponse = interp1(d.time,d.ehdr,(0:framePeriod:hdrlen)+framePeriod/2);
 case {'all'}
   d = fitTimecourse(tSeries,stimvolCanonical,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=fit1','displayFit=0','hdrlen',hdrlen);
   canonicalResponse = d.deconv.ehdr;
 otherwise
  disp(sprintf('(getCanonical) Unknown type for getting canonical'));
  return
end   

% amplitude should always be positive
if (d.amplitude <= 0) && alwaysPos
  disp(sprintf('(getCanonical) canonicalResponse amplitude (%0.2f) is less than 0. Making positive. If you do not want this behavior than set the input parmaeter alwaysPos to 0',d.amplitude));
  canonicalResponse = -canonicalResponse;
end

% reject when r2 is too low
if ~isempty(rejectr2) & (d.deconv.r2 < rejectr2)
  disp(sprintf('(getCanonical) REJECT of canonicalResponse because fit r2 (%0.2f) is less than criteria (%0.2f). Set rejectr2 parameter to something smaller if you would like to keep.',d.deconv.r2,rejectr2));
  canonicalResponse = [];
  return
end
  
% normalize area under curve
if strcmp(normType,'area')
  % Normalize the mean response. This is so that the amplitude measures
  % come out meaning something. We take the area under the curve to be equal
  % to 1% signal change. This means that the amplitude measure is in %signal change/sec
  canonicalResponse = canonicalResponse/(sum(canonicalResponse)*length(canonicalResponse)*framePeriod);
elseif strcmp(normType,'max');
  % this makes the beta equal to the peak (but only for non-temporally summed responses)
  canonicalResponse = canonicalResponse/max(canonicalResponse);
else
  disp(sprintf('(getCanonical) Unknown normType: %s',normType));
end


