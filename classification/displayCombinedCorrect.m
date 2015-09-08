% displayCombinedCorrect.m
%
%        $Id: displayCombinedCorrect.m,v 1.4 2009/02/09 02:25:15 justin Exp $ 
%      usage: displayCombinedCorrect(v,combined,combinedTitle,saveEPS)
%         by: justin gardner
%       date: 10/01/08
%    purpose: 
%
function displayCombinedCorrect(v,combined,combinedTitle,saveEPS)

% check arguments
if ~any(nargin == [2 3 4])
  help displayCombinedCorrect
  return
end

% get input arguments
if ieNotDefined('combinedTitle')
  combinedTitle = sprintf('%s',getLastDir(viewGet(v,'homedir')));
  combinedSaveName = sprintf('%s_combined',getLastDir(viewGet(v,'homedir')));
else
  combinedSaveName = fixBadChars(combinedTitle);
  combinedTitle = sprintf('%s: %s',getLastDir(viewGet(v,'homedir')),combinedTitle);
end  
if ieNotDefined('saveEPS'),saveEPS = 0;end

% see if this is a variable that has an 'all' fieldd
if isfield(combined,'contraAll')
  % get the fields
  contra = sort(combined.contraAll);
  ipsi = sort(combined.ipsiAll);
  n = size(contra,1);
  % parse into contra/ipsi
  combined.contra = mean(contra);
  combined.contraSTE = std(contra)/sqrt(n);
  combined.ipsi = mean(ipsi);
  combined.ipsiSTE = std(ipsi)/sqrt(n);
  % compute confidence intervals
  ci = 0.95;
  minN = max(round(n*((1-ci)/2)),1);maxN = min(round(n*(1-((1-ci)/2))),100);
  yMin = 100*[ipsi(minN,:);contra(minN,:)]';
  yMax = 100*[ipsi(maxN,:);contra(maxN,:)]';
  % and change title
  combinedTitle = sprintf('%s (CI=%0.2f)',combinedTitle,ci);
else
  yMin=[];yMax=[];
end

% display
smartfig(combinedTitle);clf;

% combine the values for ipsi and contra
vals =100*[combined.ipsi;combined.contra]';
valErrors = 100*[combined.ipsiSTE;combined.contraSTE]';
names = combined.name;

chanceLevel = 1/length(combined.nLeft);
mybar(vals,'groupLabels',names,'yAxisMin=0','yAxisMax=100','yError',valErrors,'withinGroupLabels',{'ipsi','contra'},'xLabelText=Area','yLabelText=Percent correct',sprintf('hline=%i',100*chanceLevel),'yMin',yMin,'yMax',yMax);

if ~isfield(combined,'nRight')
  title(sprintf('%s',combinedTitle),'Interpreter','none');
else
  % string for describing parameters
  if isfield(combined,'groupTrials') && isfield(combined,'startLag') && isfield(combined,'blockLen')
    paramsString = sprintf('\ngroupTrials: %i startLag: %i blockLen: %i',combined.groupTrials,combined.startLag,combined.blockLen);
  else
    paramsString = '';
  end
  % string for describing trials / instances
  if isfield(combined,'nLeftInstances') && isfield(combined,'nRightInstances')
    trialnumString = sprintf('\nLeft (%i trials -> %i instances) Right (%i trials -> %i instances)',round(sum(combined.nLeft)/length(combined.nLeft)),round(sum(combined.nLeftInstances)/length(combined.nLeftInstances)),round(sum(combined.nRight)/length(combined.nRight)),round(sum(combined.nRightInstances)/length(combined.nRightInstances)));
  else
    trialnumString = sprintf('\nLeft (%i trials) Right (%i trials)',round(sum(combined.nLeft)/length(combined.nLeft)),round(sum(combined.nRight)/length(combined.nRight)));
  end
  % put up title
  title(sprintf('%s%s%s',combinedTitle,paramsString,trialnumString),'Interpreter','none');
end

drawnow;
if saveEPS==1
  eval(sprintf('print -depsc2 Anal/%s.eps',combinedSaveName));
elseif saveEPS==2
  eval(sprintf('print -dpng Anal/%s.png',combinedSaveName));
end  