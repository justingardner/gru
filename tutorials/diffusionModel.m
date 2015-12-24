% diffusionModel.m
%
%        $Id:$ 
%      usage: diffusionModel()
%         by: justin gardner
%       date: 05/07/15
%    purpose: 
%
function diffusionModel(varargin)

% get arguments for model
getArgs(varargin,{'upperBound=1','lowerBound=-1','startingPoint=0','startingPointSTD=0','driftRate=0.002','driftSTD=0.1','n=10000','slowMethod=0','stopDrawingTrialsAfter=200','updateHistogramsEvery=100','fontName=Helvetica'});


% some initialization parameters
correct = 0;incorrect = 0;
correctReactionTime = []; incorrectReactionTime = [];
correctTrials = [];
incorrectTrials = [];
boundLength = 10;

% initalize the graph
mlrSmartfig('diffusionModel','reuse');clf;hold on
plot([0 boundLength],[upperBound upperBound],'k-','LineWidth',3);
plot([0 boundLength],[lowerBound lowerBound],'k-','LineWidth',3);
axis([0 boundLength lowerBound-(upperBound-lowerBound)*0.5 upperBound+(upperBound-lowerBound)*0.5]);
h = [];
xlabel('RT (samples)','fontName',fontName,'fontAngle','oblique','fontSize',14);
ylabel('Decision variable','fontName',fontName,'fontAngle','oblique','fontSize',14);
title(sprintf('Diffusion model: startingPoint=%s startingPointSTD=%s driftRate=%s driftSTD=%s n=%i',mlrnum2str(startingPoint,'sigfigs=-1'),mlrnum2str(startingPointSTD,'sigfigs=-1'),mlrnum2str(driftRate,'sigfigs=-1'),mlrnum2str(driftSTD,'sigfigs=-1'),n),'fontName',fontName,'fontAngle','oblique','fontSize',14);
maxlen = 1;
d = nan(n,boundLength);

% simulate n runs of the model
disppercent(-inf,'(diffusionModel) Simulating diffusion model');
for i = 1:n
  % calculate a path, set d to the starting point
  d(i,1) = random('norm',startingPoint,startingPointSTD);
  len = 1;
  % we could do the following, which is conceptually simpler - keep drawing
  % from a normal distribution until we hit one of the bounds, but it is a bit
  % slow since we are in a while loop
  if slowMethod
    while (d(i,len) < upperBound) && (d(i,len) > lowerBound)
      d(i,len+1) = d(i,len)+random('norm',driftRate,driftSTD);
      len = len+1;
    end
  else
    % faster method: here we guess how many samples we will need to cross the bound.
    % our best guess is basically how long it took the longest one so far to cross
    % the bound, which is in the parameter boundLength
    while (1)
      randomDrawsFromDrift = random('norm',driftRate,driftSTD,1,boundLength);
      % add these cumulatively to the drift process
      d(i,len:len+length(randomDrawsFromDrift)) = cumsum([d(i,len) randomDrawsFromDrift]);
      % see if we have crossed either bound
      boundCrossing = (d(i,:)>upperBound) | (d(i,:)<lowerBound);
      if any(boundCrossing)
	% the sequence terminated at this length then
	len = min(find(boundCrossing));
	% and break
	break;
      else
	len = len+length(randomDrawsFromDrift);
      end
    end
  end
  % fill out the rest of the sequence to nan to indicate the process has terminated
  d(i,len+1:end) = nan;
  % if we are adding to the length of the whole array then fill everyone else out to nan
  if len>maxlen
    d(1:i-1,maxlen+1:end) = nan;
    maxlen = len;
  end
  % get reaction time of trial
  rt(i) = len;
  % get correct or incorrect
  if d(i,len) > upperBound
    correct=correct+1;
    correctReactionTime(end+1) = rt(i);
    correctTrials(end+1) = i;
    plotColor = 'k:';
  else
    incorrect=incorrect+1;
    incorrectReactionTime(end+1) = rt(i);
    incorrectTrials(end+1) = i;
    plotColor = 'r:';
  end
  % draw the paths of trials for the first few (stop afterwords to make the sim go faster)
  if i < stopDrawingTrialsAfter
    plot(d(i,1:len),plotColor);
    % draw the bounds out to where we are
    if maxlen > boundLength
      plot([boundLength maxlen],[upperBound upperBound],'k-','LineWidth',3);
      plot([boundLength maxlen],[lowerBound lowerBound],'k-','LineWidth',3);
      boundLength = maxlen;
      axis([0 boundLength lowerBound-(upperBound-lowerBound)*0.5 upperBound+(upperBound-lowerBound)*0.5]);
    end
    drawnow;
  end
  % update histograms every few trials
  if mod(i,updateHistogramsEvery) == 0
    h = drawRTHistograms(correct,correctReactionTime,incorrect,incorrectReactionTime,upperBound,lowerBound,i,boundLength,h);
    drawnow;
  end
  disppercent(i/n);
end
disppercent(inf);

% compute mean trials and plot
correctTrials = d(correctTrials,1:maxlen);
correctTrials(isnan(correctTrials)) = upperBound;
plot(mean(correctTrials),'k-','LineWidth',4);
incorrectTrials = d(incorrectTrials,1:maxlen);
incorrectTrials(isnan(incorrectTrials)) = lowerBound;
plot(mean(incorrectTrials),'r-','LineWidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    drawRTHistograms    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = drawRTHistograms(correct,correctReactionTime,incorrect,incorrectReactionTime,upperBound,lowerBound,n,boundLength,h)

% don't do anything if we don't have some counts
if (correct==0) || (incorrect==0),return,end

% compute histogram of correct and incorrect counts
maxRT = max(max(correctReactionTime),max(incorrectReactionTime));
edgeWidth = maxRT/100;
edges = 0:edgeWidth:max(maxRT);
correctCounts = histc(correctReactionTime,edges);
incorrectCounts = histc(incorrectReactionTime,edges);

% initializes handles
if ~isfield(h,'correct')
  h.correct = zeros(1,length(edges));
  h.incorrect = zeros(1,length(edges));
end

% now plots those at the appropriate point on the graph
histMaxHeight = (upperBound-lowerBound)/3;
maxCounts = max(max(correctCounts),max(incorrectCounts));
for i = 1:length(edges)
  % plot correct histogram
  if correctCounts(i) > 0
    % get size of rectangle to draw
    rectanglePosition = [edges(i) upperBound edgeWidth histMaxHeight*correctCounts(i)/maxCounts];
    % draw a new one if it has not been draw before
    if h.correct(i)==0
      h.correct(i) = rectangle('Position',rectanglePosition,'FaceColor','w');
    else
      % or just update position
      set(h.correct(i),'Position',rectanglePosition);
      set(h.correct(i),'Visible','on');
    end
  else
    % make it invisible if the height should be 0
    if h.correct(i)~=0,set(h.correct(i),'Visible','off');end
  end

  % plot incorrect histogram
  if incorrectCounts(i) > 0
    incorrectHeight = histMaxHeight*incorrectCounts(i)/maxCounts;
    rectanglePosition = [edges(i) lowerBound-incorrectHeight edgeWidth incorrectHeight];
    % draw a new one if it has not been draw before
    if h.incorrect(i)==0
      h.incorrect(i) = rectangle('Position',rectanglePosition,'FaceColor','r');
    else
      % or just update position
      set(h.incorrect(i),'Position',rectanglePosition);
      set(h.incorrect(i),'Visible','on');
    end
  else
    % make it invisible if the height should be 0
    if h.incorrect(i)~=0,set(h.incorrect(i),'Visible','off');end
  end
end

% display percent corect and incorrect
correctStr = sprintf('Correct: %0.2f%%, median RT: %0.2f',100*correct/n,median(correctReactionTime));
incorrectStr = sprintf('Incorrect: %0.2f%%, median RT: %0.2f',100*incorrect/n,median(incorrectReactionTime));
if ~isfield(h,'text')
<<<<<<< .mine
  h.text(1) = text(boundLength,upperBound+histMaxHeight,correctStr,'HorizontalAlignment','right','VerticalAlignment','top','FontName',fontName,'FontAngle','oblique','FontSize',12);
  h.text(2) = text(boundLength,lowerBound-histMaxHeight,incorrectStr,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName',fontName,'FontAngle','oblique','FontSize',12);
=======
  h.text(1) = text(boundLength,upperBound+histMaxHeight,correctStr,'HorizontalAlignment','right','VerticalAlignment','top','FontName','Helvetica','FontAngle','oblique','FontSize',12,'Color','k');
  h.text(2) = text(boundLength,lowerBound-histMaxHeight,incorrectStr,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName','Helvetica','FontAngle','oblique','FontSize',12,'Color','r');
>>>>>>> .r465
else
  % just set text string if already drawn
  set(h.text(1),'String',correctStr);
  set(h.text(2),'String',incorrectStr);
  % and set there y position
  pos = get(h.text(1),'Position');pos(1) = boundLength;set(h.text(1),'Position',pos);
  pos = get(h.text(2),'Position');pos(1) = boundLength;set(h.text(2),'Position',pos);
  
end
%drawPublishAxis


%keyboard
  



