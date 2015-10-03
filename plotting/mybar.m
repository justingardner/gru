% mybar.m
%
%        $Id: mybar.m,v 1.3 2008/10/03 12:00:31 justin Exp $ 
%      usage: h = mybar(values,<varargin>)
%         by: justin gardner
%       date: 09/30/08
%    purpose: displays a bar graph. Values should be a nxk array where
%             n is the number of groups and k is the number of values
%             within each group.
%
%             h is a return structure of handles to the objects drawn
%
%             Options:
% 
%             groupLabels: Names of the groups
%               a cell array of n strings
%             withinGroupLabels: Names of values within each group
%               a cell array of k strings
%             type: Only grouped for now
%             yAxisMin: The minimum value for the y axis. 
%             yAxisMax: The maximum value for the y axis.
%             yError: nxk matrix for the error bars on each bar
%             yMin: Is an nxk matrix specifying min for each bar (displayed as white line)
%             yMax: Is an nxk matrix specifying max for each bar (displayed as white line)
%             dispValues: 0 or 1 for whether to display the values in text
%             textValueRot: Rotation of text value (defaults to 90 - set to 0 if you want horizontal)
%             hline: an array of values for where you want to have horizontal lines
%             hLineStyle: Line style for hLine, defaults to :
%             hLineColor: Line color for hLine, defaults to k
%             xLabelText: String with a label for the x-axis
%             yLabelText: String with a label for the y-axis
%             withinGroupColors: Specify the colors within each group. Should be a cell array
%               of [1x3] color values of length number of within groups
%             groupColors: only available if you have 1 item per group. Should be a cell array
%               of [1x3] color values of length number of groups (see example below).
%             pval: same size as values that contains p-values (for putting up * or **)
%             pvalThreshold: Defaults to [0.05 and 0.01] puts a * and ** depending on which threshold is passed
%             pvalThresholdStrings: Defaults to {'*','**'} specifies strings to display to indicate pvalues. SHould
%               be same length as pvalThreshold
%             supraText: Instead of pval thresholds you can specify a cell array with nGroup elements and
%               a text item for each cell element to be displayed above the bars
%             baseValue: Sets the base value from which bars are drawn (this defaults to 0 or yMinAxis)
%             labelPosCutoff: mybar tries to figure out on its own whether to draw the text values inside the 
%               bar or above the bar. If you want to set the cutoff yourself use this value. A single value
%               specifies the point at which below that value the text will be drawn above rather than inside the bar
%               an array of length two can be used for when bars are also drawn downwards (i.e. negative values or
%               values below baseValue). Then the second value specifies the cutoff spearately for those values
%             roundToInt: rounds all values to integer for display
%
%mybar(0.4+rand(2,3),'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},'yAxisMin=0','yAxisMax=1.5','yError',rand(2,3)*0.1)
%
% With custom colors:
%
%mybar(0.4+rand(2,3),'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},'yAxisMin=0','yAxisMax=1.5','yError',rand(2,3)*0.1,'withinGroupColors',{[0.3 0.7 0.2] [0.4 0.2 0.1] [0.1 0.0 0.0]})
%
% Groups having different colors
%
%mybar([-0.02 0.7 0.32 0.4]','groupColors',{[0 0.1 0],[0 0.3 0],[0 0.4 0.5],[0.7 0.2 0.1]},'yError',[0.1 0.2 0.1 0.3]','groupLabels',{'Group 1','Group 2','Group 3','Group 4'})
%
% Display text above bars
%mybar([[1.2 0.54 0.1];[-0.8 -0.1 -0.7]],'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},'yAxisMin=-1.5','yAxisMax=1.5','yError',rand(2,3)*0.3,'supraText',{{'*','n.s.','**'},{'**','*',''}})
%
% Display p-values as * and **
%mybar([[1.2 0.54 0.1];[-0.8 -0.1 -0.7]],'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},'yAxisMin=-1.5','yAxisMax=1.5','yError',rand(2,3)*0.3,'pval',[[0.1 0.2 0.02];[0.049 0.3 0.001]]);
%mybar([[1.2 0.54 0.1];[-0.8 -0.1 -0.7]],'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},'yAxisMin=-1.5','yAxisMax=1.5','yError',rand(2,3)*0.3,'pval',[[0.1 0.2 0.02];[0.049 0.3 0.001]],'pvalThreshold=[0.1 0.05 0.01]','pvalThresholdStrings',{'(*)','*','**'});
%
%
function retval = mybar(vals,varargin)

% check arguments
if any(nargin == [0])
  help mybar
  return
end

% get current axis
ax = gca;
hold on

% get arguments
withinGroupLabels=[];groupLabels=[];type=[];yAxisMin=[];yAxisMax=[];yError=[];dispValues=[];hline=[];xLabelText = [];yLabelText = [];yMin = [];yMax = [];hLineStyle = [];hLineColor = []; colorMap = [];drawAxis = [];roundToInt=[];
getArgs(varargin,{'withinGroupLabels=[]','groupLabels=[]','type=grouped','yAxisMin=[]','yAxisMax=[]','yError=[]','yMin=[]','yMax=[]','dispValues=1','hline=[]','xLabelText','','yLabelText','','hLineStyle=:','hLineColor=k','colorMap=[]','drawAxis=0','roundToInt=0','withinGroupColors=[]','clearAxis=1','groupColors=[]','supraText=[]','pval=[]','pvalThreshold=[0.05 0.01]','pvalThresholdStrings',{'*','**'},'supraTextFontSize',24,'baseValue=[]','labelPosCutoff=[]','textValueRot=90'});

% clear the axis
if clearAxis,cla('reset');end
drawPublishAxis('forceClear=1');

% get some values
nGroups = size(vals,1);
nBarsInGroup = size(vals,2);

% display the hlines
for i = 1:length(hline)
  retval.hLine(i) = line([1 nGroups],[hline(i) hline(i)]);
  set(retval.hLine(i),'Color',hLineColor);
  set(retval.hLine(i),'LineStyle',hLineStyle);
  hold on
end

% check the group colors
if ~isempty(groupColors)
  if nBarsInGroup~=1
    disp(sprintf('(mybar) Number of bars in group must be 1 if you want to set groupColors. Ignoring'));
    groupColors = [];
    % check number of colors
  elseif length(groupColors) ~= nGroups
    disp(sprintf('(mybar) groupColors must have the same number of colors (%i) as groups (%i)',length(groupColors),nGroups));
    groupColors = [];
  end
end

% check colors in withinGroupColors
if ~isempty(withinGroupColors) && (length(withinGroupColors) < nBarsInGroup)
  disp(sprintf('(mybar) withinGroupColors has %i colors, but there are %i groups',length(withinGroupColors),nBarsInGroup));
  withinGroupColors = [];
end

% set baseValue
if isempty(baseValue)
  if yAxisMin > 0
    baseValue = yAxisMin;
  else
    baseValue = 0;
  end
end

  
% display the bars
if nGroups > 1
  % special case when there is only one item per group
  % and the user wants to set groupColors
  if ~isempty(groupColors)
    for i = 1:nGroups
      retval.barHandles(i) = bar(ax,i,vals(i),'ShowBaseLine','off','FaceColor',groupColors{i},'EdgeColor',[0 0 0],'BaseValue',baseValue);
      hold on
    end
  else
    % usual case, group the bars
    retval.barHandles = bar(ax,vals,type,'ShowBaseLine','off','BaseValue',baseValue);
  end
else
  % if we have only group, add a fake second group
  % so that bar does the same thing as when it has multiple groups
  retval.barHandles = bar(ax,[vals;repmat(nan,1,nBarsInGroup)],type,'ShowBaseLine','off','BaseValue',baseValue);
  xaxis(0,1.5);
end

% set group colors
if isempty(groupColors)
  % set the withinGroupColors
  if isempty(withinGroupColors)
    for i = 1:nBarsInGroup
      withinGroupColors{i} = getSmoothColor(i,nBarsInGroup,colorMap);
    end
  end

  % set the bar colors
  for i = 1:nBarsInGroup
    set(retval.barHandles(i),'FaceColor',withinGroupColors{i});
    set(retval.barHandles(i),'EdgeColor',[0 0 0]);
  end

  for iWithinGroup = 1:nBarsInGroup
    groupHandle = get(retval.barHandles(iWithinGroup));
    if ~isempty(groupHandle) && isfield(groupHandle,'Children') && ~isempty(groupHandle.Children)
      groupHandle = get(groupHandle.Children);
      % get x position of bar
      barXpos(iWithinGroup,:) = mean(groupHandle.XData);
    else
      % the above stopped working in Matlab 2015a - this is the new way
      barXpos(iWithinGroup,:) = get(retval.barHandles(iWithinGroup),'XData')+get(retval.barHandles(iWithinGroup),'XOffset');;
    end
  end
else
  barXpos(1,1:nGroups) = 1:nGroups;
end

% display min values
if ~isempty(yMin)
  % display values
  if size(yMin,1) ~= nGroups
    disp(sprintf('(mybar) yMin does not have enough groups (%i, should be %i)',size(yMin,1),nGroups));
  elseif size(yMin,2) ~= nBarsInGroup
    disp(sprintf('(mybar) yMin does not have enough within groups (%i, should be %i)',size(yMin,2),nBarsInGroup));
  else
    for iWithinGroup = 1:nBarsInGroup
      for iGroup = 1:nGroups
	% plot the line
	retval.yMinHandle(iGroup,iWithinGroup) = line([barXpos(iWithinGroup,iGroup) barXpos(iWithinGroup,iGroup)],[vals(iGroup,iWithinGroup) yMin(iGroup,iWithinGroup)],'Color',[0.99 0.99 0.99]);
      end
    end
  end
  % get text position
  yTextPosMin = min(vals,yMin);
else
  yTextPosMin = vals;
end

% stretch the hlines to fit the left and most bar
for i = 1:length(hline)
  set(retval.hLine(i),'XData',[min(barXpos(:)) max(barXpos(:))]);
end

% display max values
if ~isempty(yMax)
  % display values
  if size(yMax,1) ~= nGroups
    disp(sprintf('(mybar) yMax does not have enough groups (%i, should be %i)',size(yMax,1),nGroups));
  elseif size(yMax,2) ~= nBarsInGroup
    disp(sprintf('(mybar) yMax does not have enough within groups (%i, should be %i)',size(yMax,2),nBarsInGroup));
  else
    for iWithinGroup = 1:nBarsInGroup
      for iGroup = 1:nGroups
	% plot the line
	retval.yMaxHandle(iGroup,iWithinGroup) = line([barXpos(iWithinGroup,iGroup) barXpos(iWithinGroup,iGroup)],[vals(iGroup,iWithinGroup) yMax(iGroup,iWithinGroup)],'Color',[0.25 0.25 0.25]);
      end
    end
  end
  % get text position
  yTextPosMax = min(vals,yMax);
else
  yTextPosMax = vals;
end

% display the error bars
if ~isempty(yError)
  if size(yError,1) ~= nGroups
    disp(sprintf('(mybar) yError does not have enough groups (%i, should be %i)',size(yError,1),nGroups));
  elseif size(yError,2) ~= nBarsInGroup
    disp(sprintf('(mybar) yError does not have enough within groups (%i, should be %i)',size(yError,2),nBarsInGroup));
  else
    for iWithinGroup = 1:nBarsInGroup
      for iGroup = 1:nGroups
	% error plots up for positive values, down for negative values
	if vals(iGroup,iWithinGroup) >= baseValue
	  thisYerror = yError(iGroup,iWithinGroup);
	else
	  thisYerror = -yError(iGroup,iWithinGroup);
	end
	% plot the line
	retval.yErrorHandle(iGroup,iWithinGroup) = line([barXpos(iWithinGroup,iGroup) barXpos(iWithinGroup,iGroup)],[vals(iGroup,iWithinGroup) vals(iGroup,iWithinGroup)+thisYerror],'Color',[0 0 0]);
      end
    end
  end
  % get text position
  yTextPosMax = nan(size(vals));
  yTextPosMax(vals>=baseValue) = vals(vals>=baseValue)+yError(vals>=baseValue);
  yTextPosMax(vals<baseValue) = vals(vals<baseValue)-yError(vals<baseValue);
else
  yTextPosMax = vals;
end

axis off

% show the within group labels
if ~isempty(withinGroupLabels)
  % check label lengths
  if length(withinGroupLabels) ~= nBarsInGroup
    disp(sprintf('(mybar) Number of within group labels (%i) is not equal to number of within groups (%i)',length(withinGroupLabels),nBarsInGroup));
  else
    for i = 1:nBarsInGroup
      withinGroupSymbols{i}{1} = 'ks';
      withinGroupSymbols{i}{2} = withinGroupColors{i};
      withinGroupSymbols{i}{3} = 'MarkerFaceColor';
      withinGroupSymbols{i}{4} = withinGroupColors{i};
    end
    retval.hLegend = mylegend(withinGroupLabels,withinGroupSymbols);
  end
end

% get axis
a = axis;
% get yaxis
if isempty(yAxisMin)
  yAxisMin = a(3);
end
if isempty(yAxisMax)
  yAxisMax = a(4);
end
yTickSize = (yAxisMax-yAxisMin)/15;
xAxisMin = a(1);
xAxisMax = a(2);
xTickSize = (xAxisMax-xAxisMin)/15;

yaxis(yAxisMin,yAxisMax);

% display values
if dispValues
  % cutoff to decide when to display a value above the bar
  if isempty(labelPosCutoff)
    labelAboveBarCutoffPos = baseValue +  (yAxisMax-baseValue)*0.1;
    labelAboveBarCutoffNeg = baseValue - (yAxisMax-baseValue)*0.1;
  else
    labelAboveBarCutoffPos = labelPosCutoff(1);
    % if specified as an array of 2 take the last value
    % if not, take the negative of the first value
    if length(labelPosCutoff) > 1
      labelAboveBarCutoffNeg = labelPosCutoff(2);
    else
      labelAboveBarCutoffNeg = -labelPosCutoff(1);
    end
  end      

  % how to align text values
  if textValueRot == 90
    hAlign = {'right','left'};
    vAlign = {'middle','middle'};
  else
    hAlign = {'center','center'};
    vAlign = {'top','bottom'};
  end
  for iWithinGroup = 1:nBarsInGroup
    for iGroup = 1:nGroups
      if roundToInt
	valText = sprintf(' %i',round(vals(iGroup,iWithinGroup)));
      else
	valText = sprintf(' %0.2f',round(100*vals(iGroup,iWithinGroup))/100);
      end
      % display the value
      if vals(iGroup,iWithinGroup) > baseValue
	if yTextPosMin(iGroup,iWithinGroup) > labelAboveBarCutoffPos
	  % display inside the bar
	  retval.hText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMin(iGroup,iWithinGroup),valText,'Color',[0.99 0.99 0.99],'VerticalAlignment',vAlign{1},'Rotation',textValueRot,'HorizontalAlignment',hAlign{1});
	else
	  retval.hText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMax(iGroup,iWithinGroup),valText,'Color',[0 0 0],'VerticalAlignment',vAlign{2},'Rotation',textValueRot,'HorizontalAlignment',hAlign{2});
	  % reset yTextPosMax, to work for supraText
	  textExtent = get(retval.hText(iGroup,iWithinGroup),'Extent');
	  yTextPosMax(iGroup,iWithinGroup) = textExtent(2)+textExtent(4);
	end
      else
	if yTextPosMin(iGroup,iWithinGroup) < labelAboveBarCutoffNeg
	  retval.hText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMin(iGroup,iWithinGroup),valText,'Color',[0.99 0.99 0.99],'VerticalAlignment',vAlign{2},'Rotation',textValueRot,'HorizontalAlignment',hAlign{2});
	  % reset yTextPosMin, to work for supraText
	  textExtent = get(retval.hText(iGroup,iWithinGroup),'Extent');
	  yTextPosMin(iGroup,iWithinGroup) = yTextPosMax(iGroup,iWithinGroup);
	else
	  retval.hText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMax(iGroup,iWithinGroup),valText,'Color',[0 0 0],'VerticalAlignment',vAlign{1},'Rotation',textValueRot,'HorizontalAlignment',hAlign{1});
	  % reset yTextPosMin, to work for supraText
	  textExtent = get(retval.hText(iGroup,iWithinGroup),'Extent');
	  yTextPosMin(iGroup,iWithinGroup) = textExtent(2);
	end
      end
    end
  end
end

% create supra text from pvals
if isempty(supraText) && ~isempty(pval)
  % check a few things
  if ~isequal(size(pval),size(vals))
    disp(sprintf('(mybar) Size of pval and vals must be the same'));
  elseif ~isequal(length(pvalThreshold),length(pvalThresholdStrings))
    disp(sprintf('(mybar) Size of pvalThreshold and pvalThresholdStrings must be the same'));
  else
    % ok actuall do it
    for iGroup = 1:nGroups
      for iWithinGroup = 1:nBarsInGroup
	pvalBelowThreshold = max(find(pval(iGroup,iWithinGroup) <= pvalThreshold));
	if ~isempty(pvalBelowThreshold)
	  supraText{iGroup}{iWithinGroup} = pvalThresholdStrings{pvalBelowThreshold};
	end
      end
    end
  end  
end

% display supra-text
if ~isempty(supraText)
  missingSupraText = false;
  for iGroup = 1:nGroups
    for iWithinGroup = 1:nBarsInGroup
      % get label
      if (length(supraText) >= iGroup) && (length(supraText{iGroup}) >= iWithinGroup)
	thisSupraText = supraText{iGroup}{iWithinGroup};
	if ~isempty(thisSupraText)
	  if vals(iGroup,iWithinGroup) > baseValue
	    retval.hSupraText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMax(iGroup,iWithinGroup),thisSupraText,'Color',[0 0 0],'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',supraTextFontSize);
	  else
	    retval.hSupraText(iGroup,iWithinGroup) = text(barXpos(iWithinGroup,iGroup),yTextPosMin(iGroup,iWithinGroup),thisSupraText,'Color',[0 0 0],'VerticalAlignment','top','HorizontalAlignment','center','FontSize',supraTextFontSize);
	  end
	end
      else
	missingSupraText = true;
      end
    end
  end
  % display if there is a missing supra text value
  if missingSupraText
    disp(sprintf('(mybar) Missing supraText values. Should be a cell array with number of groups of elements and each element having a text string for each within group value'));
  end
end

% check the group labels
if ~isempty(groupLabels)
  % check label lengths
  if length(groupLabels) ~= nGroups
    disp(sprintf('(mybar) Number of group labels (%i) is not equal to number of groups (%i)',length(groupLabels),nGroups));
    groupLabels = [];
  end
end

% parameters for the xaxis
if drawAxis 
  if nGroups > 1
    axisX.x = 1;
    axisX.y = -yTickSize/3;
    axisX.len = nGroups-1;
    axisX.ticklen = yTickSize;
    axisX.tickpos = 2;
    axisX.orient = 0;
    axisX.offset = 0;
    axisX.ticks = suggesttick(1,nGroups,1:nGroups);

    if ~isempty(groupLabels)
      axisX.ticks.labels = groupLabels;
    end
  
    % draw the x axis
    yLim = get(gca,'YLim');yLim(1) = axisX.y-axisX.ticklen;
    set(gca,'YLim',yLim);
    drawaxis(axisX);
  else
    if ~isempty(groupLabels)
      xLabelText = sprintf('%s: %s',xLabelText,groupLabels{1});
    end
  end

  % draw the y axis
  axisY.x = 0.5;
  axisY.y = 0;
  axisY.len = yAxisMax-yAxisMin;
  axisY.ticklen = xTickSize;
  axisY.tickpos = 0;
  axisY.orient = 90;
  axisY.offset = 0;
  axisY.ticks = suggesttick(yAxisMin,yAxisMax);

  xLim = get(gca,'XLim');xLim(1) = axisY.y-axisY.ticklen;
  set(gca,'XLim',xLim);
  drawaxis(axisY);

  % set labels
  xlabel(xLabelText);
  ylabel(yLabelText);

else
  set(gca,'XTick',1:nGroups);
  axis on;
  if ~isempty(groupLabels)
    set(gca,'XTickLabel',groupLabels);
  end
  % set labels
  xlabel(xLabelText);
  ylabel(yLabelText);
  drawPublishAxis('forceDisplay=1');
end

% turn the labels on
set(get(gca,'XLabel'),'Visible','on');
set(get(gca,'YLabel'),'Visible','on');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontName','Helvetica');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontName','Helvetica');
set(get(gca,'Title'),'FontName','Helcetica');
set(get(gca,'Title'),'FontSize',14);
xLabelPosition = get(get(gca,'XLabel'),'Position');
xLabelPosition(1) = median(1:nGroups);
set(get(gca,'XLabel'),'Position',xLabelPosition);

