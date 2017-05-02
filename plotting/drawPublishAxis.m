% drawPublishAxis.m
%
%        $Id:$ 
%      usage: drawPublishAxis()
%         by: justin gardner
%       date: 06/10/13
%    purpose: draws publishable axis. Note that these are not matlab axis, but lines with text drawn on the figure. To revert back to the matlab axis, call drawPublishAxis again with no arguments
%    
%             There are many settings you can change:
%
%             whichAxis (both, horizontal, vertical): which axis to draw
%             lineWidth (1): Line width of the axis
%             titleStr: title to display
%             labelFontSize (12): font size for labels
%
%             xAxisLoc (bottom, top): whether to draw the x axis on the top or bottom
%             xAxisOffset (-1/32): How much to offset the x axis. This is specified as the 
%                amount of the full range of the y axis to offset.
%             xAxisYpos: Overrides the above and sets to the exact position.
%             xAxisMajorTickLen (-1/32): Length of major tick on x axis, if negative is a percent of the
%                full size of the y axis. If positive, actual value that will be added to xAxisYpos.
%             xAxisMinorTickLen (-1/48): Like xAxisMajorTickLen but for the minor ticks
%             tickDir (in, out): Specifies direction of tick marks
%             xTick: An array that specifies the ticks (like the property of the axis);
%             xTickLabel: A cell array of labels for the ticks. Should be the same length as xTick. Also
%                 can be a string array like what is returned from the axis property. Or can be a regular
%                 array in which case mlrnum2str is used to convert values.
%             xTickLabelSigfigs (-1): Number of sigfigs to display xTickLabel with (only used in case xTickLabel
%                 is an array).
%             xAxisMin: Sets the min value of the x axis
%             xAxisMax: Sets the max value of the x axis
%             xAxisMinMaxSetByTicks (1): If 1 sets xAxisMin/Max to the min max of the ticks
%             xMinorTick: Array specifying the minor ticks. If it is of length one, it specifies
%                 spacing of minor ticks. If it is a log scale and the length is one should be an 
%                 integer specifying how many tick values you want between each major tick
%             xLabel: String that specifies xLabel
%             xLabelOffset (-3/64): Offset from x axis to display label
%             xAxisMargin (1/64): Extra little offset to increase limits of axis by to make
%                 sure everything displays correctly
%             figSize (0,0.5,1,2): Sets units of figure apropriately for exporting use print -dpdf. 0.5, 1, 2 will
%                 give half, one or two column scaling. Or as a vector of 2 gives the x and y dimensions in cm
%                 of final figure
%
%             Similar arguments are available for the vertical axis (preface with y
%                  instead of x)
%
function retval = drawPublishAxis(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get general arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get arguments
getArgs(varargin,{'whichAxis=both','tickDir=out','lineWidth=1','titleStr=[]','labelFontSize=12',...
		  'xMin=[]','xMax=[]','xAxisLoc=[]','xAxisYpos=[]','xAxisMargin=1/64',...
		  'xAxisOffset=-1/32','xScale=[]','xAxisMajorTickLen=-1/32','xAxisMinorTickLen=-1/48',...
		  'xTick=[]','xTickLabel=[]','xTickLabelHide',false,'xAxisMin=[]','xAxisMax=[]','xAxisMinMaxSetByTicks=1',...
		  'xMinorTick=[]','xLabel=[]','xLabelOffset=-4/64','xTickLabelSigfigs=-1',...
		  'yMin=[]','yMax=[]','yAxisLoc=[]','yAxisXpos=[]','yAxisMargin=1/64',...
		  'yAxisOffset=-1/32','yScale=[]','yAxisMajorTickLen=-1/32','yAxisMinorTickLen=-1/48',...
		  'yTick=[]','yTickLabel=[]','yAxisMin=[]','yAxisMax=[]','yAxisMinMaxSetByTicks=1',...
		  'yMinorTick=[]','yLabel=[]','yLabelOffset=-6/64','yTickLabelSigfigs=-1','forceDisplay=0',...
		  'forceClear=0','fontName=Helvetica','figSize=2'...
		 });
     
% set fig size
if figSize
  % borrowed from savepdf
  setFigSize(figSize);
end

% get which axis to draw
if ~validateParam('whichAxis',{'both','horizontal','vertical'});,return,end

% get which direction the ticks go in
if ~validateParam('tickDir',{'in','out'});,return,end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revert axis if this function was already called
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curAxis = gca;

% see if we have run before on this axis
dpa = get(curAxis,'UserData');
if isfield(dpa,'drawPublishAxis')
  % ist all the fields that have drawable handles
  drawableFields = {'xAxis','yAxis','axis'};
  for iDrawableFields = 1:length(drawableFields)
    drawableField = drawableFields{iDrawableFields};
    % remove all drawables
    if isfield(dpa,drawableField)
      % go through and remove everything drawn
      xAxisObjs = fieldnames(dpa.(drawableField));
      for iObj = 1:length(xAxisObjs)
	thisObj = dpa.(drawableField).(xAxisObjs{iObj});
	for j = 1:length(thisObj)
	  if ishandle(thisObj(j))
	    delete(thisObj(j));
	  end
	end
      end
      % remove the xaxis elements we have deleted
      dpa = rmfield(dpa,drawableField);
    end
  end
  % reset axis size
  if isfield(dpa,'originalAxisSize')
    axis(dpa.originalAxisSize);
  end

  % bring back title
  if isfield(dpa,'titleVisible') && dpa.titleVisible
    set(dpa.titleHandle,'Visible','on');
  end

  % reset on off
  if isfield(dpa,'axisVisible')
    set(curAxis,'Visible',dpa.axisVisible);
  end

  % turn on/off xLabel
  if isfield(dpa,'xLabelVisible') && dpa.xLabelVisible
    set(get(curAxis,'XLabel'),'Visible','on');
  else
    set(get(curAxis,'XLabel'),'Visible','off');
  end
  
  % turn on/off xLabel
  if isfield(dpa,'yLabelVisible') && dpa.yLabelVisible
    set(get(curAxis,'YLabel'),'Visible','on');
  else
    set(get(curAxis,'YLabel'),'Visible','off');
  end
  
  % set the UserData field to empty
  set(curAxis,'UserData',[]);
  
  % if passed with no arguments and there was already some axis
  % then just return
  if nargin == 0,return,end
end
if forceClear,return,end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get information about axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the title
dpa.titleHandle = get(curAxis,'Title');
if ishandle(dpa.titleHandle)
  % grab the title
  if isempty(titleStr)
    titleStr = get(dpa.titleHandle,'String');
  end
  % get whether title is currently visible
  if get(dpa.titleHandle,'Visible')
    dpa.titleVisible = true;
  end
end

% get yTick
if isempty(yTick) 
  yTick = get(curAxis,'YTick');
else
  % set the yTick appopriately
  yaxis(min(yTick),max(yTick));
  set(curAxis,'YTick',yTick);
end
yTick = yTick(:)';

% get xTick
if isempty(xTick) 
  % if not set, get from current axis
  xTick = get(curAxis,'XTick');
else
  % set the axis approriately
  xaxis(min(xTick),max(xTick));
  set(gca,'XTick',xTick);
end
xTick = xTick(:)';

% get x axis limits
xLim = get(curAxis,'XLim');
if isempty(xMin), xMin = xLim(1);end
if isempty(xMax), xMax = xLim(2);end

% get y axis limits
yLim = get(curAxis,'YLim');
if isempty(yMin), yMin = yLim(1);end
if isempty(yMax), yMax = yLim(2);end

% get direction of x-axis
xDir = get(curAxis,'xDir');
if ~validateParam('xDir',{'normal','reverse'});,return,end

% get location of x-axis
if isempty(xAxisLoc), xAxisLoc=get(curAxis,'xAxisLoc');end
if ~validateParam('xAxisLoc',{'bottom','top'});,return,end

% get xScale
if isempty(xScale), xScale=get(curAxis,'xScale');end
if ~validateParam('xScale',{'linear','log'});,return,end

% get xTickLabel
if isempty(xTickLabel)
  % get info form axis
  xTickLabelFromAxis = get(curAxis,'XTickLabel');
  xTickFromAxis = get(curAxis,'XTick');
  % try to get labels from axis
  for iTick = 1:length(xTick)
    matchFromAxis = find(xTickFromAxis==xTick(iTick));
    if isempty(matchFromAxis) || strcmp(xScale,'log')
      % use the value of xTick
      xTickLabelCellArray{iTick} = mlrnum2str(xTick(iTick),'sigfigs',xTickLabelSigfigs);
    else
      % if we found it on the axis then use that label
      if ~iscell(xTickLabelFromAxis)
        xTickLabelCellArray{iTick} = xTickLabelFromAxis(matchFromAxis,:);
      else
        xTickLabelCellArray{iTick} = xTickLabelFromAxis{matchFromAxis};
      end
    end
  end
  xTickLabel = xTickLabelCellArray;
end
% if a string then convert to cell array
if isstr(xTickLabel)
  for iLabel = 1:size(xTickLabel,1)
    xTickLabelCellArray{iLabel} = xTickLabel(iLabel,:);
  end
  xTickLabel = xTickLabelCellArray;
% if numeric convert to cell array
elseif isnumeric(xTickLabel)
  for iLabel = 1:length(xTickLabel)
    xTickLabelCellArray{iLabel} = mlrnum2str(xTickLabel(iLabel),'sigfigs',xTickLabelSigfigs);
  end
  xTickLabel = xTickLabelCellArray;
end
% say if something bad happened
if ~iscell(xTickLabel) || (length(xTickLabel) ~= length(xTick))
  disp(sprintf('(drawPublishAxis) xTickLabel should be a cell array of length %i',length(xTick)));
  return
end

% get x label
if isempty(xLabel)
  xLabel = get(curAxis,'XLabel');
  if strcmp(get(xLabel,'Visible'),'on')
    set(xLabel,'Visible','off');
    dpa.xLabelVisible = 1;
    xLabel = get(xLabel,'String');
  else
    dpa.xLabelVisible = 0;
    xLabel = '';
  end
end

% get direction of y-axis
yDir = get(curAxis,'yDir');
if ~validateParam('yDir',{'normal','reverse'});,return,end

% get location of y-axis
if isempty(yAxisLoc), yAxisLoc=get(curAxis,'yAxisLoc');end
if ~validateParam('yAxisLoc',{'left','right'});,return,end

% get yScale
if isempty(yScale), yScale=get(curAxis,'yScale');end
if ~validateParam('yScale',{'linear','log'});,return,end

% get yTickLabel
if isempty(yTickLabel)
  % get info form axis
  yTickLabelFromAxis = get(curAxis,'YTickLabel');
  yTickFromAxis = get(curAxis,'YTick');
  % try to get labels from axis
  for iTick = 1:length(yTick)
    matchFromAxis = find(yTickFromAxis==yTick(iTick));
    if isempty(matchFromAxis) || strcmp(yScale,'log')
      % use the value of xTick
      yTickLabelCellArray{iTick} = mlrnum2str(yTick(iTick),'sigfigs',yTickLabelSigfigs);
    else
      if ~iscell(yTickLabelFromAxis)
	% if we found it on the axis then use that label
	yTickLabelCellArray{iTick} = yTickLabelFromAxis(matchFromAxis,:);
      else
	% grab from cell array
	yTickLabelCellArray{iTick} = yTickLabelFromAxis{matchFromAxis};
      end
    end
  end
  yTickLabel = yTickLabelCellArray;
end
% if a string then convert to cell array
if isstr(yTickLabel)
  for iLabel = 1:size(yTickLabel,1)
    yTickLabelCellArray{iLabel} = yTickLabel(iLabel,:);
  end
  yTickLabel = yTickLabelCellArray;
% if a numeric then convert to cell array
elseif isnumeric(yTickLabel)
  for iLabel = 1:length(yTickLabel)
    yTickLabelCellArray{iLabel} = mlrnum2str(yTickLabel(iLabel),'sigfigs',yTickLabelSigfigs);
  end
  yTickLabel = yTickLabelCellArray;
end
% check if anything went wrong
if ~iscell(yTickLabel) || (length(yTickLabel) ~= length(yTick))
  disp(sprintf('(drawPublishAxis) yTickLabel should be a cell array of length %i',length(yTick)));
  return;
end

% get the y-label
if isempty(yLabel)
  yLabel = get(curAxis,'YLabel');
  if strcmp(get(yLabel,'Visible'),'on')
    set(yLabel,'Visible','off');
    dpa.yLabelVisible = 1;
    yLabel = get(yLabel,'String');
  else
    dpa.yLabelVisible = 0;
    yLabel = '';
  end
end

% this is the beginning of a variable which will be attached to the UserData field
% of the axis to specify what has been done
dpa.drawPublishAxis = true;
dpa.originalAxisSize = axis;
dpa.axisOnOff = 'on';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get positioning information for drawing horizontal axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(whichAxis,{'both','horizontal'}))
  % get the y location of the x-axis and which way ticks should go
  if (strcmp(xAxisLoc,'bottom') && strcmp(yDir,'normal')) || (strcmp(xAxisLoc,'top') && strcmp(yDir,'reverse'))
    % set xAxis y position for bottom
    if isempty(xAxisYpos)
      xAxisYpos = logScaleSum(yScale,yMin,xAxisOffset,abs(yMax-yMin));
    end
    % make tickDir a +/- value
    if strcmp(tickDir,'in') tickDirVal = -1; else tickDirVal = 1; end
  else
    % set xAxis y position for top
    if isempty(xAxisYpos)
      xAxisYpos = logScaleSum(yScale,yMax,-xAxisOffset,abs(yMax-yMin));
    end
    % make tickDir a +/- value
    if strcmp(tickDir,'in') tickDirVal = 1; else tickDirVal = -1; end
  end

  % get the bottom position of the major tick (xAxisTickYpos)
  if (xAxisMajorTickLen < 0)
    xAxisTickYpos = logScaleSum(yScale,xAxisYpos,tickDirVal*xAxisMajorTickLen,abs(yMax-yMin));
  else
    xAxisTickYpos = xAxisYpos + tickDirVal * xAxisMajorTickLen;
  end

  % get the bottom position of the minor tick (xAxisMinorTickYpos)
  if (xAxisMinorTickLen < 0)
    xAxisMinorTickYpos = logScaleSum(yScale,xAxisYpos,tickDirVal*xAxisMinorTickLen,abs(yMax-yMin));
  else
    xAxisMinorTickYpos = xAxisYpos + tickDirVal * xAxisMinorTickLen;
  end

  % get x and y text alignment
  if strcmp(xAxisLoc,'bottom')
    xTextAlignment = 'top';
  else
    xTextAlignment = 'bottom';
  end

  % get the minor ticks
  if length(xMinorTick) == 1
    xMinorTick = logScaleSpace(xScale,xTick,xMinorTick);
  end
  xMinorTick = setdiff(xMinorTick,xTick);
  
  % get the min and max of the axis
  if isempty(xAxisMin) xAxisMin = xMin; end
  if isempty(xAxisMax) xAxisMax = xMax; end
  if xAxisMinMaxSetByTicks && (length(xTick)>1), xAxisMin = min(xTick); xAxisMax = max(xTick); end

  % get the xLabelPosition
  xLabelX = logScaleMean(xScale,[xAxisMin xAxisMax]);
  xLabelY = logScaleSum(yScale,xAxisYpos,tickDirVal*xLabelOffset,abs(yMax-yMin));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get positioning information for drawing vertical axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(whichAxis,{'both','vertical'})) || ~isempty(yLabel)
  % get the y location of the x-axis and which way ticks should go
  if (strcmp(yAxisLoc,'left') && strcmp(xDir,'normal')) || (strcmp(yAxisLoc,'right') && strcmp(xDir,'reverse'))
    % set yAxis x position for left
    if isempty(yAxisXpos)
      yAxisXpos = logScaleSum(xScale,xMin,yAxisOffset,abs(xMax-xMin));
    end
    % make tickDir a +/- value
    if strcmp(tickDir,'in') tickDirVal = -1; else tickDirVal = 1; end
  else
    % set yAxis x position for right
    if isempty(yAxisXpos)
      yAxisXpos = logScaleSum(xScale,xMax,-yAxisOffset,abs(xMax-xMin));
    end
    % make tickDir a +/- value
    if strcmp(tickDir,'in') tickDirVal = 1; else tickDirVal = -1; end
  end

  % get the left position of the major tick (yAxisTickXpos)
  if (yAxisMajorTickLen < 0)
    yAxisTickXpos = logScaleSum(xScale,yAxisXpos,tickDirVal*yAxisMajorTickLen,abs(xMax-xMin));
  else
    yAxisTickXpos = yAxisXpos+tickDirVal*yAxisMajorTickLen;
  end

  % get the left position of the minor tick (yAxisMinorTickXpos)
  if (yAxisMinorTickLen < 0)
    yAxisMinorTickXpos = logScaleSum(xScale,yAxisXpos,tickDirVal*yAxisMinorTickLen,abs(xMax-xMin));
  else
    yAxisMinorTickXpos = yAxisXpos + tickDirVal * yAxisMinorTickLen;
  end
  
  % get x and y text alignment
  if strcmp(yAxisLoc,'right')
    yTextAlignment = 'left';
    yLabelTextAlignment = 'top';
  else
    yTextAlignment = 'right';
    yLabelTextAlignment = 'bottom';
  end
  
  % get the minor ticks
  if length(yMinorTick) == 1
    yMinorTick = logScaleSpace(yScale,yTick,yMinorTick);
  end
  yMinorTick = setdiff(yMinorTick,yTick);
  
  % get the min and max of the axis
  if isempty(yAxisMin) yAxisMin = yMin; end
  if isempty(yAxisMax) yAxisMax = yMax; end
  if yAxisMinMaxSetByTicks && (length(yTick)>1), yAxisMin = min(yTick); yAxisMax = max(yTick); end
  if (yAxisMin == 0) && isequal(yScale,'log'), yAxisMin = yAxisMax/1000;end

  % get the xLabelPosition
  yLabelX = logScaleSum(xScale,yAxisXpos,tickDirVal*yLabelOffset,abs(xMax-xMin));
  yLabelY = logScaleMean(yScale,[yAxisMin yAxisMax]);
end
  
% get location for axis
titleX = logScaleMean(xScale,[xAxisMin xAxisMax]);
if strcmp(yDir,'normal'), titleY = yMax;else titleY = yMin;end
titleAlignment = 'bottom';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set axis properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% and set limits
set(curAxis,'XLim',[xMin xMax]);
set(curAxis,'YLim',[yMin yMax]);

% turn matlab axis off
dpa.axisVisible = get(curAxis,'Visible');
set(curAxis,'Visible','off');

% turn title off
if ishandle(dpa.titleHandle) 
  set(dpa.titleHandle,'Visible','off');
end

% make our own title label
if ~isempty(titleStr)
  dpa.axis.title = text(titleX,titleY,titleStr,'HorizontalAlignment','center','VerticalAlignment',titleAlignment,'FontSize',labelFontSize,'FontAngle','oblique','FontName',fontName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw horizontal axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(whichAxis,{'both','horizontal'}))
  % horizontal line
  dpa.xAxis.main = line([xAxisMin xAxisMax],repmat(xAxisYpos,1,2),'Color','k','LineWidth',lineWidth,'LineStyle','-');

  % major ticks
  dpa.xAxis.majorTicks = line(repmat(xTick,2,1),repmat([xAxisYpos xAxisTickYpos]',1,length(xTick)),'Color','k','LineWidth',lineWidth,'LineStyle','-');

  % minor ticks
  if ~isempty(xMinorTick)
    dpa.xAxis.minorTicks = line(repmat(xMinorTick,2,1),repmat([xAxisYpos xAxisMinorTickYpos]',1,length(xMinorTick)),'Color','k','LineWidth',lineWidth,'LineStyle','-');
  end
  
  if ~xTickLabelHide
    % tick labels
    for iLabel = 1:length(xTickLabel)
      dpa.xAxis.tickLabel(iLabel) = text(xTick(iLabel),xAxisTickYpos,xTickLabel{iLabel},'VerticalAlignment',xTextAlignment,'HorizontalAlignment','center');
    end
  end

  % set the text label
  if ~isempty(xLabel)
    dpa.xAxis.label = text(xLabelX,xLabelY,xLabel,'VerticalAlignment',xTextAlignment,'HorizontalAlignment','center','FontAngle','oblique','FontSize',labelFontSize,'FontName',fontName);
  end

  % set the dpa variable in user data
  set(curAxis,'UserData',dpa);

  % set the y limits so that the axis displays
  if xAxisTickYpos>yMax
    set(curAxis,'YLim',[logScaleSum(yScale,yMin,-yAxisMargin,abs(yMax-yMin)) xAxisTickYpos]);
  end
  if xAxisTickYpos<yMin
    set(curAxis,'YLim',[xAxisTickYpos logScaleSum(yScale,yMax,yAxisMargin,abs(yMax-yMin))]);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw vertical axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(whichAxis,{'both','vertical'}))
  % horizontal line
  dpa.yAxis.main = line(repmat(yAxisXpos,1,2),[yAxisMin yAxisMax],'Color','k','LineWidth',lineWidth,'LineStyle','-');

  % major ticks
  dpa.yAxis.majorTicks = line(repmat([yAxisXpos yAxisTickXpos]',1,length(yTick)),repmat(yTick,2,1),'Color','k','LineWidth',lineWidth,'LineStyle','-');

  % minor ticks
  if ~isempty(yMinorTick)
    dpa.yAxis.minorTicks = line(repmat([yAxisXpos yAxisMinorTickXpos]',1,length(yMinorTick)),repmat(yMinorTick,2,1),'Color','k','LineWidth',lineWidth,'LineStyle','-');
  end
  
  % tick labels
  for iLabel = 1:length(yTickLabel)
    dpa.yAxis.tickLabel(iLabel) = text(yAxisTickXpos,yTick(iLabel),yTickLabel{iLabel},'VerticalAlignment','middle','HorizontalAlignment',yTextAlignment);
  end
  

  % set the dpa variable in user data
  set(curAxis,'UserData',dpa);

  % set the y limits so that the axis displays
  if yAxisTickXpos>xMax
    set(curAxis,'XLim',[logScaleSum(xScale,xMin,-xAxisMargin,abs(xMax-xMin)) yAxisTickXpos]);
  end
  if yAxisTickXpos<xMin
    set(curAxis,'XLim',[yAxisTickXpos logScaleSum(xScale,xMax,xAxisMargin,abs(xMax-xMin))]);
  end
end

% set the text label
if ~isempty(yLabel)
  dpa.yAxis.label = text(yLabelX,yLabelY,yLabel,'VerticalAlignment',yLabelTextAlignment,'HorizontalAlignment','center','FontAngle','oblique','FontSize',labelFontSize,'Rotation',90,'FontName',fontName);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    validateParam    %
%%%%%%%%%%%%%%%%%%%%%%%
function tf = validateParam(paramName,paramVals)

paramValsStr = '';
for iVal = 1:length(paramVals)
  paramValsStr = sprintf('%s %s',paramValsStr,paramVals{iVal});
end

% get parameter - this is a bit sneaky - we are getting the variable
% by its name from the caller environment
paramOriginalValue = evalin('caller',paramName);
param = lower(paramOriginalValue);


% set default tf
tf = true;

% check if the parameter was set to a single character value that we can parse
if length(param) == 1
  for iVal = 1:length(paramVals)
    if lower(paramVals{iVal}(1)) == param
      param = paramVals{iVal};
      break;
    end
  end
end

% now see if there is a match
if isempty(strcmp(param,paramVals))
  % if not display error
  disp(sprintf('(drawPublishAxis) %s (%s) should be one of:%s',paramName,paramOriginalValue,paramValsStr));
  tf = false;
  return
end

% set the variable in the caller to what we have changed it to
% again, this is a bit sneaky but makes for slightly cleaner code above
assignin('caller',paramName,param);

%%%%%%%%%%%%%%%%%%%%%
%    logScaleSum    %
%%%%%%%%%%%%%%%%%%%%%
function val = logScaleSum(scale,val1,mult2,val2);

% helper function which sums val1 and val2 (multiplied by mult2)
% in log or linear space depending on the value of scale

% check scale
if strcmp(scale,'linear')
  % linear sum
  val = val1+mult2*val2;
else
  % log sum - deal with sign of val2 properly (that is, if it is negative
  % then we still want mult2 to move left or right depending on the sign
  % so have to reverse sign
  if log10(val2) < 0
    val = log10(val1)-2*mult2*log10(val2);
  else
    val = log10(val1)+2*mult2*log10(val2);
  end
  val = 10^val;
end

%%%%%%%%%%%%%%%%%%%%%%
%    logScaleMean    %
%%%%%%%%%%%%%%%%%%%%%%
function val = logScaleMean(scale,vals)

n = length(vals);
if strcmp(scale,'linear')
  val = sum(vals)/n;
else
  val = 10^(sum(log10(vals))/n);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    logScaleSpace    %
%%%%%%%%%%%%%%%%%%%%%%%
function val = logScaleSpace(scale,ticks,spaceVal)

if strcmp(scale,'linear')
  val = min(ticks):spaceVal:max(ticks);
else
  % for log scale spaceVal is the number of ticks per interval
  val = [];
  for iTick = 1:length(ticks)-1
    thisInterval = (ticks(iTick+1)-ticks(iTick))/(spaceVal+1);
    val = [val ticks(iTick)+thisInterval*(1:spaceVal)];
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    setFigSize    %
%%%%%%%%%%%%%%%%%%%%%
function setFigSize(figSize)

% get figure handle
h = gcf;

% set units and get position
set(h,'Units','Centimeters');
pos = get(h,'Position');

% set colors
%set(gcf,'Color',[1 1 1]);
%set(gca,'Color',[1 1 1]);

% set paper position mode and units
set(h,'PaperPositionMode','Auto','PaperUnits','Centimeters');

% two-dimensional figSize
if length(figSize)==2
  pos(4) = figSize(2);
  pos(3) = figSize(1);
  figSize = figSize(1);
end

% set fig size according to JN style
switch figSize
 case 0
  % do nothing
  figSize = pos(3);
 case 1
  % single column
  figSize = 8.5;
 case 1.5
  % full column
  figSize = 11.6;
 case 2 
  % double column
  figSize = 17.6;
end

% set figure position
set(h,'Position',[pos(1), pos(2), figSize, pos(4)*figSize/pos(3)]);


