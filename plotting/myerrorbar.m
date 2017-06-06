% myerrorbar.m
%
%      usage: myerrorbar(x,y,varargin)
%         by: justin gardner
%       date: 06/24/07
%    purpose: draw plots with error bars
%       e.g.: y error bars
%             myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10));
%
%             x error bars
%             myerrorbar(1:10,rand(1,10),'xError',0.5*rand(1,10));
%
%             x and yerror bars
%             myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10),'xError',0.5*rand(1,10));
%
%             different lower and upper bounds
%             myerrorbar(1:10,rand(1,10),'yLow',2*rand(1,10),'yHigh',0.5*rand(1,10));
%
%             draws a filled polygon rather than error bars
%             myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10),'yErrorBarType=fill','fillAlpha=0.2');
%   
%    options: Symbol = symbol to use, default 'o-'
%             Color = symbol color, default 'k'
%             MarkerFaceColor = symbol face color, defaults to Color
%             MarkerEdgeColor = symbol edge color, defaults to Color
%             MarkerSize = symbol size, default 8
%             myerrorbar(1:10,rand(1,10),'yError',rand(1,10)/2,'Symbol=s-','Color=[1 0.5 0]');
%             tee = draw tees or not, default 0
%             yTeelen = length of tee on y error, default to 1/10 of x spacing
%             xTeelen = length of tee on x error, default to 1/10 of y spacing
%             yErrorBarType = type of error bar can be: both, lower, upper, logy or fill 
%                             defaults to both which shows a vertical line between lower and upper
%                             error. Fill draws a transparent background with thickness of the error.
%             fillAlpha = 0.2 Only used for fill yErrorBarType. Sets the alpha of the fill region
%
% 
%
function retval = myerrorbar(x,y,varargin)

% check arguments
if nargin < 2
  help myerrorbar
  return
end
 
% check for old style usage
if (nargout == 1) || ((length(varargin) >= 1) && isnumeric(varargin{1}))
  retval = myerrorbarold(x,y,varargin);
  return
end

% get arguments
getArgs(varargin);

% no passed in x
if ieNotDefined('x');x = 1:length(y);end
  
% get length of x
n = length(x);

% get y upper and lower bounds
if ~ieNotDefined('yError'),yLow=yError;yHigh = yError;end
if ieNotDefined('yLow'),yLow = zeros(1,n);end
if ieNotDefined('yHigh'),yHigh = yLow;end
if ieNotDefined('yErrorBarType') yErrorBarType = 'both';end

% get x upper and lower bounds
if ~ieNotDefined('xError'),xLow=xError;xHigh = xError;end
if ieNotDefined('xLow'),xLow = zeros(1,n);end
if ieNotDefined('xHigh'),xHigh = xLow;end

% colors and symbols
if ieNotDefined('Symbol'),Symbol = 'o-';end
if ieNotDefined('Color')
  if ~ieNotDefined('MarkerFaceColor')
    Color = MarkerFaceColor;
  else
    Color = 'k';
  end
end
if ieNotDefined('MarkerEdgeColor'),MarkerEdgeColor='w';end
if ieNotDefined('MarkerFaceColor'),MarkerFaceColor=Color;end
if ieNotDefined('MarkerSize'),MarkerSize=8;end
if ieNotDefined('LineWidth'),LineWidth = 0.5;end
if ieNotDefined('fillAlpha'),fillAlpha = 0.2;end

% whether to draw tees or not
if ieNotDefined('tee'),tee = 0;end
if tee
  if ieNotDefined('yTeelen'),yTeelen = mean(diff(x))/10;end
  if ieNotDefined('xTeelen'),xTeelen = mean(diff(y))/10;end
end
hold on

% polar plot
if ~ieNotDefined('polarPlot') && polarPlot
  % x and y are angle and magnitude in degrees, so convert
  ang = pi*x(:)/180;mag = y(:);magMin = yLow(:);magMax = yHigh(:);
  % convert to radians
  x = cos(ang).*mag;
  y = sin(ang).*mag;
  % convert to x erorr 
  xLow = cos(ang).*(mag-magMin);
  xHigh = cos(ang).*(mag+magMax);
  yLow = sin(ang).*(mag-magMin);
  yHigh = sin(ang).*(mag+magMax);
  
  % plot the magnitude error bars
  for iError = 1:length(x)
    plot([xLow(iError) xHigh(iError)],[yLow(iError) yHigh(iError)],'-','Color',Color','LineWidth',LineWidth);
  end
  % wrap last point
  x(end+1) = x(1);
  y(end+1) = y(1);
% cartesian plot
else  
  % plot the y error bars
  if any(any(yLow ~= 0)) || any(any(yHigh ~= 0))
    % if we are doing a polygon fill of the error window then 
    if isequal(lower(yErrorBarType),'fill')
      x = x(:);y = y(:);yLow = yLow(:);yHigh = yHigh(:);
      fill([x;flipud(x)],[y-yLow;flipud(y+yHigh)],Color,'FaceAlpha',fillAlpha,'LineStyle','none');
      % draw regular error bars
    else
      for i = 1:length(x)
	switch yErrorBarType
	 case {'both','b'}
	  plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	 case {'lower','lo','l','bottom','bot'}
	  plot([x(i) x(i)],[y(i)-yLow(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
	 case {'higher','upper','up','top','hi'}
	  plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	 case {'logy'}
	  if (y(i)-yLow(i)) > 0
	    plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	  else
	    disp(sprintf('(myerrorbar) Dropping lower errorbar on %i which goes to %f',i,y(i)-yLow(i)));
	    plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	  end	 
	 otherwise
	  plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	end      
	% draw the tees if necessary
	if tee
	  plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)-yLow(i) y(i)-yLow(i)],'-','Color',Color,'LineWidth',LineWidth);
	  plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)+yHigh(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
	end
      end
    end
  end

  % plot the x error bars
  if any(any(xLow ~= 0)) || any(any(xHigh ~= 0))
    for i = 1:length(x)
      plot([x(i)-xLow(i) x(i)+xHigh(i)],[y(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
      % draw the tees if necessary
      if tee
	plot([x(i)-xLow(i) x(i)-xLow(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
	plot([x(i)+xHigh(i) x(i)+xHigh(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
      end
    end
  end
end

% plot the symbols
plot(x,y,Symbol,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'Color',Color,'MarkerSize',MarkerSize,'LineWidth',LineWidth);





