function h = errbar(x,y,varargin)
% errbar  plots vertical or horizontal errorbars. Error bars may be
% symmetic or asymmetric about their data points. 
% 
%% Syntax 
% 
%  errbar(x,y,E)
%  errbar(x,y,L,U)
%  errbar(...,'horiz')
%  errbar(...,LineSpec)
%  h = errbar(...)
% 
%% Description 
% 
% errbar(x,y,E) plots vertical lines of length 2*E centered on x,y. E
% may be a scalar or a vector matching the size of x. 
% 
% errbar(x,y,L,U) specifies different lower and upper bounds of error
% bars. L and U may be mixed vector and scalar. 
%
% errbar(...,'horiz') plots error bars horizontally. 
% 
% errbar(...,LineSpec) specifies line type. 
%
% h = errbar(...) returns a handle h of plotted error bars. 
% 
%% Usage note
% This function is very similar to Matlab's built-in errorbar function,
% but the two functions differ in the following ways: 
% 
% * errbar supports horizontal error bars with the simple inclusion of
% the 'horiz' tag. 
% 
% * errbar does not plot perpendicular lines at the ends of error bars,
% because in my experience, those lines tend to be garrish and distracting,
% and they do not automatically rescale when zooming. If you'd like
% perpendicular lines at the ends of your errbars, you can get a similar
% effect by specifying a LineStyle '+-'. 
% 
% * errbar does not plot markers at datapoints x,y. This is intentional,
% because it's rare that I want my error bars to have the same linewidth
% and bold coloring 
% 
% * errbar allows error values to be entered as a scalar, whereas Matlab's
% built-in errorbar function forces users to populate an array of redundant
% values. 
% 
%% Example 1: Simple vertical error bars
% % Stealing this example from the documentation for Matlab's built-in
% % errorbar function: 
% 
% X = 0:pi/10:pi;
% Y = sin(X);
% E = std(Y);
% 
% plot(X,Y,'bo-','linewidth',3)
% hold on
% 
% errbar(X,Y,E) % <- Here's the errbar function
% 
%% Example 2: Formatted error bars 
% % Formatting with errbar is as simple as any LineSpec formatting. Each
% % error bar is made of only two endpoints and a line connecting them, so
% % any marker properties will be applied to the ends of the error bars.
% % Here's a simple red error bar with plus signs at its ends: 
% 
% x = 1:5; 
% y = sin(x); 
% 
% figure
% plot(x,y,'bo','markersize',10)
% hold on
% 
% errbar(x,y,.2,'r+-')
% 
%% Example 3: Horizontal error bar
% % Now we'll add horizontal error bars to the plot above.  For the lower
% % bound we'll specify a scalar, and for the upper bound some random
% % numbers.  Make the error bars thick and black with 'k','linewidth',3 : 
% 
% xLowerUncertainty = .1; 
% xUpperUncertainty = (2*rand(size(x))).^2; 
% 
% errbar(x,y,xLowerUncertainty,xUpperUncertainty,'k','linewidth',3,'horiz')
% 
%% Author Info 
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG), April 2015. http://www.chadagreene.com.
% 
% Updated June 12, 2015 to turn off error bar handle visibility, thereby preventing
% error bars from showing up in a plot legend.  Thanks for this suggestion, FEX 
% user David! 
% 
% Updated July 17, 2015 to allow NaN inputs.  Thanks for this suggestion,
% FEX user Ian Craig! 
% 
% See also errorbar, plot, and annotation.

%% Initial error checks: 

assert(nargin>2,'errbar function requires at least three inputs.')
assert(isnumeric(x)==1,'errbar function input x must be numeric.')
assert(isnumeric(y)==1,'errbar function input y numeric.') 
assert(isnumeric(varargin{1})==1,'The third input to the errbar function must be numeric.') 
assert(numel(x)==numel(y),'errbar inputs x and y must be the same size.') 

%% Set defaults: 

% Set default errobar direction: 
PlotVertical = true; 

%% Parse inputs: 

% An array containing ones for any formatting-related input data:  
formatting = ones(length(varargin),1); 

% Set lower errorbar limit: 
L = varargin{1}; 
formatting(1) = 0; 

% Assume lower and upper limits are equal unless otherwise specified by user:  
U = L;

if nargin>3
    if isnumeric(varargin{2})
        U = varargin{2}; 
        formatting(2) = 0; 
    end
    
    tmp = strncmpi(varargin,'horiz',3); 
    if any(tmp) 
        PlotVertical = false; 
        formatting(tmp) = 0; 
    end
end

% Clip inputs, assuming everything leftover is linespec formatting:  
varargin = varargin(logical(formatting)); 

%% Check errorbar values 

assert(all(L(isfinite(L))>=0)==1,'Errorbar lower limits must be greater than or equal to zero.') 
assert(all(U(isfinite(U))>=0)==1,'Errorbar lower limits must be greater than or equal to zero.') 

%% Plot errorbars: 


if PlotVertical
    if isempty(varargin)
        h = plot([x(:) x(:)]',[y(:)-L(:) y(:)+U(:)]','b',varargin{:}); 
    else
        h = plot([x(:) x(:)]',[y(:)-L(:) y(:)+U(:)]',varargin{:}); 
    end
else
    if isempty(varargin)
        h = plot([x(:)-L(:) x(:)+U(:)]',[y(:) y(:)]','b',varargin{:}); 
    else
        h = plot([x(:)-L(:) x(:)+U(:)]',[y(:) y(:)]',varargin{:}); 
    end
end

% Turn off handle visibilty so error bars won't show up in legend: (Thanks for this suggestion, FEX user David.): 
set(h(:),'handlevisibility','off');


%% Clean up:

if nargout==0
    clear h
end

% To errbar is human. 
        
    
