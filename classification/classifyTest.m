% svmtest.m
%
%      usage: [x1 x2 svm] = classifyTest;
%         by: justin gardner
%       date: 08/01/05
%    purpose: Test program for classifiers. Can be used to generate two data sets of points,
%             display points and clasisifcation and run various forms of classification - see below.
%
%       e.g.: create two distributions and display results of svm classifier. Note that you
%             click on the graph to make points in the first distribution. Then you 
%             double click to make the last point. You next make a second distribution by
%             clicking a new set of points and double clicking. The return arguments are
%             the two distributions and a structure holding info about the classification
%
%             [x1 x2 c] = classifyTest;
%
%             You can also run on an already created distribution of points. Note that x1 and x2
%             points are n x 2 (containing n points with x and y as returned from classifyTest)
%
%             [x1 x2 c] = classifyTest(x1,x2);
%
function [x1,x2,classifier] = classifyTest(varargin)

% parse arguments
[tf x1 x2 classifyArgs] = parseArgs(varargin);
if ~tf,return,end

% set up figure
f = dispClassifyFigure(x1,x2,classifyArgs);

% choose some points for the second distribution if not passed in
if isempty(x1)
  title(sprintf('Click to add points to first distribution\ndouble-click to end'));
  [x1 y1] = getpts(f.fignum);x1 = [x1 y1];
end

% plot the first distribution
plotPoints(x1);

% choose some points for the second distribution if not passed in
if isempty(x2)
  title(sprintf('Click to add points to second distribution\ndouble-click to end'));
  [x2 y2] = getpts(f.fignum);x2 = [x2 y2];
end

% plot the points
plotPoints([],x2);

% if we are just defining points
if strcmp(classifyArgs.type,'none'),return,end

% build a classifier
classifier = buildClassifier({x1 x2},'type',classifyArgs.type,'kernelfun',classifyArgs.kernelFun,'kernelargs',classifyArgs.kernelArgs,'C',classifyArgs.C,'projectionLine',classifyArgs.projectionLine,'biasPoint',classifyArgs.biasPoint,'w',classifyArgs.w);

% calculate decision surface
decisionSurface = getDecisionSurface(classifier,f);

% dispay the decision surface
dispDecisionSurface(decisionSurface,f);

% plot the points again with the svm 
plotPoints(x1,x2,classifier,classifyArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display decision surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispDecisionSurface(decisionSurface,f)

% the decisionSurface is a number for each point (which is a projection
% on to the linear classifier). We now normalize the points so that
% we can display as going from red to blue according to how close they
% are to one or the other distribition. First normalize to a value between
% 1 and 256
maxcolor = max(max(max(decisionSurface)),min(min(decisionSurface)));
decisionSurface = floor(255*(decisionSurface+maxcolor)/(2*maxcolor))+1;

% display as an image
image(f.x,f.y,decisionSurface');

% set the colormap to "cool"
cool256 = cool(256);
colormap(cool256);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get decision surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decisionSurface = getDecisionSurface(classifier,f)

% calculate decision surface
disppercent(-inf,'Calculating decision surface');
for i = 1:length(f.x)
  disppercent(i/length(f.y));
  for j = 1:length(f.y)
    % calculate the classification for each point in the display region
    % this should be a value between -inf and inf indicating how far along
    % the projection each point is. The more negative the more likely the point
    % is to have come from the first distribution and the more positive the
    % more likely to have come from the other distribution. We set to the
    % raw classification output matrix: classMatrix
    [class classValue classMatrix] = classifyInstance(classifier,[f.x(i) f.y(j)]);
    decisionSurface(i,j) = -classMatrix(1,1);
  end
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPoints(x1,x2,classifier,classifyArgs)

if (exist('x1') && ~isempty(x1))
  % plot distributions
  plot(x1(:,1),x1(:,2),'ro','MarkerFaceColor','r');
end
hold on
if (exist('x2') && ~isempty(x2))
  plot(x2(:,1),x2(:,2),'bo','MarkerFaceColor','b');
  % expand axis
end

% if svm has been passed in then display the support vectors, title and decision boundary
if nargin > 2

  % display misclassified points
  misClassifiedPointsX1 = find(classifyInstance(classifier,x1) ~= 1);
  plot(x1(misClassifiedPointsX1,1),x1(misClassifiedPointsX1,2),'yx');
  misClassifiedPointsX2 = find(classifyInstance(classifier,x2) ~= 2);
  plot(x2(misClassifiedPointsX2,1),x2(misClassifiedPointsX2,2),'yx');
  nMisClassify = length(misClassifiedPointsX1) + length(misClassifiedPointsX2);
  
  % plot support vectors
  if strcmp(classifier.type,'svm')
    for i = 1:classifier.svm(2).n
      % only plot support vectors that are not pegged
      % at the C value
      if abs(classifier.svm(2).alpha(i)) ~= classifier.svm(2).C
	plot(classifier.svm(2).sv(i,1),classifier.svm(2).sv(i,2),'wo');
      end
    end
  end
  
  % display title
  if strcmp(classifier.type,'svm')
    svm = classifier.svm(2);
    title(sprintf('%s classifer\n(%i/%i support vectors)\narg = %s C = %s %i misclassified points',svm.kernelfun,svm.n-sum(abs(svm.alpha)==svm.C),length(x1)+length(x2),num2str(svm.kernelargs),num2str(svm.C),nMisClassify));
  elseif strcmp(classifier.type,'userDefinedLinear')
    title(sprintf('%s classifer: %i misclassified points',classifier.type,nMisClassify));
  else
    title(sprintf('%s %s classifer: %i misclassified points',classifier.type,classifier.svm(2).kernelfun,nMisClassify));
  end

  % display the decision boundary for linear machines
  % Decision boundary for linear kernel
  a = axis;
  % get weights and bias for the classification of 1 vs 2
  [w bias svm] = getClassifierWeights(classifier,1,2);

    if ~isempty(w)
    % get x points
    x = [a(1) a(2)];

    % draw decision boundary
    y = -(w(1)/w(2))*x-bias/w(2);
    plot(x,y,'w-');

    % draw projection line for userDefinedLinear
    if strcmp(classifier.type,'userDefinedLinear') && isfield(classifier,'projectionLine') && ~isempty(classifier.projectionLine)
      p1 = classifier.projectionLine(1,:);
      p2 = classifier.projectionLine(2,:);
      x = [-1.5 1.5];
      y = ((p2(2)-p1(2))/(p2(1)-p1(1)))*(x-p1(1))+p1(2);
      plot(x,y,'w--');
    end
    
    % draw margin boundaries for svm
    if strcmp(classifier.type,'svm') && isequal(svm.kernelfun,'linear')
      y = -(w(1)/w(2))*x-(1+bias)/w(2);
      plot(x,y,'y:');
      y = -(w(1)/w(2))*x-(-1+bias)/w(2);
      plot(x,y,'y:');
    end
  end
end



%%%%%%%%%%%%%%%%%
%    getArgs    %
%%%%%%%%%%%%%%%%%
function [tf x1 x2 classifyArgs] = parseArgs(args)

% default return values
tf = true;
x1 = [];x2 = [];
classifyArgs = [];

% first see if we are passed in data points for the first distribution
if (length(args) > 0) && isnumeric(args{1})
  x1 = args{1};
  args = {args{2:end}};
end

% see if we are passed in points for the second distribition
if (length(args) > 0) && isnumeric(args{1})
  x2 = args{1};
  args = {args{2:end}};
end

% check that the first distribution is good.
if ~isempty(x1)
  warnStr = '(classifyTest:parseArgs) First distribution of points must be a n x 2 vector';
  % check if it has two dimensions
  if (ndims(x1) ~= 2),disp(warnStr);,tf = false;return, end
  % and the second dimension is 2
  if (size(x1,2) ~= 2)
    x1 = x1';
    if size(x1,2) ~= 2,disp(sprintf(warnStr));tf = false;return,end
  end
end

% check that the second distribution is good.
if ~isempty(x2)
  warnStr = '(classifyTest:parseArgs) Second distribution of points must be a n x 2 vector';
  % check if it has two dimensions
  if (ndims(x2) ~= 2),disp(warnStr);,tf = false;return, end
  % and the second dimension is 2
  if (size(x2,2) ~= 2)
    x2 = x2';
    if size(x2,2) ~= 2,disp(sprintf(warnStr));tf = false;return,end
  end
end
    
% deal with other arguments
argDefaults = {'type=svm','kernelFun=linear','kernelArgs=0','C=1000','projectionLine=[]','biasPoint=[]','w=[]','tightAxis=0'};
getArgs(args,argDefaults);

% check kernel function
validKernelFuns = {'linear' 'polynomial' 'radialbasis' 'sigmoid' 'fisher' 'fisherpillow'};
if ~any(strcmp(lower(kernelFun),validKernelFuns))
  validKernelFunsStr = '';
  for i = 1:length(validKernelFuns)
    validKernelFunsStr = sprintf('%s %s',validKernelFunsStr,validKernelFuns{i});
  end
  disp(sprintf('(classifyTest) Kernel function should be one of %s',validKernelFunsStr));
  tf = false;
  return
end

% check for valid projectionLine
if ~isempty(projectionLine)
  if ~isnumeric(projectionLine) || ~isequal(size(projectionLine),[2 2])
    disp(sprintf('(classifyTest) ProjectionLine must be a 2x2 matrix with each row containing one point [x y] that defines the projection line'));
    tf = false;
    return
  end
  type = 'userDefinedLinear';
end

% check for valid 
if ~isempty(w)
  type = 'userDefinedLinear';
end

% wrap up arguments into a returnable structure
for iArg = 1:length(argDefaults)
  thisArg = strtok(argDefaults{iArg},'=');
  classifyArgs.(thisArg) = eval(thisArg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispClassifyFigure    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = dispClassifyFigure(x1,x2,classifyArgs)

% start the figure
f.fignum = mlrSmartfig('classifyTest','reuse');
clf;

% default axis limits
defaultXLim = [-1.5 1.5];
defaultYLim = [-1.5 1.5];
xLim = defaultXLim;
yLim = defaultYLim;

% set x and y axis limits according to the passed in input
allx = [x1 ; x2];
if ~isempty(allx)
  xLim(1) = min(allx(:,1));
  xLim(2) = max(allx(:,1));
  yLim(1) = min(allx(:,2));
  yLim(2) = max(allx(:,2));
end

if ~classifyArgs.tightAxis
  % make at least the size of the defauls
  xLim(1) = min(xLim(1),defaultXLim(1));
  xLim(2) = max(xLim(2),defaultXLim(2));
  yLim(1) = min(yLim(1),defaultYLim(1));
  yLim(2) = max(yLim(2),defaultYLim(2));
end

% set up axis  
xaxis(xLim(1),xLim(2));
yaxis(yLim(1),yLim(2));
axis off
hold on

% get some limits for displaying decision surfaces
f.x = xLim(1):(xLim(2)-xLim(1))/30:xLim(2);
f.y = yLim(1):(yLim(2)-yLim(1))/30:yLim(2);
