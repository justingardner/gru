% myhist.m
%
%      usage: myhist(data,bins,color,doNormalize)
%             draws histogram using data and supplied bins
%             sets color to color and allows for data normalization
%             data can be a cell array of data to histogram
%             in different colors
%             - bins can be passed as a 'number of bins' or
%               as the actual bins
%         by: justin gardner modified by franco pestilli
%       date: 10/30/00 - 10/30/2009
%
function H = myhist(data,bins,color,doNormalize)

if ~any(nargin == [1:4])
  help myhist;
  return
end

% make sure that data is a cell array,
% this is for stacking multiple data sets
% into one histogram
data = cellArray(data);

% default arguments
if (exist('bins')~=1)||isempty(bins)
  alldata = [];
  for i = 1:length(data)
    alldata = union(data{i},alldata);
  end
  bins = min(alldata):(max(alldata)-min(alldata))/20:max(alldata);
end
if (length(bins) == 1) && (round(bins)==bins)
  alldata = [];
  for i = 1:length(data)
    alldata = union(data{i},alldata);
  end
  bins = min(alldata):(max(alldata)-min(alldata))/bins:max(alldata);
end

if ieNotDefined('color'),color='gray';,end
if ieNotDefined('offset'),offset = 0;,end

% init some variables
histval = [];
histvalTotal = zeros(1,length(bins)-1);
binsize = bins(2)-bins(1);
bins = bins(1:(length(bins)-1));

for dnum = 1:length(data);
  % get length of data.
  total = length(data{dnum});

  if ieNotDefined('doNormalize')
   doNormalize = 0;
  else
   normalizationFactor = total;
  end

  % set up variables
  histval(dnum,:) = 0; k = 1;

  % set colors
  if isequal(color,'gray')
    H.colors{dnum} = 1-(dnum/length(data))*[1 1 1];
  else
    H.colors{dnum} = (dnum/length(data))*color2RGB(color);
  end

  % cycle over each bin
  for i = bins;
    % find the number of data points in each bin
    count = (data{dnum} >= i);
    count = and(count,data{dnum} < (i+binsize));
    
    % keep the count in histval
    % normalize the histogram by the number of 
    % of total occurrences, this transforms the 
    % plot in a probability plot
    if doNormalize
     histval(dnum,k) = sum(count)/normalizationFactor;
    else
     histval(dnum,k) = sum(count);
    end
    
    % plot the rectangle
    if histval(dnum,k) > 0
      hRect = rectangle('Position',[bins(k) histvalTotal(k) binsize histval(dnum,k)]);
      set(hRect,'FaceColor',H.colors{dnum});
      set(hRect,'EdgeColor','w');
    end

    % update counter
    k = k + 1;
  end

  % update the height
  histvalTotal = histval(dnum,:)+histvalTotal;
  
end

% set axis accordingly
%minbin = bins(first(find(sum(histval))));
%maxbin = bins(last(find(sum(histval))));
minbin = min(bins);
maxbin = max(bins);
xaxis(minbin-binsize/2,maxbin+3*binsize/2);
yaxis(min(sum(histval)),max(sum(histval)));
set(gca,'XTick',(minbin+binsize/2):binsize:(maxbin+binsize/2));

% pack up variables into an output structure
if nargout>0
  H.bins = bins;
  H.histval = histval;
end