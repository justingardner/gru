% dispmatrix
%
%      usage: dispmatrix()
%         by: justin gardner
%       date: 06/04/03
%    purpose: displays a matrix by drawing the contents
%             by column or row 
%      usage: dispmatrix(matrix,'parameter',value,...);
% parameters: m         = matrix to display. if m is a
%                         string then display text rather
%                         than matrix
%             ndispcols = number of columns to display
%                         defaults to maximum  
%             ndisprows = number of rows to display            
%                         defaults to maximum
%             disptype  = 'bycol' or 'byrow' to display by
%                         column or rows, defaults to columns
%             rowsize   = size of rows
%             colsize   = size of cols
%             rowoffset = offset for rows, in pixels to get
%                         to columns set as e.g. '10*colsize'
%             coloffset = offset for column
%             rowlabel  = label for rows, defaults to nrows 
%             collabel  = label for cols, defaults to ncols 
%            tracecolor = color of traces, defaults to 'k' 
%       e.g.: 
%          m = rand(100,100);
%          dispmatrix(m,'disptype','byrow','rowoffset',3,'tracecolor','r');  
%
%
function dispmatrix(m,var1,val1,var2,val2,var3,val3,var4,val4,var5,val5,var6,val6,var7,val7,var8,val8,var9,val9,var10,val10)

% check input arguments
if ((nargin < 1) || (floor(nargin/2) == nargin/2))
  help dispmatrix;
  return
end

%  number of rows and columns
nrows = size(m,1);
ncols = size(m,2);

% number of columns or rows to dislay
ndispcols = ncols;
ndisprows = nrows;

% disp by rows or cols
disptype = 'bycol';

% size of rows and columns
rowsize = 10;
colsize = 10;

% number of rows and columns to offset display by
rowoffset = 0;
coloffset = 0;

% default labels
rowlabel = sprintf('%i',nrows);
collabel = sprintf('%i',ncols);

% color of traces
tracecolor = 'k';

% set user specified parameters
for i = 1:(nargin-1)/2
  if (isstr(eval(sprintf('val%i',i))))
    eval(sprintf('%s = ''%s'';',eval(sprintf('var%i',i)),eval(sprintf('val%i',i))));
  else
    eval(sprintf('%s = %i;',eval(sprintf('var%i',i)),eval(sprintf('val%i',i))));
  end
end

% display string instead of matrix
if (isstr(m))
  thandle = text(coloffset,-rowoffset+rowsize/2,m);
  set(thandle,'HorizontalAlignment','center');
  return
end
% position and size of matrix
mpos = [coloffset-colsize/2 -rowoffset (ndispcols+1)*colsize -(ndisprows+1)*rowsize];

% left side of matrix
plot([mpos(1) mpos(1)],[mpos(2) mpos(2)+mpos(4)],'k','LineWidth',2);hold on
plot([mpos(1) mpos(1)+colsize/3],[mpos(2) mpos(2)],'k','LineWidth',2);hold on
plot([mpos(1) mpos(1)+colsize/3],[mpos(2) mpos(2)]+mpos(4),'k','LineWidth',2);hold on

% right side of matrix
plot([mpos(1) mpos(1)]+mpos(3),[mpos(2) mpos(2)+mpos(4)],'k','LineWidth',2);
plot([mpos(1) mpos(1)-colsize/3]+mpos(3),[mpos(2) mpos(2)],'k','LineWidth',2);hold on
plot([mpos(1) mpos(1)-colsize/3]+mpos(3),[mpos(2) mpos(2)]+mpos(4),'k','LineWidth',2);hold on

% display row count
thandle = text(mpos(1),mpos(2)-rowsize,rowlabel);
set(thandle,'HorizontalAlignment','right');
set(thandle,'VerticalAlignment','bottom');
set(thandle,'Rotation',90);

% display col count
thandle = text(mpos(1),mpos(2),collabel);
set(thandle,'HorizontalAlignment','left');
set(thandle,'VerticalAlignment','bottom');


% display by cols
if (findstr(disptype,'bycols'))
  % reduce the number of display rows by two if
  % possible to make room for ellipses, but only
  % if we are trying to display less than the ful
  % matrix size
  nmaxrows = ndisprows;
  if (ndisprows < nrows)
    if (ndisprows > 2)
      ndisprows = ndisprows - 2;
    end
  end
  nmaxcols = ndispcols;
  if (ndispcols < ncols)
    if (ndispcols > 2)
      ndispcols = ndispcols - 2;
    end
  end
  % now plot data, column by column
  for i = 1:nmaxcols
    if ((nmaxcols < ncols) && (i >= (nmaxcols-2)))
      % just plot ellipses
      plot(coloffset+colsize*i+colsize/2,-(rowoffset+rowsize*ndisprows/2),'k.')
    else
      % set y as center of each row
      y = -(rowoffset-rowsize/2+(1:ndisprows)*rowsize);
      % get data in matrix
      x = m(1:ndisprows,i);
      % normalize to column width
      x = x-min(x);
      if (max(x) ~= 0)
	x = x/max(x);
      end
      x = (colsize-2*colsize/5)*x+coloffset-colsize+colsize/5+colsize*i;
      % plot the data
      plot(x,y,tracecolor);
    end
  end
  % make ellipises if necessary for missing columns
  if (ndisprows ~= nrows)
    while (ndisprows <= nmaxrows)
      plot(coloffset-colsize+colsize*(ndispcols-1)/2,-(rowoffset+rowsize*ndisprows+rowsize/2),'k.')
      ndisprows = ndisprows + 1;
    end
  end
% display by rows
else
  % reduce the number of display columns by two if
  % possible to make room for ellipses, but only
  % if we are trying to display less than the ful
  % matrix size
  nmaxcols = ndispcols;
  if (ndispcols < ncols)
    if (ndispcols > 3)
      ndispcols = ndispcols - 3;
    end
  end
  nmaxrows = ndisprows;
  if (ndisprows < nrows)
    if (ndisprows > 2)
      ndisprows = ndisprows - 2;
    end
  end
  % now plot data, column by column
  for i = 1:nmaxrows
    if ((nmaxrows < nrows) && (i >= (nmaxrows-2)))
      % just plot ellipses
      plot(coloffset-colsize+colsize*ndispcols/2,-(rowoffset+rowsize*i+rowsize/2),'k.')
    else
      % set x as center of each col
      x = coloffset-colsize/2+(1:ndispcols)*colsize;
      % get data in matrix
      y = m(i,1:ndispcols);
      % normalize to row height
      y = y-min(y);
      if (max(y) ~= 0)
	y = -y/max(y);
      end
      y = -((rowsize-2*rowsize/5)*y+rowoffset-(1/2)*rowsize+rowsize/5+rowsize*i);
      % plot the data
      plot(x,y,tracecolor);
    end
  end
  % make ellipses if necessary for missing columns
  ndispcols = ndispcols+1;
  if (ndispcols ~= ncols)
    while (ndispcols <= nmaxcols)
      plot(coloffset-colsize+colsize*ndispcols+colsize/2,-(rowoffset+rowsize*ndisprows/2),'k.')
      ndispcols = ndispcols + 1;
    end
  end
end
axis off

