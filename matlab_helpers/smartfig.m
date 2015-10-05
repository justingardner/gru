% smartfig.m
%
%      usage: smartfig
%         by: justin gardner & i-han chou
%       date: 06/11/01
%    purpose: figure that remembers where it was last sized to.
%       e.g.: smartfig or smartfig('mtsac')
%             To keep reusing the same figure do:
%             smartfig('mtsac','reuse');
%
function f = smartfig(event,fignum)

% set path
global initpath;
initpath = '~/';

% check arguments
if (nargin == 0)
  event = 'init';
  figtype = 'normal';
elseif (nargin == 1)
  figtype = fixBadChars(event);
  event = 'init';
elseif (nargin ~= 2)
  help smartfig;
  return;
end

% global variable to hold common information
% like the figure handle
global gVar;
global gFignum;

% check for reuse. this means to bring up the figure
% if it has not been closed
if exist('fignum','var') && isstr(fignum) && strcmp(fignum,'reuse')
  event = fixBadChars(event);
  for i = 1:length(gVar)
    if isfield(gVar(i),'fig') && ~isempty(gVar(i).fig) && (gVar(i).fig ~= 0) && strcmp(gVar(i).figtype,event)
      f = gVar(i).fig;
      figure(f);
      return
    end
  end
  % if we didn't find it, then we have to open it  
  figtype = event;
  event = 'init';
end

% if fignum is set to a character string
% it means the user wants to set the figure
% up for landscape / portrait printing
figorient = 'none';
if ((nargin == 2) && isstr(fignum) && ~strcmp(fignum,'reuse'))
  figtype = event;
  event = 'init';
  if (findstr(fignum,'landscape'))
    figorient = 'landscape';
  else
    figorient = 'portrait';
  end
end


% event loop, farms off events to proper handlers
switch lower(event)
  case {'init'}
    initHandler(figtype);
    f = gVar(gFignum).fig;
    if (strcmp('landscape',figorient))
      set(f,'PaperPosition', [0.25 0.25 10.5 8]);
      set(f,'PaperOrientation','landscape');
    elseif (strcmp('portrait',figorient))
      set(f,'PaperPosition', [0.25 0.25 8 10.5]);
      set(f,'PaperOrientation','portrait');
    end
  case {'close'}
    closeHandler(fignum);
  otherwise
    disp(sprintf('Unknown event %s',event));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handler to intialize the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(figtype)

% common variables
global gVar;
global gFignum;
SCREENXMAX = 4800;
SCREENXMIN = -1600;
SCREENYMAX = 1600;
SCREENYMIN = 0;

% get next figure number
if (isempty(gFignum)) gFignum = 1;, else gFignum = gFignum + 1;, end

% open the figure
gVar(gFignum).fig = figure;

% set figtype
gVar(gFignum).figtype = figtype;

% set the close handler functions for the figure
set(gVar(gFignum).fig,'CloseRequestFcn',sprintf('smartfig(''close'',%i)',gFignum));

% see if the init file has an initial position
initpos = readInit(figtype,'initpos');

% if it does then reset the position of the figure
if (~isempty(initpos))
  % check position to make sure it is not off screen
  if (initpos(1) > SCREENXMAX)
    initpos(1) = initpos(1)-SCREENXMAX;
  end
  if ((initpos(1)+initpos(3)) < SCREENXMIN)
    initpos(1) = initpos(1)+SCREENXMAX;
  end
  if (initpos(2) > SCREENYMAX)
    initpos(2) = initpos(2)-SCREENYMAX;
  end
  if ((initpos(2)+initpos(4)) < SCREENYMIN)
    initpos(2) = initpos(2)+SCREENYMAX;
  end
  % see if we are running on unix, and displaying on yoyodyne
  % since there is some annoying drift of the size
%  if (isunix)
%    [s w] = system('setenv | grep DISPLAY');
%    if (findstr('yoyodyne',w))
%      barheight = 22;
%      initpos(2) = initpos(2)+barheight;
%      initpos(4) = initpos(4)-barheight;
%    end
%  end
  % set the positio
  set(gVar(gFignum).fig,'Position',initpos);
end
% turn off menus
%set(gVar(gFignum).fig,'MenuBar','none');

% set up buttons
gVar(gFignum).buttons.close = ...
uicontrol('Style','pushbutton',...
          'String','Close',...
          'Callback',...
	  sprintf('smartfig(''close'',%i)',gFignum),...
	  'Position',buttonpos(1,1));

set(gVar(gFignum).fig,'Name',figtype);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the position of button on the page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos] = buttonpos(row,col)

bwidth = 50;
bheight = 20;
marginsize = 5;
buttontop=10;

% left position
pos(1) = marginsize + (bwidth+marginsize)*(col-1);
pos(2) = buttontop + (bheight+marginsize)*(row-1);
pos(3) = bwidth;
pos(4) = bheight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handler to close the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeHandler(fignum)

% common variables
global gVar;

% if gVar does not have fig number, then it means
% that we must have cleared all variables, so there
% is no info remaining about the figure. just close it
% without saving then.
if length(gVar) < fignum
  closereq
  return
end

% make sure the fig number is valid
if gVar(fignum).fig == 0
  closereq
  return
end

% make sure we are switched to the figure
figure(gVar(fignum).fig);

huh = get(gVar(fignum).fig);
if ~isfield(huh,'Position')
  disp(sprintf('HUH: Figure %i does not have position',gVar(fignum).fig));
  return
end

% get position and save in init file
writeInit(gVar(fignum).figtype,'initpos',get(gVar(fignum).fig,'Position'));

% reset figure number
gVar(fignum).fig = 0;
gVar(fignum).closetime = datestr(now);
gVar(fignum).closepos = huh.Position;
% close the figure
closereq;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to write a value to the init file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeInit(figtype,name,value)

% initfile path and name
global initpath;
initname = 'smartfig';

% just call it .smartfig.mat
fullinitname = sprintf('.smartfig.mat');

%disp(sprintf('Saving to: %s',[initpath fullinitname]));

% see if we have a file
finit = fopen([initpath fullinitname],'r');
if (finit ~= -1)
  % close the file
  fclose(finit);
  % get init variable since it already exisits
  load([initpath fullinitname])
end

% set the name,value pair
eval(sprintf('%s.%s.%s = value;',initname,figtype,name));

% write it to the file
if (str2num(first(version)) < 7)
  eval(sprintf('save %s %s;',[initpath fullinitname],initname));
else
  eval(sprintf('save %s %s -V6;',[initpath fullinitname],initname));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to read a value to the init file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value] = readInit(figtype,name)

% initfile path and name
global initpath;
initname = 'smartfig';

% just call it .smartfig.mat
fullinitname = sprintf('.smartfig.mat');

%disp(sprintf('smartfig %s',fullinitname));
% see if we have a file
try
  finit = fopen([initpath fullinitname],'r');
catch
  disp(sprintf('ERROR %s',lasterr));
  disp(sprintf('UHOH: ERROR OPENING %s',fullinitname));
  value = [];
  return;
end

% check if the init var exists
if (finit ~= -1)
  % close the file
  fclose(finit);
  % get init variable since it already exisits
  load([initpath fullinitname])
  % convert the name to initvar
  eval(sprintf('initvar = %s;',initname));
else
  % no init variable so the field does not exist
  value = [];
  return;
end

% get the name,value pair
if (isfield(initvar,figtype) & ...
    (isfield(eval(sprintf('initvar.%s',figtype)),name)))
  eval(sprintf('value = initvar.%s.%s;',figtype,name));
else
  value = [];
end

