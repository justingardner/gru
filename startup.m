dbstop('if','error');

pathNames = {'~/proj/mrTools','~/proj/mgl','~/proj/gru','~/proj/grustim','~/proj/matlab/plugins','~/proj/steeve','~/proj/vistasoft','~/proj/mba'};

addpath('~/proj/gru/mac_helpers');
for i = 1:length(pathNames)
  if isdir(pathNames{i})
    addpath(genpath_exclude(pathNames{i},{'.git','.svn'}));
  end
end

if usejava('jvm')
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging', false);
end
format long g
beep off

% some preferences for mrTools
mrSetPref('site','Stanford');
mrSetPref('magnet',{'GE Discovery MR750','other'} ); 
mrSetPref('coil',{'Nova Medical 32-channel','Nova Medical 16-channel', 'other'} ); 
mrSetPref('pulseSequence',{'epi','tsense','gre', 'se', 'sense','other'} );
mrSetPref('pluginPaths','~/proj/gru/plugins');

% select some plugins
selectedPlugins = mrGetPref('selectedPlugins');
addedPlugin = false;
if ~any(strcmp('pRF',selectedPlugins))
  selectedPlugins{end+1} = 'pRF';
  addedPlugin = true;
end
if ~any(strcmp('gru',selectedPlugins))
  selectedPlugins{end+1} = 'gru';
  addedPlugin = true;
end
if addedPlugin
  mrSetPref('selectedPlugins',selectedPlugins);
end

% run user specific matlab startup. This will look for a matlab
% script like startupkenji to run, if it finds it on the path
% will run
username = getusername;
userStartup = sprintf('startup%s',username);
if exist(userStartup) == 2
  disp(sprintf('(startup) Running user specific startup: %s',userStartup));
  eval(userStartup);
else
  disp(sprintf('(startup) No user specific startup scrip found. If you need to, create one called: %s',userStartup));
end  


