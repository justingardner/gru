% switchEditor.m
%
%        $Id:$ 
%      usage: switchEditor(whichEditor)
%         by: justin gardner
%       date: 12/11/09
%    purpose: Set whichEditor = 1 if you want emacs or whichEditor = 0 for emacs. If you
%             omit, should switch back and forth
%
function retval = switchEditor(whichEditor)

% check arguments
if ~any(nargin == [0 1])
  help switchEditor
  return
end

% chosse which emacs to us (prefer Aquamacs if it is there)
aquamacs = '/Applications/Aquamacs.app';
if isdir(aquamacs)
  whichEmacs = aquamacs;
else
  whichEmacs = '/Applications/Emacs.app';
end

if nargin == 0
  % check which editor we think we are using
  gCurrentBuiltInEditor = com.mathworks.services.Prefs.getBooleanPref('EditorBuiltinEditor');
  if gCurrentBuiltInEditor
    whichEditor = 1;
  else
    whichEditor = 0;
  end
end

if whichEditor == 1
  com.mathworks.services.Prefs.setStringPref('EditorOtherEditor',whichEmacs);
  com.mathworks.services.Prefs.setBooleanPref('EditorBuiltinEditor',false);
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false);
  disp(sprintf('(switchEditor) Switching to Emacs'));
elseif whichEditor == 2
  com.mathworks.services.Prefs.setStringPref('EditorOtherEditor','/Applications/TextMate.app');
  com.mathworks.services.Prefs.setBooleanPref('EditorBuiltinEditor',false);
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false);
  disp(sprintf('(switchEditor) Switching to TextMate'));
elseif whichEditor == 3
  % To use this option, create a shortcut to Sublime Text 3: 
  % ln -s /Applications/Sublime\ Text\ 3.app/Contents/SharedSupport/bin/subl /usr/local/bin/subl3
  com.mathworks.services.Prefs.setStringPref('EditorOtherEditor','/usr/local/bin/subl3');
  com.mathworks.services.Prefs.setBooleanPref('EditorBuiltinEditor',false);
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false);
  disp(sprintf('(switchEditor) Switching to Sublime Text 3'));
else
  com.mathworks.services.Prefs.setStringPref('EditorOtherEditor','');
  com.mathworks.services.Prefs.setBooleanPref('EditorBuiltinEditor',true);
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',true);
  disp(sprintf('(switchEditor) Switching to built-in editor'));
end

%edit(fullfile(prefdir,'matlab.prf'));
