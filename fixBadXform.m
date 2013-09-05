% fixBadXform.m
%
%        $Id:$ 
%      usage: fixBadXform()
%         by: justin gardner
%       date: 09/04/13
%    purpose: This should be run from an mrTools session directory and will check to see if there is
%             a bad xform problem for the data set. i.e. After the console upgrade it seems that the
%             xforms computed for epis were incorrect (off by one in readout and phase enocde). This
%             was usually compensated by hand by moving the epis in mrAlign. This program will compute
%             how far from optimal that fix was. It can also replace all the xforms with the correct
%             ones. After doing that though, it is recommended to go check in mrAlign that everything
%             still looks good.
%
function retval = fixBadXform()

% check arguments
if ~any(nargin == [0])
  help fixBadXform
  return
end

mrQuit(0);

% get session info
[tf s] = getSessionInfo;
if ~tf,return,end

% display session info
s = dispSessionInfo(s);

% fix it
if any(s.needsFix)
  % display what the fix will do
  [tf s] = dispFixedVsOld(s)
  if ~tf,return,end
  % ask user if they want to fix
  if askuser('Fix')
    fixSession(s);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    dispFixedVsOld    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [tf s] = dispFixedVsOld(s)

tf = false;
% read the epi
filename = fullfile(s.sessionDir,'Raw','TSeries',s.tseriesFilename{1});
disppercent(-inf,sprintf('Reading %s',filename));
[data hdr] = mlrImageLoad(filename);

% remove all but one volume
data = data(:,:,:,1);

% write it out
epiOriginal = 'epiOriginal.hdr';
mlrImageSave(epiOriginal,data,hdr);

% fix the header
epiFixed = 'epiFixed.hdr';
hdr.sform = s.newSform{1};
mlrImageSave(epiFixed,data,hdr);

% check the old
anatFilename = fullfile(s.sessionDir,'Anatomy',s.anatFiles{s.anatNum});
if askuser(sprintf('CHECK original epi: %s alignment to session anatomy: %s',epiOriginal,anatFilename))
  mlrVol(anatFilename,epiFixed);
end

% check the fixed
if askuser(sprintf('CHECK fixed epi: %s alignment to session anatomy: %s',epiFixed,anatFilename))
  mlrVol(anatFilename,epiFixed);
end

tf = true;

%%%%%%%%%%%%%%%%%%%%
%    fixSession    %
%%%%%%%%%%%%%%%%%%%%
function fixSession(s)

% if not all the sforms are the same, then complain
needsFix = find(s.needsFix);
sform = s.newSform{needsFix(1)};
for i = 2:length(needsFix)
  if ~isequal(sform,s.newSform{needsFix(i)})
    disp(sprintf('(fixBadXform) Not all the fixed sforms are equal'));
    keyboard
  end
end

% save the sforms
saveSform(sform);

% if we got here then store what we did
s.date = now;
s.dateStr = datestr(s.date);
save fixBadXformDone s

%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispSessionInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function s = dispSessionInfo(s)

% display anatomy
dispHeader('Session Anatomy');
disp(sprintf('Filename: %s',s.anatFiles{s.anatNum}));
disp(sprintf('Fidname: %s',s.anatFidname));

% check that anat header matches
sigfigs = 1;
if ~all(all(nearlyequal(s.anatHeader.qform,s.anatXform,sigfigs)))
  disp(sprintf('(fixBadXform) Stored qform and qform recomputed from fid are mismatched by more than %i decimal place',sigfigs));
  disp(sprintf('Stored qform:\n%s',mlrnum2str(s.anatHeader.qform,'sigfigs',sigfigs,'compact',false)));
  disp(sprintf('Recomputed qform:\n%s',mlrnum2str(s.anatXform,'sigfigs',sigfigs,'compact',false)));
  % get what the sform should be
  disp(sprintf('(fixBadXform) Need to recompute sform based on *new* qform'));
  % not implemented yet
  keyboard
else
  disp(sprintf('Stored qform and qform recomputed from fid MATCH'));
  s.anatQform = s.anatHeader.qform;
  s.anatSform = s.anatHeader.sform;
end

% now display scans in raw
for iScan = 1:s.nScans
  dispHeader(sprintf('Scan Raw:%i',iScan));
  disp(sprintf('Filename: %s',s.tseriesFilename{iScan}));
  disp(sprintf('Fidname: %s',s.fidFilename{iScan}));
  % check for match of qform
  fidXform = s.fidXform{iScan};
  fidXformOriginal = s.fidXformOriginal{iScan};
  qform = s.qform{iScan};
  sform = s.sform{iScan};
  if ~nearlyequal(qform,fidXform,sigfigs)
    % check if it matches original xform
    if ~nearlyequal(qform,fidXformOriginal,sigfigs)
      disp(sprintf('Qform is mismatched by more than %i decimal place (as compared to what it should have been w/out the fix for the console bug)',sigfigs));
      disp(sprintf('Stored qform:\n%s',mlrnum2str(qform,'sigfigs',sigfigs,'compact',false)));
      disp(sprintf('Recomputed qform:\n%s',mlrnum2str(fidXformOriginal,'sigfigs',sigfigs,'compact',false)));
    else
      disp(sprintf('Qform is not up-to-date'));
    end
    s.needsFix(iScan) = true;
  else
    disp(sprintf('Qform OK'));
    s.needsFix(iScan) = false;
  end
  
  % compute recomputedSform
  newSform = s.anatSform*inv(s.anatQform)*fidXform;

  % check if they are the same as what is already there or not
  if nearlyequal(sform,newSform,sigfigs)
    disp(sprintf('Recomputed sform matches. No need to fix'));
    s.needsFix(iScan) = false;
  else
    % check rotation
    sigfigs = 0;
    if ~nearlyequal(sform(1:3,1:3),newSform(1:3,1:3),sigfigs)
      disp(sprintf('Recomputed sform rotation is off by more than %i decimal places',sigfigs))
      disp(sprintf('Stored sform:\n%s',mlrnum2str(sform,'sigfigs',sigfigs,'compact',false)));
      disp(sprintf('Recomputed sform:\n%s',mlrnum2str(newSform,'sigfigs',sigfigs,'compact',false)));
    end
    % check offsets
    offsetShift = round(sform(1:3,4)-newSform(1:3,4));
    if any(offsetShift)
      disp(sprintf('Offset of new sForm from old sForm: [%i %i %i]',offsetShift(1),offsetShift(2),offsetShift(3)));
    end
  end

  % store
  s.newSform{iScan} = newSform;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    getSessionInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [tf s] = getSessionInfo

tf = false;

% session directory is directory we started in
s.sessionDir = pwd;

% init the view
s.v = newView;
if isempty(s.v),return,end

% check for anat directory
s.anatDir = fullfile(s.sessionDir,'Anatomy');
if ~isdir(s.anatDir)
  disp(sprintf('(fixBadXform) Could not find dir: %s',s.anatDir));
  return
end

% get anat files
s.anatFiles = mlrImageGetAllFilenames(s.anatDir);

% check for at least one anat file
if length(s.anatFiles) < 1
  disp(sprintf('(fixBadXform) Could not find any session anatomy in dir: %s',s.anatDir));
  return
end
% if more than one anat file have user select which anat file
s.anatNum = 1;
if length(s.anatFiles) > 1
  for iAnat = 1:length(s.anatFiles)
    disp(sprintf('%i: %s',iAnat,s.anatFiles{iAnat}));
  end
  s.anatNum = getnum('(fixBadXform) Choose the session anatomy (0 to abort)',0:length(s.anatFiles));
  if s.anatNum == 0,return,end
end

% load the session anatomy header
s.anatHeader = mlrImageHeaderLoad(fullfile(s.anatDir,s.anatFiles{s.anatNum}));
if isempty(s.anatHeader),return,end

% go look for the fid that this anatomy was created from
s.preDir = fullfile(s.sessionDir,'Pre');
[tfGetFid s.anatFidname s.anatXform s.anatInfo] = getFidFile(s,s.anatFiles{s.anatNum});
if ~tfGetFid,return,end

% now switch to Raw group
rawGroupNum = viewGet(s.v,'groupNum','Raw');
if isempty(rawGroupNum),return,end
s.v = viewSet(s.v,'curGroup',rawGroupNum);

% see how many scans there are
s.nScans = viewGet(s.v,'nScans');
if (s.nScans < 1)
  disp(sprintf('(fixBadXform) No scans found in group Raw'));
  return
end

% load info for each scan (actually, only need to load for one scan - because we
% check that all qforms / sforms match later, so this has been changed to only load the first)
disppercent(-inf,'(fixBadXform) Loading scan info');
s.nScans = 1;
for iScan = 1:s.nScans
  % get tseries filename
  s.tseriesFilename{iScan} = viewGet(s.v,'tseriesfile',iScan);
  if isempty(s.tseriesFilename{iScan})
    disp(sprintf('(fixBadXform) Could not find tseries filename for scan Raw:%i',iScan));
    return
  end
  
  % get fid name
  s.fidFilename{iScan} = viewGet(s.v,'fidFilename',iScan);
  if isempty(s.fidFilename{iScan})
    disp(sprintf('(fixBadXform) Could not find tseries filename for scan Raw:%i',iScan));
    return
  elseif length(s.fidFilename{iScan}) > 1
    disp(sprintf('(fixBadXform) Scan Raw:%i has multiple fids associated with it',iScan));
    return
  else
    s.fidFilename{iScan} = s.fidFilename{iScan}{1};
  end

  % check that there is a valid fid
  if ~isdir(s.fidFilename{iScan})
    disp(sprintf('(fixBadXform) Could not find fid: %s',s.fidFilename{iScan}));
    return
  end
  
  % try to load the fidxform
  [s.fidXform{iScan} s.fidInfo{iScan}] = fid2xform(s.fidFilename{iScan});
  if isempty(s.fidXform{iScan}),return,end

  % try to load the fidxform w/out fix
  s.fidXformOriginal{iScan} = fid2xform(s.fidFilename{iScan},0,'fixConsoleUpgradeBug=0');
  if isempty(s.fidXform{iScan}),return,end
  
  % get current sform and qform
  s.qform{iScan} = viewGet(s.v,'scanqform',iScan);
  s.sform{iScan} = viewGet(s.v,'scansform',iScan);
  disppercent(iScan/s.nScans);
end
disppercent(inf);

qformRef = [];sformRef = [];
% now go through and read all the qforms and sforms
s.nGroups = viewGet(s.v,'nGroups');
for iGroup = 1:s.nGroups
  nScans = viewGet(s.v,'nScans',iGroup);
  for iScan = 1:nScans
    qform = viewGet(s.v,'scanQform',iScan,iGroup);
    sform = viewGet(s.v,'scanQform',iScan,iGroup);
    % on first qform/sform just store as reference
    if isempty(qformRef)
      qformRef = qform;
      sformRef = sform;
    else
      % on each other one, make sure there is a match
      if ~isequal(qform,qformRef)
	disp(sprintf('(fixBadXform) Qform for %s:%i does not match. All qforms must be the same to run this program',viewGet(s.v,'groupName',iGroup),iScan));
	return
      end
      if ~isequal(sform,sformRef)
	disp(sprintf('(fixBadXform) Sform for %s:%i does not match. All sforms must be the same to run this program',viewGet(s.v,'groupName',iGroup),iScan));
	return
      end
    end
  end
end
disp(sprintf('(fixBadXform) All scans match qform/sform'));
tf = true;


%%%%%%%%%%%%%%%%%%%%
%    getFidFile    %
%%%%%%%%%%%%%%%%%%%%
function [tf fidname fidXform fidInfo] = getFidFile(s,filename)

tf = false;fidname = [];fidXform = [];
fidname = setext(filename,'fid');
fidFullname = fullfile(s.preDir,fidname);

% check for file
if ~isdir(fidFullname)
  disp(sprintf('(fixBadXform) Could not find fid file: %s',fidFullname));
  return
end

% get xform
[fidXform fidInfo] = fid2xform(fidFullname);
if isempty(fidXform),return,end

% done
tf = true;
  


