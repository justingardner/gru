% dofmricns.m
%
%        $Id:$ 
%      usage: dofmricns
%         by: justin gardner
%       date: 03/17/2015
%    purpose: do initial processing stream for Stanford data brought down
%             from NIMS. Based on dofmrigru code.
%
%             To use this, first make a directory where you want
%             all the data saved on your local machine. This should
%             be of the format sXXXXYYYYMMDD eg. for subject s0012
%             taken on 20150312 you would do
%
%             cd ~/data
%             mkdir s001220150312
%             cd s001220150312
%
%                'stimfileDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    stimfiles form a specific directory.
%                'pdfDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    pdf files form a specific directory.
%                'numMotionComp=1': Set to 0 if you don't want to run MLR motion comp. Set to > 1 if you want
%                    to set multiple motionComp parameters (e.g. for motionComping two sets of scans taken at
%                    different resolutions)
%                'anatFilename=[]': Set this to a filename if you want to specify a particular name for the
%                    anatomy file (e.g. 'anatFilename=myanat'), you can also set to a cell array of filenames
%
%             This will sort through the specified datadir and copy
%             all of these files in the correct directories on your local
%             computer.
%
function retval = dofmricni(varargin)

% Default arguments
pdfDir = [];
stimfileDir = [];
numMotionComp = [];
getArgs(varargin,{'pdfDir=[]','stimfileDir=[]','numMotionComp=1','anatFilename=[]'});

% check to make sure we have the computer setup correctly to run epibsi, postproc and sense
if checkCommands == 0,return,end

% make sure that MLR is not running
mrQuit;

% check correct directory
if strcmp(getLastDir(pwd),'Pre') || isfile(fullfile(fileparts(pwd),'mrSession.mat'))
  if askuser(sprintf('(dofmricni) Current path is %s. Did you want to start in %s',pwd,fileparts(pwd)));
    cd('..');
  end
end

% see if this is the first preprocessing or the second one
if ~isdir('Pre')
  disp(sprintf('(dofmricni) Running Initial dofmricni process'));
  dofmricni1(stimfileDir,pdfDir,numMotionComp,anatFilename);
else
  disp(sprintf('(dofmricni) Running Second dofmricni process'));
  dofmricni2(anatFilename);
end

%%%%%%%%%%%%%%%%%%%%
%%   dofmricni1   %%
%%%%%%%%%%%%%%%%%%%%
function dofmricni1(fidDir,carextDir,stimfileDir,pdfDir,numMotionComp,tsense,anatFilename)

global dataDir;

% get the current directory
expdir = getLastDir(pwd);

% location of all directories. Should be a directory
% in there that has the same name as the current directory
% i.e. s00120090706 that contains the fid files and the
% car/ext files
if isempty(fidDir),fidDir = fullfile(dataDir,expdir,'raw');end
if isempty(carextDir),carextDir = fullfile(dataDir,expdir,'aux');end
if isempty(stimfileDir), stimfileDir = fullfile(dataDir,expdir,'aux');end
if isempty(pdfDir),pdfDir = fullfile(dataDir,expdir,'aux');end

% check the fiddir
if ~isdir(fidDir)
  disp(sprintf('(dofmricni1) Could not find fid directory %s',fidDir));
  return
end

% now get info about fids
[fidList fidListArray] = getFileList(fidDir,'fid');
fidList = getFidInfo(fidList);
fidList = sortFidList(fidList);

% get list of epi scans
epiNums = getEpiScanNums(fidList);
if isempty(epiNums)
  disp(sprintf('(dofmricni1) Could not find any epi files in %s',fidDir));
  return
end

% check to see if any of the epi scans have acceleration greater than 1
senseProcessing = isSenseProcessing(fidList,epiNums);

% set the tsense array
tsense = setTsense(tsense,length(epiNums));

% check the tsense scans, for ones in which we are
% accelerating make sure the setting's are correct
% otherwise if the scan looks like it should be accelerated
% then give warning
[tf fidList tsense volTrigRatio] = checkTsense(fidList,epiNums,tsense);
if ~tf && ~askuser('(dofmricni) Continue?')
  return
end

% if this is not a tsense run, then get volTrigRatio independently
if isempty(volTrigRatio)
  for iEPI = 1:length(epiNums)
    volTrigRatio(iEPI) = getVolTrigRatio(fidList{epiNums(iEPI)},0);
  end
end

% update fidList so that it will display getVolTrigRatio
for iEPI = 1:length(epiNums)
  fidList{epiNums(iEPI)}.dispstr = sprintf('%s volTrigRatio=%s',fidList{epiNums(iEPI)}.dispstr,mlrnum2str(volTrigRatio(iEPI),'compact',true));
end

% get anatomy scan nums
anatNums = getAnatScanNums(fidList,anatFilename);
if isempty(anatNums)
  disp(sprintf('(dofmricni1) Could not find any non-raw anatomies in %s',fidDir));
end

% get sense noise/ref scans
[senseNoiseNums senseRefNums] = getSenseNums(fidList);
if isempty(senseRefNums) && senseProcessing
  disp(sprintf('(dofmricni1) Could not find any sense ref files in %s',fidDir));
  return
end
if isempty(senseNoiseNums) && senseProcessing
  disp(sprintf('(dofmricni1) Could not find any sense noise files in %s',fidDir));
  return
end

% find car/ext files
carList = getFileList(carextDir,'car','ext');
carList = getCarInfo(carList);
if isempty(carList)
  disp(sprintf('(dofmricni1) Could not find any car/ext files in %s',carextDir));
end

% find pdf files
pdfList = getFileList(pdfDir,'pdf');

% get stimfile list
stimfileList = getFileList(stimfileDir,'mat');
stimfileList = setStimfileListDispStr(stimfileList);

% display what we found
dispList(fidList,epiNums,sprintf('Epi scans: %s',fidDir));
dispList(fidList,anatNums,sprintf('Anatomy scans: %s',fidDir));
if senseProcessing
  dispList(fidList,senseRefNums,sprintf('Sense reference scan: %s',fidDir));
  dispList(fidList,senseNoiseNums,sprintf('Sense noise scan: %s',fidDir));
end
dispList(pdfList,nan,sprintf('PDF files: %s',pdfDir));
dispList(stimfileList,nan,sprintf('Stimfiles: %s',stimfileDir));
dispList(carList,nan,sprintf('Car/Ext files: %s',carextDir));

% go find the matching car files for each scan
if ~isempty(carList)
  [carMatchNum carList] = getCarMatch(carList,fidList,epiNums);
  if isempty(carMatchNum),return,end
else
  carMatchNum = [];
end

% setup directories etc.
doMoveFiles(1,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense);

% make mask dir
if senseProcessing,makeMaskDir(1,fidList,anatNums,senseRefNums);end

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
fidList = doMoveFiles(0,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense);

% make mask dir
if senseProcessing,makeMaskDir(0,fidList,anatNums,senseRefNums);end

% save tsense settings
if ~isempty(tsense)
  save tsense tsense
end

% get subject id
exptname = getLastDir(pwd);
if (length(exptname) > 4) && (exptname(1) == 's')
  subjectID = exptname(1:4);
else
  subjectID = 'sxxx';
end

% now run mrInit
disp(sprintf('(dofmricni1) Setup mrInit for your directory'));
mrInit([],[],sprintf('subject=%s',subjectID),'makeReadme=0');

% set some info in the auxParams about the scane
v = newView;
if ~isempty(v)
  for iScan = 1:length(epiNums)
    % set fid name
    v = viewSet(v,'auxParam','fidFilename',fidList{epiNums(iScan)}.filename,iScan);
    % set volTrigRatio (only if different from 1)
    if volTrigRatio(iScan) ~= 1
      v = viewSet(v,'auxParam','volTrigRatio',volTrigRatio(iScan),iScan);
    end
    % set tsense field
    if ~isempty(tsense) && (length(tsense) >= iScan) && ~isempty(tsense{iScan})
      v = viewSet(v,'auxParam','tSense',tsense{iScan},iScan);
    end
    % set framePeriod as recorded in stimfile
    stimfile = viewGet(v,'stimfile',iScan);
    if ~isempty(stimfile)
      if strcmp(stimfile{1}.filetype,'mgl')
	% get all the volume events
	volEvents = find(stimfile{1}.myscreen.events.tracenum == 1);
	if length(volEvents) > 1
	  % get the framePeriod
	  framePeriod = median(diff(stimfile{1}.myscreen.events.time(volEvents)));
	  % now see if there are more volumes than acquisition triggers
	  if volTrigRatio(iScan)>1
	    framePeriod = framePeriod/volTrigRatio(iScan);
	  end
	  % round to nearest 1/1000 of a second
	  framePeriod = round(framePeriod*1000)/1000;
	  disp(sprintf('(dofmricni) Frame period as recorded in stimfile is: %0.3f',framePeriod));
  	  % set the frame period
	  scanParams = viewGet(v,'scanParams',iScan);
	  scanParams.framePeriod = framePeriod;
	  v = viewSet(v,'scanParams',scanParams,iScan);
	end
      end
    end
  end
  saveSession(0);
  deleteView(v);
end

% set up motion comp parameters
for i = 1:numMotionComp
  v = newView;
  if ~isempty(v)
    [v params] = motionComp(v,[],'justGetParams=1');
    deleteView(v);
    eval(sprintf('save motionCompParams%i params',i));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    removeTempFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function removeTempFiles(fidList)

% find all temp files to delete
deleteList = {};
fileExt = {'hdr','img','sdt','spr','edt','epr'};
for i = 1:length(fidList)
  filestem = stripext(fidList{i}.fullfile);
  for j = 1:length(fileExt)
    % check for files with that extension
    filename = setext(filestem,fileExt{j});
    if isfile(filename)
      deleteList{end+1} = filename;
    end
  end
end

% check for some other files
otherFiles = {'dofmricni2.log','logbook'};
for i = 1:length(otherFiles)
  filename = fullfile('Pre',otherFiles{i});
  if isfile(filename)
    deleteList{end+1} = filename;
  end
end

if ~isempty(deleteList)
  disp(sprintf('=============================================='));
  disp(sprintf('Temporary files'));
  disp(sprintf('=============================================='));
  for i = 1:length(deleteList)
    disp(sprintf('%i: %s',i,deleteList{i}));
  end
  disp(sprintf('=============================================='));

  if askuser('Found temporary files. Ok to remove them and place them in the directory Pre/deleteme')
    % make the directory if necessary
    if ~isdir('Pre/deleteme')
      mkdir(fullfile('Pre/deleteme'));
    end
    % move the files
    for i = 1:length(deleteList)
      moveToFilename = fullfile('Pre','deleteme',getLastDir(deleteList{i}));
      disp(sprintf('Moving %s -> %s',deleteList{i},moveToFilename));
      movefile(deleteList{i},moveToFilename);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%%   dofmricni2   %%
%%%%%%%%%%%%%%%%%%%%
function dofmricni2(movepro,tsense,getFiles,anatFilename,tsenseUseMask,tsenseUseNoise,dcCorrect)

% dc correct defaults to falst
if isempty(dcCorrect),dcCorrect = 0;end

% get fid file list
fiddir = 'Pre';
fidList = getFileList(fiddir,'fid');
fidList = getFidInfo(fidList);
fidList = sortFidList(fidList);

% get list of epi scans
epiNums = getEpiScanNums(fidList);

% remove any already existing temp files
removeTempFiles({fidList{epiNums}});

% get anatomy scan nums
anatNums = getAnatScanNums(fidList,anatFilename);

% check for something to do
if isempty(epiNumsWithCarExt)
  if ~askuser('(dofmricni) No epi scans with car/ext files. Process anyway?'),return,end
  epiNumsWithCarExt = epiNums;
% check to see if we should quit given that some of the peak files are missing
elseif length(epiNumsWithPeaks) ~= length(epiNums)
  if ~askuser('You are missing some peak files. Do you still want to continue'),return,end
end

% disp the scans
dispList(fidList,epiNumsWithCarExt,sprintf('Epi scans to process'));
dispList(fidList,anatNums,sprintf('Anatomy files'));

% display the commands for processing
processFiles(1,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect);

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
processFiles(0,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect);

%%%%%%%%%%%%%%%%%%%%%%
%%   processFiles   %%
%%%%%%%%%%%%%%%%%%%%%%
function processFiles(justDisplay,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect)

global postproc;
global senseCommand;
global tsenseCommand;

command = sprintf('cd Pre');
if justDisplay,disp(command),else,eval(command),end

% open the logfile
if ~justDisplay
  openLogfile('dofmricni2.log');
end

% copy anatomy files
for i = 1:length(anatNums)
  disp(sprintf('=============================================='));
  disp(sprintf('Copy processed anatomy into Anatomy directory'));
  disp(sprintf('=============================================='));
  if isfile(setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'hdr'))
    command = sprintf('copyfile %s ../Anatomy',setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'hdr'));
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('copyfile %s ../Anatomy',setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'img'));
    if justDisplay,disp(command),else,eval(command),end
  end
  if isfile(setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'nii'))
    command = sprintf('copyfile %s ../Anatomy',setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'nii'));
    if justDisplay,disp(command),else,eval(command),end
  end
end

% and delete headers
disp(sprintf('=============================================='));
disp(sprintf('Remove temporary nifti files'));
disp(sprintf('=============================================='));
command = 'mysystem(''rm -f ../Raw/TSeries/*.hdr'');';
if justDisplay,disp(command),else,eval(command);end
command = 'mysystem(''rm -f ../Raw/TSeries/*.img'');';
if justDisplay,disp(command),else,eval(command);end

disp(sprintf('=============================================='));
disp(sprintf('Run motion comp'));
disp(sprintf('=============================================='));
if ~justDisplay
  motionCompFilenameNum = 1;
  motionCompFilename = sprintf('motionCompParams1.mat');
  while isfile(motionCompFilename)
    % load the params
    eval(sprintf('load %s',motionCompFilename));
		       % run the motion comp
		       v = newView;
		       v = motionComp(v,params);
		       deleteView(v);
		       % get the next motionCompParams to run
		       motionCompFilenameNum = motionCompFilenameNum+1;
		       motionCompFilename = sprintf('motionCompParams%i.mat',motionCompFilenameNum);
  end
end

disp(sprintf('=============================================='));
disp(sprintf('DONE'));
disp(sprintf('=============================================='));

if ~justDisplay
  closeLogfile;
end

%%%%%%%%%%%%
%% myeval %%
%%%%%%%%%%%%
function myeval(command,justDisplay)

if justDisplay
  disp(command);
else
  eval(command);
end  
%%%%%%%%%%%%%%%%%%%%%
%%   doMoveFiles   %%
%%%%%%%%%%%%%%%%%%%%%
function fidList = doMoveFiles(justDisplay,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense)

disp(sprintf('=============================================='));
disp(sprintf('Making directories'));
disp(sprintf('=============================================='));

% list of directories to make
dirList = {'Etc','Pre','Doc','Pre/Aux','Raw','Raw/TSeries','Anatomy','Anal'};

% make them
for i = 1:length(dirList)
  if ~isdir(dirList{i})
    command = sprintf('mkdir(''%s'');',dirList{i});
    if justDisplay,disp(command),else,eval(command);,end
  end
end

% open the logfile
if ~justDisplay
  cd('Pre');
  openLogfile('dofmricni.log');
  cd('..');
  % write info about various scans
  dispList(fidList,epiNums,'Epi scans',true);
  dispList(fidList,anatNums,'Anatomy scans',true);
  if ~isempty(senseRefNums) dispList(fidList,senseRefNums,'Sense reference scan',true);end
  if ~isempty(senseNoiseNums) dispList(fidList,senseNoiseNums,'Sense noise scan',true);end
  dispList(carList,nan,'Car/Ext files',true);
  dispList(pdfList,nan,'PDF files',true);
  dispList(stimfileList,nan,'Stimfiles',true);
end

disp(sprintf('=============================================='));
disp(sprintf('Copying pdf files'));
disp(sprintf('=============================================='));

% move all the pdf files into the directory
for i = 1:length(pdfList)
  command = sprintf('copyfile %s %s',pdfList{i}.fullfile,fullfile('Doc',pdfList{i}.filename));
  if justDisplay,disp(command),else,eval(command);,end
end

disp(sprintf('=============================================='));
disp(sprintf('Copying stimfiles'));
disp(sprintf('=============================================='));

% move stimfiles
for i = 1:length(stimfileList)
  command = sprintf('copyfile %s %s',stimfileList{i}.fullfile,fullfile('Etc',stimfileList{i}.filename));
  if justDisplay,disp(command),else,eval(command);disp(command);,end
end

if justDisplay,disp(sprintf('+++++++'));end

disp(sprintf('=============================================='));
disp(sprintf('Copying fid files'));
disp(sprintf('=============================================='));

% move fid directories
allUsefulFids = [epiNums anatNums senseNoiseNums senseRefNums];
for i = 1:length(fidList)
  % fix name in our fidList as any . will be replaced by _ in our Pre
  % directory but not in /usr1 or wherever the files were copied from,
  fidList{i}.filename = setext(fixBadChars(stripext(fidList{i}.filename),{'.','_'}),'fid');
  % not a useful scan, put it in Pre/Aux
  if ~any(i==allUsefulFids)
    command = sprintf('copyfile %s %s',fidList{i}.fullfile,fullfile('Pre/Aux',fidList{i}.filename));
    if justDisplay,disp(command);else,eval(command);,end
    %otherwise put it in Pre
  else
    command = sprintf('copyfile %s %s',fidList{i}.fullfile,fullfile('Pre',fidList{i}.filename));
    if justDisplay,disp(command),else,eval(command);,end
  end
end

disp(sprintf('=============================================='));
disp(sprintf('DONE moving files.'));
disp(sprintf('=============================================='));

% convert anatomy scans
for i = 1:length(anatNums)
  disp(sprintf('=========================='));
  disp(sprintf('Convert %s to nifti  ',fidList{anatNums(i)}.filename));
  disp(sprintf('=========================='));
  % convert anatomy to nifti
  command = sprintf('fid2nifti %s %s;',fullfile('Pre',fidList{anatNums(i)}.filename),setext(fidList{anatNums(i)}.filename,'hdr'));
  if justDisplay,disp(command),else,eval(command);end
end


if ~justDisplay
  closeLogfile
end

%%%%%%%%%%%%%%%%%%%%%%
%%   dispScanlist   %%
%%%%%%%%%%%%%%%%%%%%%%
function dispStr = dispList(fidList,nums,name,toLog)

if nargin < 4,toLog = false;end

dispStr = {};
dispConOrLog(sprintf('============================='),toLog);
dispConOrLog(sprintf('%s',name),toLog);
dispConOrLog(sprintf('============================='),toLog);

  
% if nums is nan, show all files
if isnan(nums),nums = 1:length(fidList);end

for i = 1:length(nums)
  % get display string
  if isfield(fidList{nums(i)},'dispstr')
    dispstr = fidList{nums(i)}.dispstr;
  else
    dispstr = fidList{nums(i)}.filename;
  end
  % display the string
  dispConOrLog(dispstr,toLog);
end

% empty nums means to display all
if isempty(nums)
  disp('NO MATCHING FILES');
end

%%%%%%%%%%%%%%%%%%%%
%%   sortFidList  %%
%%%%%%%%%%%%%%%%%%%%
function fidList = sortFidList(fidList)

% sort by time stamp from log file
for i = 1:length(fidList)
  for j = 1:length(fidList)-1
    if fidList{j}.startTime > fidList{j+1}.startTime
      temp = fidList{j};
      fidList{j} = fidList{j+1};
      fidList{j+1} = temp;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getDatenumFromLogLine   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thisDatenum thisDatestr]= getDatenumFromLogLine(thisline)

thisDatenum = 0;thisDatestr = '';
if ~isempty(thisline)
  splitThisLine = strfind(thisline,': ');
  if ~isempty(splitThisLine)
    thisline = thisline(1:splitThisLine-1);
    [dow month day hhmmss year] = strread(thisline,'%s %s %s %s %s');
    thisDatestr = sprintf('%s-%s-%s %s',day{1},month{1},year{1},hhmmss{1});
    thisDatenum = datenum(thisDatestr);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   getEpiScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function epiNums = getEpiScanNums(fidList,nVolsCutoff)

if ieNotDefined('nVolsCutoff'),nVolsCutoff = 10;end

epiNums = [];
for i = 1:length(fidList)
  if ~isempty(fidList{i}.info)
    if fidList{i}.info.dim(4) > nVolsCutoff
      epiNums(end+1) = i;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getAnatScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function anatNums = getAnatScanNums(fidList,anatFilename)

% set the extension of the anatFilenames to search for to fid
anatFilename = cellArray(anatFilename);
for i = 1:length(anatFilename)
  anatFilename{i} = lower(setext(anatFilename{i},'fid'));
end

anatNums = [];
% look for 3D scans or one that matches anatFilename
for i = 1:length(fidList)
  % check for matching name
  if ~isempty(anatFilename) && any(strcmp(lower(fidList{i}.filename),anatFilename))
    anatNums(end+1) = i;
  elseif ~isempty(fidList{i}.info)
    % check for 3d scan
    if fidList{i}.info.acq3d
      % make sure it is not a raw scan
      if fidList{i}.info.dim(3) == length(fidList{i}.info.pss)
        anatNums(end+1) = i;
      end
    else
      % check for an anatomy sounding name (i.e. has the word "anat" in it
      if ~isempty(strfind(lower(fidList{i}.filename),'anat'))
	%          keyboard
	%          %make sure it is not a raw scan (has been epibsi5 processed)
	%          command = sprintf('thisanatprocpar = readprocpar(''%s'')',fidList{i}.fullfile);
	%          eval(command);
	%          if thisanatprocpar.ni ~= 0
	anatNums(end+1) = i;
	%         end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%
%%   getFileList   %%
%%%%%%%%%%%%%%%%%%%%%
function [fileList fileListArray] = getFileList(dirname,extList,matchExtList,filenameMatch)

fileList = {};fileListArray = {};

% make into a cell array
extList = cellArray(extList);

if ~ieNotDefined('matchExtList')
  matchExtList = cellArray(matchExtList);
else
  matchExtList = [];
end

% default to no filename matching 
if ieNotDefined('filenameMatch')
  filenameMatch = '';
else
  filenameMatch = cellArray(filenameMatch);
end

% open the directory
dirList = dir(dirname);
if isempty(dirList),return,end

% now go through the directory looking for matches
for i = 1:length(dirList)
  match = 0;
  % skip all . files
  if dirList(i).name(1) == '.',continue,end
  % skip all files that don't match the filename match if specified
  if ~isempty(filenameMatch)
    noFilenameMatch = 0;
    for j = 1:length(filenameMatch)
      if isempty(strfind(lower(dirList(i).name),lower(filenameMatch{j}))),noFilenameMatch=1;,end
    end
    if noFilenameMatch,continue,end
  end
  % get the file extension
  thisExt = getext(dirList(i).name);
  % check for match
  for j = 1:length(extList)
    if strcmp(extList{j},thisExt)
      % found the match
      if ~match
	% keep the name of the file
	fileList{end+1}.filename = dirList(i).name;
	fileListArray{end+1} = dirList(i).name;
	fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	fileList{end}.date = dirList(i).date;
	fileList{end}.datenum = dirList(i).datenum;
	% if we need to find matches
	if ~isempty(matchExtList)
	  fileList{end}.filename = dirList(i).name;
	  fileListArray{end+1} = dirList(i).name;
	  fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	  fileList{end}.path = dirname;
	  fileList{end}.(sprintf('%sfilename',extList{j})) = dirList(i).name;
	  % check for matching file
	  stemName = stripext(fullfile(dirname,dirList(i).name));
	  for k = 1:length(matchExtList)
	    if ~isfile(setext(stemName,matchExtList{j}))
	      % new system does not require matching ext files (everything stored in car)
	      %disp(sprintf('(dofmricni1:getFileList) No matching %s file for %s',matchExtList{j},dirList(i).name));
	    else
	      fileList{end}.(sprintf('%sfilename',matchExtList{j})) = setext(dirList(i).name,matchExtList{j});
	    end
	  end
	  
	end
      end
      match = 1;
    end
  end
end

%%%%%%%%%%%%%%%%%%
%%   mysystem   %%
%%%%%%%%%%%%%%%%%%
function status = mysystem(commandName)

% display to buffer
disp('=================================================================');
disp(sprintf('%s',datestr(now)));
disp(commandName);
disp('=================================================================');
% run the command
[status result] = system(commandName);
disp(result);


% write into the logfile
writeLogFile('\n=================================================================\n');
writeLogFile(sprintf('%s\n',datestr(now)));
writeLogFile(sprintf('%s\n',commandName));
writeLogFile('=================================================================\n');
writeLogFile(result);

%%%%%%%%%%%%%%%%%%%%%
%%   openLogfile   %%
%%%%%%%%%%%%%%%%%%%%%
function openLogfile(filename)

global gLogfile;

% open logfile
gLogfile.fid = fopen(filename,'w');
if gLogfile.fid == -1
  disp(sprintf('(dofmricni) Could not open logfile %s',filename))
  return
end

% remember filename
gLogfile.filename = filename;

%%%%%%%%%%%%%%%%%%%%%%
%%   closeLogfile   %%
%%%%%%%%%%%%%%%%%%%%%%
function closeLogfile

global gLogfile;

% close logfile
if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fclose(gLogfile.fid);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   writeLogFile   %%
%%%%%%%%%%%%%%%%%%%%%%
function writeLogFile(text)

global gLogfile;

if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fprintf(gLogfile.fid,text);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispConeOrLog    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispConOrLog(textStr,toLog)

if nargin < 2,toLog = false;end
if toLog
  writeLogFile(sprintf('%s\n',textStr));
else
  disp(textStr);
end
%%%%%%%%%%%%%%%%%%%%%%%
%    checkCommands    %
%%%%%%%%%%%%%%%%%%%%%%%
function retval = checkCommands

% commands to check
commandNames = {};
helpFlag = {};
for i = 1:length(commandNames)
  % suse which to tell if we have the command
  [commandStatus commandRetval] = system(sprintf('which %s',eval(commandNames{i})));
  % check for commandStatus error
  if commandStatus~=0
    disp(sprintf('(dofmricni) Could not find %s command: %s',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.stanford.edu/doku.php/gruprivate/stanford#computer_setup for help setting up your computer'));
    retval = 0;
    return
  end
  % run the command to see what happens
  [commandStatus commandRetval] = system(sprintf('%s %s',eval(commandNames{i}),helpFlag{i}));
  % check for commandStatus error
  if commandStatus>1
    disp(commandRetval);
    disp(sprintf('(dofmricni) Found %s command: %s, but could not run (possibly missing fink library?)',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmricni for help setting up your computer'));
    retval = 0;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setStimfileListDispStr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimfileList = setStimfileListDispStr(stimfileList)

for i = 1:length(stimfileList);
  % try to load it
  if isfile(stimfileList{i}.fullfile)
    stimfile = load(stimfileList{i}.fullfile);
    if ~isfield(stimfile,'myscreen')
      stimfileList{i}.dispstr = sprintf('%s (!!!No myscreen variable!!!)',stimfileList{i}.filename);
      continue;
    end
    myscreen = stimfile.myscreen;
    % make a string of some info myscreen
    stimfileStr = stimfileList{i}.filename;
    if isfield(myscreen,'volnum')
      stimfileStr = sprintf('%s [%i vols] ',stimfileStr,myscreen.volnum);
    end
    if isfield(myscreen,'starttime')
      stimfileStr = sprintf('%s%s ',stimfileStr,myscreen.starttime);
    end
    if isfield(myscreen,'endtime')
      stimfileStr = sprintf('%s(End: %s) ',stimfileStr,myscreen.endtime);
    end
    stimfileList{i}.dispstr = stimfileStr;
  end
end

