% dofmrigru.m
%
%        $Id:$ 
%      usage: dofmrigru
%         by: justin gardner
%       date: 07/06/09
%    purpose: do initial processing stream from fid files collected on the Varian magnet
%             to MLR nifti file format. You will run this program twice to
%             process this data. The first time files will get moved into
%             your local computer. You then have to manually generate a
%             mask and run peak. Then you run it a second time to do the
%             sense reconstruction and run physiofix, etc.
%
%             To use this, first make a directory where you want
%             all the data saved on your local machine. Use the same directory name
%             as the one on /usr1 or /usr4 which contains the raw data. For example:
%
%             cd ~/data
%             mkdir S00120091216
%             cd S00120091216
%
%             Now, make sure that your raw data is setup as follows:
%
%             /usr1/justin/data/S00120091216 should contain all of your data. Including
%             all of the .fid directories)
%             all of the car/ext files
%             all of the stimfiles (if you have them)
%             any pdf files (documents if you have them)
% 
%             Now you can run dofmrigru as follows:
%
%             dofmrigru('dataDir=/usr1/justin/data')
%
%             Other options:
%                'epibsi=epibsi': set which epibsi processing function to use
%                'navCorrectMag=1': use navigator magnitude correction with epibsi
%                'navCorrectPhase=1': use navigator phase correction with epibsi
%                'dcCorrect=[]': use DC correction with epibsi and postproc. Set to 1 to use, defaults to not use
%                'refScan=1': do reference scan correction
%                'postproc=pp': set which postproc program to use
%                'tsenseCommand=tsense_test': set which tsense recon program to use
%                'tsense=1': set to 0 to turn off using tsense. If you want to
%                   set particular acceleration factor set the acc factor
%                   and shot order (e.g. [3 3 2 1]). Or set a cell array
%                   with the acc factor/shot order array for each epi in your scan. If
%                   set to 1, will use the ideal tsense acceleration factor 
%                   according to the number of shots and interleaving
%                'notchFilter=1': set to 0 if you don't want to notch filter the tsense data
%                'senseCommand=sense_mac_intel': set which sense reconstruction to use
%                'fidDir=/usr1/justin/data/s00620101001/Pre': Set this if you want to load the first pass
%                    fid files form a specific directory.
%                'carextDir=/usr1/justin/data/s00620101001/carext': Set this if you want to load the first pass
%                    car/ext files form a specific directory.
%                'stimfileDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    stimfiles form a specific directory.
%                'pdfDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    pdf files form a specific directory.
%                'numMotionComp=1': Set to 0 if you don't want to run MLR motion comp. Set to > 1 if you want
%                    to set multiple motionComp parameters (e.g. for motionComping two sets of scans taken at
%                    different resolutions)
%                'anatFilename=[]': Set this to a filename if you want to specify a particular name for the
%                    anatomy file (e.g. 'anatFilename=myanat'), you can also set to a cell array of filenames
%                'tsenseUseMask=1': Set this to 0 if you want to suppress using the mask with tSense processing
%                'tsenseUseNoise=1': Set this to 0 if you want to suppress using the noise with tSense processing
%                'tsenseRef=-1': set to a number to specify which volume to calculate the sensitivity profile for
%                     tSense to. Default (-1) will compute the average sensitivity profile across all volumes. Set
%                     to 0 if you want the sensitivity profile computed every full volume.
%
%             First pass will sort through the specified datadir and copy
%             all of these files in the correct directories on your local
%             computer.
%
%             If you are running the second pass again and want to grab the
%             acq/respir/cardio peak bit files and the mask files from an
%             already run directory (this usually happens if you need to change
%             something in the processing and are running again from scratch)
%             You can specify getFiles, and it will copy over those files
%             on the second pass:
%
%             dofmrigru('getFiles=~/s00620120216_old')
%
function retval = dofmrigru(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14])
  help dofmrigru
  return
end

% Default arguments
global dataDir;
fidDir = [];
carextDir = [];
pdfDir = [];
stimfileDir = [];
numMotionComp = [];
movepro=[];
global epibsi;
global epibsiArgs;
global postproc;
global senseCommand;
global tsenseCommand;
getArgs(varargin,{'dataDir=/usr1/justin/data','fidDir=[]','carextDir=[]','pdfDir=[]','stimfileDir=[]','epibsi=epibsi','postproc=/usr4/local/mac_bin2/pp','tsenseCommand=/usr4/local/mac_bin2/tsense_test','senseCommand=/usr1/mauro/SenseProj/command_line/current/executables/sense_mac_intel','numMotionComp=1','movepro=0','tsense=1','dcCorrect=[]','navCorrectMag=[]','navCorrectPhase=[]','refScan=[]','getFiles=[]','anatFilename=[]','tsenseUseMask=1','tsenseUseNoise=1','notchFilter=1','tsenseRef=-1'});

% interpert the arguments for epibsiArgs
epibsiArgs = setEpibsiArgs(navCorrectMag,navCorrectPhase,dcCorrect,refScan);

% check to make sure we have the computer setup correctly to run epibsi, postproc and sense
if checkfMRISupportUnitCommands == 0,return,end

% make sure that MLR is not running
mrQuit;

% check correct directory
if strcmp(getLastDir(pwd),'Pre') || isfile(fullfile(fileparts(pwd),'mrSession.mat'))
  if askuser(sprintf('(dofmrigru) Current path is %s. Did you want to start in %s',pwd,fileparts(pwd)));
    cd('..');
  end
end

% see if this is the first preprocessing or the second one
if ~isdir('Pre')
  disp(sprintf('(dofmrigru) Running Initial dofmrigru process'));
  dofmrigru1(fidDir,carextDir,stimfileDir,pdfDir,numMotionComp,tsense,anatFilename,notchFilter);
else
  disp(sprintf('(dofmrigru) Running Second dofmrigru process'));
  dofmrigru2(movepro,tsense,getFiles,anatFilename,tsenseUseMask,tsenseUseNoise,dcCorrect,notchFilter,tsenseRef);
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
otherFiles = {'dofmrigru2.log','logbook'};
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
%%   dofmrigru2   %%
%%%%%%%%%%%%%%%%%%%%
function dofmrigru2(movepro,tsense,getFiles,anatFilename,tsenseUseMask,tsenseUseNoise,dcCorrect,notchFilter,tsenseRef)

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

% set the tsense array
if tsense == 1
  tsense = setTsense(tsense,length(epiNums));
else
  tsense = setTsense(tsense,length(epiNums));
  % check the tsense settings, to fill out the parameters if user passed in
  [tf fidList tsense] = checkTsense(fidList,epiNums,tsense);
  if ~tf && ~askuser('(dofmrigru) Continue?')
    return
  end
end

% load up tsense info if it exist
if isfile('tsense.mat')
  tsenseLoaded = load('tsense');
  % check for match
  if ~isempty(tsense) && ~all(cell2mat(tsense)==1) && ~isequal(tsense,tsenseLoaded.tsense)
    disp(sprintf('(dofmrigru2) Tsense acceleration factor from pass 1 (%s) are different from pass 2 (%s)',num2str(cell2mat(tsenseLoaded.tsense),'%i '),num2str(cell2mat(tsense),'%i ')));
    if ~askuser('Continue'),return,end
  end
  tsense = tsenseLoaded.tsense;
end

% check the tsense settings again, to make sure
[tf fidList tsense] = checkTsense(fidList,epiNums,tsense);
if ~tf && ~askuser('(dofmrigru) Continue?')
  return
end

% check for sense processing
senseProcessing = isSenseProcessing(fidList,epiNums);

% if getFiles was set, then check in the directory specified by getFiles
% for the peak files and the mask files (if any any copy them.
if ~isempty(getFiles) 
  getFiles = mlrExplicitPath(fullfile(getFiles,'Pre'));
  dispHeader(sprintf('Get files from %s',getFiles));
  if ~isdir(getFiles)
    disp(sprintf('(dofmrigru) Could not find direcotry %s from which to get files',getFiles));
    return
  end
  % now check for epis
  getFilesFidList = getFileList(getFiles,'fid');
  getFilesFidList = getFidInfo(getFilesFidList);
  getFilesFidList = sortFidList(getFilesFidList);
  % get list of epi scans
  getFilesEpiNums = getEpiScanNums(getFilesFidList);
  % check to make sure that all epi runs are processed
  [getFilesEpiNumsWithPeaks getFilesEpiNumsWithCarExt] = checkForPeaks(getFilesFidList,getFilesEpiNums);
  % now check to see if we have matching epis
  toEPI = {fidList{epiNums}};
  fromEPI = {getFilesFidList{getFilesEpiNums}};
  % check length
  if length(toEPI) ~= length(fromEPI)
    disp(sprintf('(dofmrigru) Found %i epi directories in %s, but there should be %i',length(fromEPI),getFiles,length(toEPI)));
    return
  end
  peakFileNames = {'acq.peak.bit','cardio.peak.bit','respir.peak.bit'};
  % display what bit files are going to be copied
  for i = 1:length(fromEPI)
    % check matching epi name
    if ~strcmp(fromEPI{i}.filename,toEPI{i}.filename)
      disp(sprintf('(dofmrigru) !!!! EPI filename from %s is %s, does not match %s !!!!',getFiles,fromEPI{i}.filename,toEPI{i}.filename));
      return
    end
    % look for bit files
    if sum(fromEPI{i}.hasPeakFiles) == 0
      disp(sprintf('(dofmrigru) No peak files in %s',fullfile(getFiles,fromEPI{i}.filename)));
    else
      % has some bit files, display what will be copied
      for iPeak = 1:3
	if fromEPI{i}.hasPeakFiles(iPeak)
	  disp(sprintf('copyfile %s %s',fullfile(getFiles,fromEPI{i}.filename,peakFileNames{iPeak}),toEPI{i}.fullfile));
	end
      end
    end
  end
  % look for mask files
  getFilesMaskList = getFileList(getFiles,'hdr','img','mask'); 
  for i = 1:length(getFilesMaskList)
    disp(sprintf('copyfile %s Pre',fullfile(getFiles,getFilesMaskList{i}.hdrfilename)));
    disp(sprintf('copyfile %s Pre',fullfile(getFiles,getFilesMaskList{i}.imgfilename)));
  end
  % now do it
  if askuser('Copy above files')
    for i = 1:length(fromEPI)
      % has some bit files, display what will be copied
      for iPeak = 1:3
	if fromEPI{i}.hasPeakFiles(iPeak)
	  eval(sprintf('copyfile %s %s',fullfile(getFiles,fromEPI{i}.filename,peakFileNames{iPeak}),toEPI{i}.fullfile));
	end
      end
    end
    % copy mask files
    for i = 1:length(getFilesMaskList)
      eval(sprintf('copyfile %s Pre',fullfile(getFiles,getFilesMaskList{i}.hdrfilename)));
      eval(sprintf('copyfile %s Pre',fullfile(getFiles,getFilesMaskList{i}.imgfilename)));
    end
    % now reget fid list
    fidList = getFidInfo(fidList);
    fidList = sortFidList(fidList);
  else
    return
  end
end

% get mask and noise for tsense
tSenseMaskList = [];tSenseNoiseList = [];
if ~isempty(tsense)
  % go look for masks if we are using masks
  if tsenseUseMask
    tSenseMaskList = getFileList(fiddir,'hdr','img','mask'); 
    if isempty(tSenseMaskList)
      disp(sprintf('(dofrmigru2) Could not find any mask files for tSense processing'));
    else
      % load headers for each mask
      for i = 1:length(tSenseMaskList)
	% load header
	tSenseMaskList{i}.h = mlrImageHeaderLoad(tSenseMaskList{i}.fullfile);
	tSenseMaskList{i}.dispstr = sprintf('%s [%s]',tSenseMaskList{i}.filename,num2str(tSenseMaskList{i}.h.dim,'%i '));
      end
      % display list of mask files
      dispList(tSenseMaskList,1:length(tSenseMaskList),'Mask files for tSense');
    end
  end
  % go look for noise scans if we are using noise
  if tsenseUseNoise
    tSenseNoiseList = getFileList(fiddir,'edt','epr','noise');
    if isempty(tSenseNoiseList)
      disp(sprintf('(dofrmigru2) Could not find any noise files for tSense processing'));
    else
      dispList(tSenseNoiseList,1:length(tSenseNoiseList),'Noise files for tSense');
    end
  end
end

% get mask/noise/ref for sense processing
maskList = [];noiseList = [];refList = [];
if senseProcessing
  % look for mask in nifti form
  %the mask is made manually before dofmrigru1 from noise/ref in maskdir.
  %copy the correct maskfile into Pre MANUALLY
  maskList = getFileList(fiddir,'hdr','img','mask'); 
  if isempty(maskList)
    disp(sprintf('(dofrmigru2) Could not find any mask files'));
    return
  end
  noiseList = getFileList(fiddir,'edt','epr','noise');
  if isempty(noiseList)
    disp(sprintf('(dofrmigru2) Could not find any noise files'));
    return
  end
  refList = getFileList(fiddir,'edt','epr','ref');
  if isempty(refList)
    disp(sprintf('(dofrmigru2) Could not find any ref files'));
    return
  end
end

% get anatomy scan nums
anatNums = getAnatScanNums(fidList,anatFilename);

% check to make sure that all epi runs are processed
[epiNumsWithPeaks epiNumsWithCarExt] = checkForPeaks(fidList,epiNums);

% get some info for sense processing
if senseProcessing
  senseScanNums = getSenseScanNums(fidList,epiNums);
  % get which mask goes for each sense file
  [maskNums noiseNums refNums] = chooseMaskForSenseFiles(fidList,maskList,noiseList,refList,epiNums);
else
  % no sense processing defaults
  senseScanNums = [];maskNums = [];noiseNums = [];refNums = [];
end

% get info for tsense processing
if ~isempty(tsense)
  % get the associated mask and noise nums
  [tSenseMaskNums tSenseNoiseNums] = chooseMaskForSenseFiles(fidList,tSenseMaskList,tSenseNoiseList,[],epiNums);
else
  tSenseMaskNums = [];
  tSenseNoiseNums = [];
end

% check for something to do
if isempty(epiNumsWithCarExt)
  if ~askuser('(dofmrigru) No epi scans with car/ext files. Process anyway?'),return,end
  epiNumsWithCarExt = epiNums;
% check to see if we should quit given that some of the peak files are missing
elseif length(epiNumsWithPeaks) ~= length(epiNums)
  if ~askuser('You are missing some peak files. Do you still want to continue'),return,end
end

% disp the scans
dispList(fidList,epiNumsWithCarExt,sprintf('Epi scans to process'));
dispList(fidList,anatNums,sprintf('Anatomy files'));

% display the commands for processing
processFiles(1,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect,notchFilter,tsenseRef);

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
processFiles(0,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect,notchFilter,tsenseRef);

%%%%%%%%%%%%%%%%%%%%%%
%%   processFiles   %%
%%%%%%%%%%%%%%%%%%%%%%
function processFiles(justDisplay,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro,tsense,tSenseMaskList,tSenseMaskNums,tSenseNoiseList,tSenseNoiseNums,dcCorrect,notchFilter,tsenseRef)

global postproc;
global senseCommand;
global tsenseCommand;

command = sprintf('cd Pre');
if justDisplay,disp(command),else,eval(command),end

% open the logfile
if ~justDisplay
  openLogfile('dofmrigru2.log');
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

if senseProcessing
  for i = 1:length(maskList)
    disp(sprintf('=============================================='));
    disp(sprintf('Convert masks from nifti to sdt/spr'));
    disp(sprintf('=============================================='));
    % convert mask from nifti to to sdt/spr form
    maskname = maskList{i}.filename; % still in nifti form
    command = sprintf('[d h] = cbiReadNifti(''%s'');',maskname);
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('writesdt(''%s'',d);',setext(maskname,'sdt'));
    if justDisplay,disp(command),else,eval(command),end
  end
end

if ~isempty(tsense)
  for i = 1:length(tSenseMaskList)
    disp(sprintf('=============================================='));
    disp(sprintf('Convert tSense masks from nifti to sdt/spr'));
    disp(sprintf('=============================================='));
    % convert mask from nifti to to sdt/spr form
    maskname = tSenseMaskList{i}.filename; % still in nifti form
    command = sprintf('[d h] = cbiReadNifti(''%s'');',maskname);
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('writesdt(''%s'',d);',setext(maskname,'sdt'));
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
disp(sprintf('Physiofix, sense process and convert to nifti epi files'));
disp(sprintf('=============================================='));
% physiofix the files

for i = 1:length(epiNumsWithCarExt)
  % get some names for things
  fidname = fidList{epiNumsWithCarExt(i)}.filename;
  edtname = setext(fidList{epiNumsWithCarExt(i)}.filename,'edt');
  sdtname = setext(fidList{epiNumsWithCarExt(i)}.filename,'sdt');
  hdrname = setext(fidList{epiNumsWithCarExt(i)}.filename,'hdr');
  imgname = setext(fidList{epiNumsWithCarExt(i)}.filename,'img');

  % check if we should run with pp with dc correction
  ppoptions = '';
  if dcCorrect,ppoptions = '-dc';end
  %check to see if this epi has peaks. if so run postproc with physiofix
  if any(epiNumsWithCarExt(i) == epiNumsWithPeaks)
    % run with physiofix correction
    ppoptions = sprintf('%s -physiofix',ppoptions);
  end
  
  if movepro~=0,disp('(dofmrigru) movepro not implemented - need to read fid directly using fid2nifiti and movepro');keyboard;end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %    sense processing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if any(i==senseScanNums)
    % noise, ref and mask numbers
    noiseNum = noiseNums(find(i==senseScanNums));
    refNum = refNums(find(i==senseScanNums));
    maskNum = maskNums(find(i==senseScanNums));
    %covarNum = covarNums(find(i==senseScanNums));
    maskname = setext(maskList{maskNum}.filename,'sdt');
    %covarname = covarList{covarNum}.filename;
    noisename = noiseList{noiseNum}.filename;
    refname = refList{refNum}.filename;
    % run postproc
    command = sprintf('mysystem(''%s -outtype 1 %s %s %s'');',postproc,ppoptions,fidname,edtname);
    if justDisplay,disp(command),else,eval(command),end
    % run the sense processing for this file
    command = sprintf('mysystem(''%s -data %s -full %s -noise %s -mask %s -remove -recon %s'');',senseCommand,stripext(fidname),stripext(refname),stripext(noisename),stripext(maskname),stripext(fidname));
    if justDisplay,disp(command),else,eval(command),end
    % then convert the sdt file into a nifti
    command = sprintf('[hdr] = fid2niftihdr(''%s'');',fidname);
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('data = readsdt(''%s.sdt'');',stripext(sdtname));
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('cbiWriteNifti(''%s.hdr'',data.data,hdr);',stripext(fullfile('..','Raw','TSeries',hdrname)));
    if justDisplay,disp(command),else,eval(command),end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %    tsense
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif ~isempty(tsense) && tsense{i}(1)>1
    % convert to epr
    command = sprintf('mysystem(''%s %s %s -outtype 1 %s.edt'');',postproc,ppoptions,fidname,stripext(fidname));
    if justDisplay,disp(command),else,eval(command),end
    if ~justDisplay && ~isfile(edtname)
      disp(sprintf('(dofmrigru) File %s not exists. Postproc failed?',edtname));
      keyboard
    end
    % see if there is a mask file
    if ~isempty(tSenseMaskNums) && (tSenseMaskNums(i) ~= 0)
      maskStr = sprintf('-mask %s',stripext(tSenseMaskList{tSenseMaskNums(i)}.filename));
    else
      maskStr = '';
    end
    % see if there is a noise file
    if ~isempty(tSenseNoiseNums) 
      noiseStr = sprintf('-noise %s',stripext(tSenseNoiseList{tSenseNoiseNums(i)}.filename));
    else
      noiseStr = '';
    end
    % set the reference for tsense (this usually defaults to 0 which means use each volume
    % but you can pass in the argument tsenseRef to set which volume to use for the reference
    % for caluclating the sensitivity map
    tsenseRefStr = '';
    if ~isempty(tsenseRef)
      tsenseRefStr = sprintf('-ref %i',tsenseRef);
    end
    % run tsense 
    if length(tsense{i}) > 1
      command = sprintf('mysystem(''%s -full %s -recon %s -accf %s -remove %s %s %s'');',tsenseCommand,stripext(fidname),stripext(fidname),num2str(tsense{i},'%i '),maskStr,noiseStr,tsenseRefStr);
    else
      command = sprintf('mysystem(''%s -full %s -recon %s -remove %s %s %s'');',tsenseCommand,stripext(fidname),stripext(fidname),maskStr,noiseStr,tsenseRefStr);
    end
    if justDisplay,disp(command),else,eval(command),end
    % load the created file
    command = sprintf('d = readsdt(''%s'');',sdtname);
    if justDisplay,disp(command),else,eval(command),end
    % then convert the output of tsense to a valid nifti file, by
    % pasting on the header from fid2niftihdr
    command = sprintf('[hdr info] = fid2niftihdr(''%s'',1);',fidname);
    if justDisplay,disp(command),else,eval(command),end
    % double number of volumes
    command = sprintf('hdr.dim(5) = hdr.dim(5)*%i;hdr.pixdim(5) = hdr.pixdim(5)/%i',tsense{i}(1),tsense{i}(1));
    if justDisplay,disp(command),else,eval(command),end
    % save header
    command = sprintf('cbiWriteNifti(''%s'',d.data,hdr);',fullfile('..','Raw','TSeries',hdrname));
    if justDisplay,disp(command),else,eval(command),end
  else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %    normal processing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use postproc to read the fid, run physiofix and dc correction and output to sdt
    command = sprintf('mysystem(''%s -outtype 2 %s %s %s'');',postproc,ppoptions,fidname,sdtname);
    if justDisplay,disp(command),else,eval(command),end
    % load sdt
    command = sprintf('sdt = readsdt(''%s'');',sdtname);
    if justDisplay,disp(command),else,eval(command),end
    % convert the sdt file into a nifti by getting the header info from fid2niftihdr
    command = sprintf('[hdr info] = fid2niftihdr(''%s'',1,''movepro=%f'');',fidname,movepro);
    if justDisplay,disp(command),else,eval(command),end
    % combine coils
    disp(sprintf('(dofmrigru) Take sum-of-squares of coils and remove ref scans'));
    if ~justDisplay && info.numReceivers > 1
      % reshape data to make receivers last dimension
      sdt.data = reshape(sdt.data,info.dim(1),info.dim(2),info.dim(3),info.numReceivers,info.dim(4)+info.nRefVolumes);
      % take sum of squares
      sdt.data = sqrt(sum(sdt.data(:,:,:,:,:).^2,4)),
      % squeeze out the receivers dimension - but dont use squeeze
      % cause this would squeeze out the slice dim if you only have one slice
      sdt.data = reshape(sdt.data,info.dim(1),info.dim(2),info.dim(3),info.dim(4)+info.nRefVolumes);
      % remove reference
      sdt.data = sdt.data(:,:,:,info.nRefVolumes+1:end);
    end
    % save nifti to Raw/TSeries
    command = sprintf('cbiWriteNifti(''%s'',sdt.data,hdr);',fullfile('..','Raw','TSeries',hdrname));
    if justDisplay,disp(command),else,eval(command),end
  end
end

command = sprintf('cd ..');
if justDisplay,disp(command),else,eval(command),end

if ~isempty(tsense) && notchFilter
  disp(sprintf('=============================================='));
  disp(sprintf('Run notch filtering for tSense data'));
  disp(sprintf('=============================================='));
  if ~justDisplay
    % load notch params
    notchParams = load('tSenseNotchParams');
    if ~isfield(notchParams,'params') || isempty(notchParams.params)
      disp(sprintf('(dofmrigru:processFiles) No notch settings, skipping notch filtering'));
      notchFilter = 0;
    else
      % run the notch filter.
      v = newView;
      v = tSenseNotch(v,notchParams.params);
      deleteView(v);
    end
  end
end

disp(sprintf('=============================================='));
if ~isempty(tsense) && notchFilter
  disp(sprintf('Run motion comp on notch filtered data'));
else
  disp(sprintf('Run motion comp'));
end
disp(sprintf('=============================================='));
if ~justDisplay
  motionCompFilenameNum = 1;
  motionCompFilename = sprintf('motionCompParams1.mat');
  while isfile(motionCompFilename)
    % load the params
    eval(sprintf('load %s',motionCompFilename));
    % see if we need to modify the parameters if notch filtering
    % has been run
    if ~isempty(tsense) && notchFilter
      % go into notch group to find out which scans correspond to the ones
      % initially set in the raw group
      v = newView;
      notchGroupNum = viewGet(v,'groupNum',notchParams.params.newGroupName);
      % if the notch group has been created
      if ~isempty(notchGroupNum)
	% set the group to the notch group
	v = viewSet(v,'curGroup',notchParams.params.newGroupName);
	nScans = viewGet(v,'numScans');
	% get the original scan names
	notchScan = [];
	for iScan = 1:nScans
	  thisOriginalGroup = viewGet(v,'originalGroupName',iScan);
	  thisOriginalFilename = viewGet(v,'originalFilename',iScan);
	  % check if the scan comes from single group and scan (if it doesn't then we can't do motionComp, so will skip)
	  % Also must come from Raw group
	  if iscell(thisOriginalGroup) && (length(thisOriginalGroup)==1) && iscell(thisOriginalFilename) && (length(thisOriginalFilename) == 1) && isequal(thisOriginalGroup{1},'Raw')
	    notchScan(end+1).scanNum = iScan;
	    notchScan(end).originalFilename = thisOriginalFilename{1};
	    notchScan(end).tseriesFilename = viewGet(v,'tseriesfile',iScan);
	    notchScan(end).description = viewGet(v,'description',iScan);
	  end
	end
	% now make a new structure to run the motin comp from the notch group
	if ~isempty(notchScan)
	  % set up a motion comp params structure same as the desired one but with no scans
	  notchMotionParams = params;
	  notchMotionParams.targetScans = [];
	  notchMotionParams.tseriesfiles =  {};
	  notchMotionParams.descriptions = {};
	  notchMotionParams.groupName = notchParams.params.newGroupName;
	  % for each scan if it matches a notch scan then move it into the notchMotionParams
	  for iNotchScan = 1:length(notchScan)
	    % see if it exists in the motionComp params
	    notchMatchScan = find(strcmp(stripext(notchScan(iNotchScan).originalFilename),stripext(params.tseriesfiles)));
	    if length(notchMatchScan) == 1
	      % move it to the notch params
	      notchMotionParams.targetScans(end+1) = notchScan(iNotchScan).scanNum;
	      notchMotionParams.tseriesfiles{end+1} = notchScan(iNotchScan).tseriesFilename;
	      notchMotionParams.descriptions{end+1} = sprintf('motionComp of %s',notchScan(iNotchScan).description);
	      % remove it from the regular motionComp
	      keepScans = setdiff(1:length(params.targetScans),notchMatchScan);
	      params.targetScans = params.targetScans(keepScans);
	      params.tseriesfiles = {params.tseriesfiles{keepScans}};
	      params.descriptions = {params.descriptions{keepScans}};
	    end
	  end
	  % run it if not empty
	  if length(notchMotionParams.targetScans) > 0
	    v = motionComp(v,notchMotionParams);
	    deleteView(v);
	    writeLogFile(sprintf('(dofmrigru) motionComp run on Notch for scans: %s',num2str(notchMotionParams.targetScans)))
	  end
	  if length(params.targetScans) > 0
	    % warn if there are still left over scans in params
	    dispMessage = sprintf('(dofmrigru) Not all scans have been notch filtered. Running motionComp on %s',num2str(params.targetScans));
	    disp(dispMessage);
	    writeLogFile(dispMessage);
	  end
	end
      end
    end
    % end changing of parameters of motionComp for running with notch filtered data

    % run the motion comp
    if isfield(params,'targetScans') && (length(params.targetScans)>0)
      v = newView;
      v = motionComp(v,params);
      deleteView(v);
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   chooseMaskForSenseFiles   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maskNums noiseNums refNums] = chooseMaskForSenseFiles(fidList,maskList,noiseList,refList,epiNums)
% set up variables
maskNums = []; noiseNums = []; refNums = [];
if isempty(epiNums),return,end

% no specified files
if isempty(noiseList)&&isempty(maskList)&&isempty(refList),return,end

% get the names of the mask files
if ~isempty(maskList)
  for i = 1:length(maskList)
    maskNames{i} = maskList{i}.filename;
  end
end

% get the names of the noise files
if ~isempty(noiseList)
  for i = 1:length(noiseList)
    noiseNames{i} = noiseList{i}.filename;
  end
end

% get the names of the covar files
if ~isempty(refList)
  for i = 1:length(refList)
    refNames{i} = refList{i}.filename;
  end
end

% get info from scans
for i = 1:length(epiNums)
  scanNames{i} = sprintf('%s [%s]',stripext(fidList{epiNums(i)}.filename),num2str(fidList{epiNums(i)}.info.dim,'%i '));
  if ~isempty(maskList)
    maskNamesMatch{i} = {'None' maskNames{:}};
  end
  if ~isempty(noiseList)
    noiseNamesMatch{i} = noiseNames;
  end
  if ~isempty(refList)
    refNamesMatch{i} = refNames;
  end
  startTime{i} = fidList{epiNums(i)}.startTimeStr;
  endTime{i} = fidList{epiNums(i)}.endTimeStr;
end

% check dimensions of mask
epiHasMatchingMask(1:length(epiNums)) = 0;
for i = 1:length(maskList)
  % check dimensions to be 3
  if maskList{i}.h.nDim == 3
    % compare against that of each epi scan
    for j = 1:length(epiNums)
      if isequal(fidList{epiNums(j)}.info.dim(1:3),maskList{i}.h.dim(1:3))
	% matching dimensions, promote this mask to the top of list as first choise
	maskNamesMatch{j} = putOnTopOfList(maskList{i}.filename,maskNamesMatch{j});
	epiHasMatchingMask(j) = 1;
      end
    end
  else
    disp(sprintf('(dofmrigru:chooseMaskForSenseFiles) Mask %s has strange dimensions: %s',maskList{i}.fullfile,num2str(maskList{i}.h.dim,'%i')));
  end
end

% display warning for mismatch
for i= 1:length(epiHasMatchingMask)
  if ~epiHasMatchingMask(i)
    disp(sprintf('(dofmrigru:chooseMaskForSenseFiles) EPI scan %s with dims [%s] does not match dimensions with any of the masks',fidList{epiNums(i)}.filename,num2str(fidList{epiNums(i)}.info.dim,'%i ')));
  end
end

% make an mrParamsDialog where you can select
paramsInfo{1} = {'scanNum',1,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',length(epiNums))};
paramsInfo{2} = {'scanName',scanNames,'group=scanNum','type=String','editable=0'};
paramsInfo{3} = {'startTime',startTime,'group=scanNum','type=String','editable=0'};
paramsInfo{4} = {'endTime',endTime,'group=scanNum','type=String','editable=0'};
if ~isempty(maskList)
  paramsInfo{end+1} = {'maskFile',maskNamesMatch,'group=scanNum'};
end
if ~isempty(noiseList)
  paramsInfo{end+1} = {'noiseFile',noiseNamesMatch,'group=scanNum'};
end
if ~isempty(refList)
  paramsInfo{end+1} = {'refFile',refNamesMatch,'group=scanNum'};
end

% put up the dialog
if (isempty(maskList) || (length(maskNames) == 1)) && (isempty(refList) || (length(refNames) == 1)) && (isempty(noiseList) || (length(noiseNames) == 1))
  % if there is nothing to choose (i.e. there are not multiple mask and covar, then just select the default)
  params = mrParamsDefault(paramsInfo);
else
  % otherwise have the user choose
  params = mrParamsDialog(paramsInfo,'Select which mask to use');
end
if isempty(params),return;end

% get which mask,noise and ref to use
for i = 1:length(epiNums)
  if ~isempty(maskList)
    if strcmp(params.maskFile{i},'None')
      maskNums(i) = 0;
    else
      maskNums(i) = find(strcmp(params.maskFile{i},maskNames));
    end
  end
  if ~isempty(noiseList)
    noiseNums(i) = find(strcmp(params.noiseFile{i},noiseNames));
  end
  if ~isempty(refList)
    refNums(i) = find(strcmp(params.refFile{i},refNames));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   checkForPeaks   %%
%%%%%%%%%%%%%%%%%%%%%%%
function [epiNumsWithPeaks epiNumsWithCarExt] = checkForPeaks(fidList,epiNums)

passedCheckCarExt = ones(1,length(epiNums));
passedCheckPeaks = ones(1,length(epiNums));
runForSixtySecs = zeros(1,length(epiNums));

for i = 1:length(epiNums)
  % shortcut
  fid = fidList{epiNums(i)};
  % check for car/ext
  if ~any(fid.hasCarExt)
    passedCheckCarExt(i) = 0;
    % display what is missing
    filenames = {'car','ext'};
    for j = 1:length(filenames)
      if ~fid.hasCarExt(j)
	disp(sprintf('%s: Missing %s file',fid.filename,filenames{j}));
      end
    end
  end
  % check for peak files
  if ~all(fid.hasPeakFiles)
    passedCheckPeaks(i) = 0;
    % display what is missing
    filenames = {'acq.peak.bit','cardio.peak.bit','respir.peak.bit'};
    for j = 1:length(filenames)
      if ~fid.hasPeakFiles(j)
	disp(sprintf('%s: Missing %s file (need to run peak)',fid.filename,filenames{j}));
      end
    end
  end
  % check, if it was run for more than 60 seconds, assume it is real
  if fid.info.elapsedSecs > 60
    runForSixtySecs(i) = true;
  end
end

scansThatPassedPeaks = find(passedCheckPeaks);
epiNumsWithPeaks = [];
for i = 1:length(scansThatPassedPeaks)
  epiNumsWithPeaks(end+1) = epiNums(scansThatPassedPeaks(i));
end

% if it had no car/ext, but is longer than sixty seconds, then add here
scansThatPassedCarExt = find(passedCheckCarExt|runForSixtySecs);
epiNumsWithCarExt = [];
for i = 1:length(scansThatPassedCarExt)
  epiNumsWithCarExt(end+1) = epiNums(scansThatPassedCarExt(i));
end

%%%%%%%%%%%%%%%%%%%%
%%   dofmrigru1   %%
%%%%%%%%%%%%%%%%%%%%
function dofmrigru1(fidDir,carextDir,stimfileDir,pdfDir,numMotionComp,tsense,anatFilename,notchFilter)

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
  disp(sprintf('(dofmrigru1) Could not find fid directory %s',fidDir));
  return
end

% now get info about fids
[fidList fidListArray] = getFileList(fidDir,'fid');
fidList = getFidInfo(fidList);
fidList = sortFidList(fidList);

% get list of epi scans
epiNums = getEpiScanNums(fidList);
if isempty(epiNums)
  disp(sprintf('(dofmrigru1) Could not find any epi files in %s',fidDir));
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
if ~tf && ~askuser('(dofmrigru) Continue?')
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
  disp(sprintf('(dofmrigru1) Could not find any non-raw anatomies in %s',fidDir));
end

% get sense noise/ref scans
[senseNoiseNums senseRefNums] = getSenseNums(fidList);
if isempty(senseRefNums) && senseProcessing
  disp(sprintf('(dofmrigru1) Could not find any sense ref files in %s',fidDir));
  return
end
if isempty(senseNoiseNums) && (senseProcessing || ~isempty(tsense))
  disp(sprintf('(dofmrigru1) Could not find any sense noise files in %s',fidDir));
  return
end

% find car/ext files
carList = getFileList(carextDir,'car','ext');
carList = getCarInfo(carList);
if isempty(carList)
  disp(sprintf('(dofmrigru1) Could not find any car/ext files in %s',carextDir));
end

% find pdf files
pdfList = getFileList(pdfDir,'pdf');

% find edf files
edffileList = getFileList(stimfileDir,'edf');

% get stimfile list
stimfileList = getFileList(stimfileDir,'mat');
stimfileList = setStimfileListDispStr(stimfileList);

% display what we found
dispList(fidList,epiNums,sprintf('Epi scans: %s',fidDir));
dispList(fidList,anatNums,sprintf('Anatomy scans: %s',fidDir));
if senseProcessing
  dispList(fidList,senseRefNums,sprintf('Sense reference scan: %s',fidDir));
end
if senseProcessing || ~isempty(tsense)
  dispList(fidList,senseNoiseNums,sprintf('Noise scan: %s',fidDir));
end
dispList(pdfList,nan,sprintf('PDF files: %s',pdfDir));
dispList(stimfileList,nan,sprintf('Stimfiles: %s',stimfileDir));
dispList(edffileList,nan,sprintf('EDFfiles: %s',stimfileDir));
dispList(carList,nan,sprintf('Car/Ext files: %s',carextDir));

% go find the matching car files for each scan
if ~isempty(carList)
  [carMatchNum carList] = getCarMatch(carList,fidList,epiNums);
  if isempty(carMatchNum),return,end
else
  carMatchNum = [];
end

% setup directories etc.
doMoveFiles(1,fidList,carList,pdfList,stimfileList,edffileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense);

% make mask dir
if senseProcessing,makeMaskDir(1,fidList,anatNums,senseRefNums);end

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
fidList = doMoveFiles(0,fidList,carList,pdfList,stimfileList,edffileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense);

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
disp(sprintf('(dofmrigru1) Setup mrInit for your directory'));
mrInit([],[],sprintf('subject=%s',subjectID),'makeReadme=0');

% set some info in the auxParams about the scan
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
	  framePeriod = framePeriod/volTrigRatio(iScan);
	  % round to nearest 1/1000 of a second
	  framePeriod = round(framePeriod*1000)/1000;
	  disp(sprintf('(dofmrigru) Frame period as recorded in stimfile is: %0.3f',framePeriod));
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

% set up notch filter parameters
if ~isempty(tsense) && notchFilter
  v = newView;
  if ~isempty(v);
    [v params] = tSenseNotch(v,[],'justGetParams=1');
    deleteView(v);
    save('tSenseNotchParams','params');
  end
end

% set up motion comp parameters
for i = 1:numMotionComp
  v = newView;
  if ~isempty(v)
    % tell user that motion comp will be run from notch filter
    if ~isempty(tsense) && notchFilter
      disp(sprintf('(dofmrigru) Note that motion comp will be run from the Notch group (i.e. after notch filtering has been processed for tsense data). This is because we dont want the motionComp to get confused by the tsense artifacts that are present in the time series before notch filtering has been run.'));
    end
    [v params] = motionComp(v,[],'justGetParams=1');
    deleteView(v);
    eval(sprintf('save motionCompParams%i params',i));
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    makeMaskDir    %
%%%%%%%%%%%%%%%%%%%%%
function makeMaskDir(justDisplay,fidList,anatNums,senseRefNums)

disp(sprintf('=============================================='));
disp(sprintf('Make empty MLR directory for mask creation    '));
disp(sprintf('=============================================='));

% process ref and anatomy fid to make the mask
command = sprintf('cd Pre;'); myeval(command,justDisplay);% move into Pre
global maskdir;
% location of masks
% make new empty mask directory in Pre
if ~justDisplay
  makeEmptyMLRDir('Mask','description=Empty session for making a mask','defaultParams=1');
else
  disp('makeEmptyMLRDisplay');
end

maskdir = fullfile(pwd,'Mask');

for i = 1:length(senseRefNums)
  disp(sprintf('=============================================='));
  disp(sprintf('Convert ref to nifti and copy into maskdir'));
  disp(sprintf('=============================================='));
  % convert ref to nifti
  if ~justDisplay
    fid2nifti(fidList{senseRefNums(i)}.filename,setext(fidList{senseRefNums(i)}.filename,'hdr'));
  end
  
  % move into emptyMLR directory for mask
  command = sprintf('movefile %s %s;',setext(fidList{senseRefNums(i)}.filename,'hdr'),fullfile(maskdir,'Anatomy'));
  myeval(command,justDisplay);
  command = sprintf('movefile %s %s;',setext(fidList{senseRefNums(i)}.filename,'img'),fullfile(maskdir,'Anatomy'));
  myeval(command,justDisplay);

end  
% move out of Pre
command = sprintf('cd ..;'); myeval(command,justDisplay);

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
function fidList = doMoveFiles(justDisplay,fidList,carList,pdfList,stimfileList,edffileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums,tsense)

% some command names
global epibsi;
global epibsiArgs;
global postproc;

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
  openLogfile('dofmrigru.log');
  cd('..');
  % write info about various scans
  dispList(fidList,epiNums,'Epi scans',true);
  dispList(fidList,anatNums,'Anatomy scans',true);
  if ~isempty(senseRefNums) dispList(fidList,senseRefNums,'Sense reference scan',true);end
  if ~isempty(senseNoiseNums) dispList(fidList,senseNoiseNums,'Sense noise scan',true);end
  dispList(carList,nan,'Car/Ext files',true);
  dispList(pdfList,nan,'PDF files',true);
  dispList(stimfileList,nan,'Stimfiles',true);
  dispList(edffileList,nan,'EDFfiles',true);
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
disp(sprintf('Copying edf files'));
disp(sprintf('=============================================='));

% move edffiles
if ~isempty(edffileList)
    for i = 1:length(edffileList)
        command = sprintf('copyfile %s %s',edffileList{i}.fullfile,fullfile('Etc',edffileList{i}.filename));
        if justDisplay,disp(command),else,eval(command);disp(command);,end
    end
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
disp(sprintf('Copying car/ext files'));
disp(sprintf('=============================================='));

% move car/ext files into fid directories
for i = 1:length(carMatchNum)
  if carMatchNum(i) ~= -1
    command = sprintf('copyfile %s %s',carList{carMatchNum(i)}.fullfile,fullfile('Pre',fidList{epiNums(i)}.filename));
    if justDisplay,disp(command),else,eval(command),end
    if isfield(carList{carMatchNum(i)},'extfilename')
      disp('ext not found. not being copied...');
      command = sprintf('copyfile %s %s',fullfile(carList{carMatchNum(i)}.path,carList{carMatchNum(i)}.extfilename),fullfile('Pre',fidList{epiNums(i)}.filename));
    end
    if justDisplay,disp(command),else,eval(command);,end
  end
end

disp(sprintf('=============================================='));
disp(sprintf('Processing epi files'));
disp(sprintf('=============================================='));

command = sprintf('cd Pre');
if justDisplay,disp(command),else,eval(command),end

% run epibsi on all scans and sense noise/ref 
allScanNums = [epiNums senseNoiseNums senseRefNums];
for i = 1:length(allScanNums)
  command = sprintf('status = mysystem(''%s %s %s'');',epibsi,fidList{allScanNums(i)}.filename,epibsiArgs);
  if justDisplay
    disp(command)
  else
    eval(command);
    if status > 1
      disp(sprintf('(dofmrigru) %s with args %s generated error %i when running',epibsi,epibsiArgs,status));
      keyboard
    end    
  end
end

% convert sense noise/ref into edt files
senseNums = [senseNoiseNums senseRefNums];
for i = 1:length(senseNums)
  %make fid to edt
  command = sprintf('mysystem(''%s -intype 0 -outtype 1 %s %s'');',postproc,setext(fidList{senseNums(i)}.filename,'fid'),setext(fidList{senseNums(i)}.filename,'edt'));
  if justDisplay,disp(command),else,eval(command);,end
end

command = sprintf('cd ..');
if justDisplay,disp(command),else,eval(command);,end

disp(sprintf('=============================================='));
disp(sprintf('Making temporary (empty) nifti files for running mrInit'));
disp(sprintf('=============================================='));
for i = 1:length(epiNums)
  srcName = fullfile('Pre',fidList{epiNums(i)}.filename);
  destName = fullfile('Raw/TSeries',setext(fidList{epiNums(i)}.filename,'hdr'));
  if ~justDisplay
    % make a nifti header
    [h info] = fid2niftihdr(srcName);
    % now, grab the data
    if ~isempty(tsense) && (tsense{i}(1) > 1)
      % load the file
      d = fid2nifti(srcName,1);
      % now we have to create data as if it has gone through
      % tsense processing
      d = repmat(d,[1 1 1 tsense{i}(1)]);
      % fix up the header
      h.pixdim(5) = h.pixdim(5)/tsense{i}(1);
      h.dim(5) = h.dim(5)*tsense{i}(1);
    elseif info.accFactor == 1
      d = fid2nifti(srcName,1);
      % for a sense file we are going to have to grab the reference scan data
    else
      d = fid2nifti(fullfile('Pre',fidList{senseRefNums(1)}.filename),1);
      % now make the data info a single volume
      d = mean(d(:,:,:,2:end),4);
      % and replicate for the correct number of frames
      d = repmat(d,[1 1 1 h.dim(5)]);
    end
    % draw x's into image to indicate they are only temporary
    xSize = 3;
    x = zeros(min(xSize,size(d,1)),min(xSize,size(d,2)));
    for i = 1:xSize
      x(i,i) = 1;
      x(i,xSize-i+1) = 1;
    end
    % draw x, scaling to each slice, each volume
    for iSlice = 1:size(d,3)
      for iVol = 2:size(d,4)
	% get min and max
	minVal = min(min(d(:,:,iSlice,iVol)));
	maxVal = max(max(d(:,:,iSlice,iVol)));
	d(1:xSize,1:xSize,iSlice,iVol) = x*(maxVal-minVal)+minVal;
      end
    end
    % write as a nifti file
    cbiWriteNifti(destName,d,h);
  else
    disp(sprintf('Make temporary nifti file for %s in %s',srcName,destName));
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

%%%%%%%%%%%%%%%%%%%%%
%%   getCarmatch   %%
%%%%%%%%%%%%%%%%%%%%%
function [carMatchNum carList] = getCarMatch(carList,fidList,epiNums);

carMatchNum = [];

% get the names of the car files
for i = 1:length(carList)
  carNames{i} = carList{i}.dispstr;
end
carNames{end+1} = 'None';
noneNum = length(carNames);

% get info from scans
for i = 1:length(epiNums)
  scanNames{i} = stripext(fidList{epiNums(i)}.filename);
  if i <= length(carNames)
    carMatchNames{i} = putOnTopOfList(carNames{i},carNames);
  else
    carMatchNames{i} = putOnTopOfList('None',carNames);
  end
  startTime{i} = fidList{epiNums(i)}.startTimeStr;
  endTime{i} = fidList{epiNums(i)}.endTimeStr;
end

% make an mrParamsDialog where you can select
paramsInfo{1} = {'scanNum',1,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',length(epiNums))};
paramsInfo{2} = {'scanName',scanNames,'group=scanNum','type=String','editable=0'};
paramsInfo{3} = {'startTime',startTime,'group=scanNum','type=String','editable=0'};
paramsInfo{4} = {'endTime',endTime,'group=scanNum','type=String','editable=0'};
paramsInfo{5} = {'carFile',carMatchNames,'group=scanNum'};

% put out the dialog
params = mrParamsDialog(paramsInfo,'Match the car file with the scan');
if isempty(params),return;end

% now get the matching numbers
for i = 1:length(params.carFile)
  carMatchNum(i) = find(strcmp(params.carFile{i},carNames));
end
carMatchNum(carMatchNum==noneNum) = -1;

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
%%%%%%%%%%%%%%%%%%%%
%%   getCarInfo   %%
%%%%%%%%%%%%%%%%%%%%
function carList = getCarInfo(carList)

% channel which has the stimulus trigger trace on it
trigChannel = 14;

disppercent(-inf,'(dofmrigru1) Get car info');
for i = 1:length(carList)
  disppercent(i/length(carList));
  carList{i}.car = readcar(carList{i}.fullfile);
  if isfield(carList{i},'extfilename')
    disp('ext not found. skipping...');
    carList{i}.ext = readext(fullfile(carList{i}.path,carList{i}.extfilename));
  end
  
  % make the display string
  dispstr = sprintf('%s',carList{i}.filename);

  % get number of triggers
  if isequal(exist('getedges'),2) && isfield(carList{i},'car')  && isfield(carList{i}.car,'channels') && (size(carList{i}.car.channels,1) >= trigChannel)
    numTrig = getedges(carList{i}.car.channels(trigChannel,:),max(carList{i}.car.channels(trigChannel,:))/2);
    dispstr = sprintf('%s %i',dispstr,length(setdiff([numTrig.rising numTrig.falling],1)));
  end
  
  % get save date
  if isfield(carList{i},'date')
    dispstr = sprintf('%s (%s end)',dispstr,carList{i}.date);
  end
  % get number of acq pulses
  if isequal(exist('getedges'),2) && isfield(carList{i},'car')  && isfield(carList{i}.car,'acq')
    numAcq = getedges(carList{i}.car.acq,0.5);
    dispstr = sprintf('%s %i',dispstr,numAcq.n);
  end
  % add this to the carList so it gets printed out later
  carList{i}.dispstr = dispstr;
end

disppercent(inf);
%%%%%%%%%%%%%%%%%%%%
%%   getFidInfo   %%
%%%%%%%%%%%%%%%%%%%%
function fidList = getFidInfo(fidList)

nVols = [];
disppercent(-inf,'(dofmrigru1) Getting fid info');
for i = 1:length(fidList)
  disppercent(i/length(fidList));
  % read the procpar
  fidList{i}.procpar = readprocpar(fidList{i}.fullfile);
  % get info using fid2xform if we can
  if isfield(fidList{i}.procpar,'ni')
    [fidList{i}.xform fidList{i}.info] = fid2xform(fidList{i}.procpar,-1);
  else
    fidList{i}.xform = [];fidList{i}.info = [];
  end
  % load the logfile
  log = 0;
  if isfile(fullfile(fidList{i}.fullfile,'log'))
    f = fopen(fullfile(fidList{i}.fullfile,'log'));
    thisline = fgets(f);
    [fidList{i}.startTime fidList{i}.startTimeStr] = getDatenumFromLogLine(thisline);
    thisline = fgets(f);
    [fidList{i}.endTime fidList{i}.endTimeStr] = getDatenumFromLogLine(thisline);
    fclose(f);
    log = 1;
  else
    disp(sprintf('(dofmrigru:getFidInfo) Missing log in %s',fidList{i}.filename));
  end
  % make the display string
  fidList{i}.dispstr = sprintf('%s: ',fidList{i}.filename);
  if isfield(fidList{i},'startTimeStr') && isfield(fidList{i},'endTimeStr')
    fidList{i}.dispstr = sprintf('%s%s->%s',fidList{i}.dispstr,fidList{i}.startTimeStr,fidList{i}.endTimeStr);
  else
    if ~isfield(fidList{i},'startTime'),fidList{i}.startTime = [];end
    fidList{i}.startTimeStr = '';
  end
  if ~isempty(fidList{i}.info)
    if fidList{i}.info.isepi && fidList{i}.info.compressedFid
      % unprocessed epi (the phase encodes are set wrong)
      fidList{i}.dispstr = sprintf('%s [%i ? %i %i] [%0.1f ? %0.1f]',fidList{i}.dispstr,fidList{i}.info.dim(1),fidList{i}.info.dim(3),fidList{i}.info.dim(4),fidList{i}.info.voxsize(1),fidList{i}.info.voxsize(3));
    else
      % otherwise
      fidList{i}.dispstr = sprintf('%s [%i %i %i %i] [%0.1f %0.1f %0.1f]',fidList{i}.dispstr,fidList{i}.info.dim(1),fidList{i}.info.dim(2),fidList{i}.info.dim(3),fidList{i}.info.dim(4),fidList{i}.info.voxsize(1),fidList{i}.info.voxsize(2),fidList{i}.info.voxsize(3));
    end
    if ~isempty(fidList{i}.info.accFactor)
      fidList{i}.dispstr = sprintf('%s x%i',fidList{i}.dispstr,fidList{i}.info.accFactor);
    end
  end
  fidList{i}.hasCarExt = [0 0];
  fidList{i}.hasPeakFiles = [0 0 0];
  if ~isempty(dir(fullfile(fidList{i}.fullfile,'*.car')))
    fidList{i}.hasCarExt(1) = 1;
  end
  if ~isempty(dir(fullfile(fidList{i}.fullfile,'*.ext')))
    fidList{i}.hasCarExt(2) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'acq.peak.bit'))
    fidList{i}.hasPeakFiles(1) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'cardio.peak.bit'))
    fidList{i}.hasPeakFiles(2) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'respir.peak.bit'))
    fidList{i}.hasPeakFiles(3) = 1;
  end
end
disppercent(inf);

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

%%%%%%%%%%%%%%%%%%%%%%
%%   getSenseNums   %%
%%%%%%%%%%%%%%%%%%%%%%
function [senseNoiseNums senseRefNums] = getSenseNums(fidList)

senseNoiseNums = [];senseRefNums = [];

for i = 1:length(fidList)
  if ~isempty(fidList{i}.info)
    % make sure this is not a 3d acq or a sense
    if ~fidList{i}.info.acq3d && (fidList{i}.info.accFactor == 1)
      % make sure it's not raw
      if fidList{i}.info.dim(3) == length(fidList{i}.info.pss)
        if ~isempty(strfind(fidList{i}.filename,'noise'))
	  senseNoiseNums(end+1) = i;
        elseif ~isempty(strfind(fidList{i}.filename,'ref'))
	  senseRefNums(end+1) = i;
        end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getSenseScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function senseNums = getSenseScanNums(fidList,epiNums)

senseNums = [];
for i = 1:length(epiNums)
  if ~isempty(fidList{epiNums(i)}.info)
    if ~isempty(fidList{epiNums(i)}.info.accFactor)
      senseNums(end+1) = i;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%   getEpiScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function epiNums = getEpiScanNums(fidList,nVolsCutoff)

if ieNotDefined('nVolsCutoff'),nVolsCutoff = 3;end

epiNums = [];
for i = 1:length(fidList)
  if ~isempty(fidList{i}.info)
    % check that it is an epi and that it has more than the cutoff number of volumes
    % cutoff is set just to discard any false-starts
    if fidList{i}.info.isepi && (fidList{i}.info.dim(4) > nVolsCutoff)
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
	      %disp(sprintf('(dofmrigru1:getFileList) No matching %s file for %s',matchExtList{j},dirList(i).name));
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
  disp(sprintf('(dofmrigru) Could not open logfile %s',filename))
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
function retval = checkfMRISupportUnitCommands

global epibsi;
global postproc;
global senseCommand;
global tsenseCommand;

retval = 1;

% commands to check
commandNames = {'epibsi','postproc','senseCommand','tsenseCommand'};
helpFlag = {'','-help','',''};
for i = 1:length(commandNames)
  % suse which to tell if we have the command
  [commandStatus commandRetval] = system(sprintf('which %s',eval(commandNames{i})));
  % check for commandStatus error
  if commandStatus~=0
    disp(sprintf('(dofmrigru) Could not find %s command: %s',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmrigru for help setting up your computer'));
    if strcmp(commandNames{i},'senseCommand')
      % just warn
      continue;
    end
    retval = 0;
    return
  end
  % run the command to see what happens
  [commandStatus commandRetval] = system(sprintf('%s %s',eval(commandNames{i}),helpFlag{i}));
  % check for commandStatus error
  if commandStatus>1
    disp(commandRetval);
    disp(sprintf('(dofmrigru) Found %s command: %s, but could not run (possibly missing fink library?)',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmrigru for help setting up your computer'));
    retval = 0;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%
% isSenseProcessing %
%%%%%%%%%%%%%%%%%%%%%
function senseProcessing = isSenseProcessing(fidList,epiNums)
senseProcessing = 0;

for i = epiNums
  if isfield(fidList{i}.info,'accFactor') && (fidList{i}.info.accFactor>1)
    senseProcessing = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    getVolTrigRatio    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function volTrigRatio = getVolTrigRatio(fidList,tsenseAcc)

% default return value
volTrigRatio = 1;

% if no acceleration then just set this to one so the code below
% calculates volTrigRation appropriately
if isempty(tsenseAcc)||(tsenseAcc==0),tsenseAcc = 1;end

% get number of shots
numshots = fidList.procpar.numshots;

if ~isfield(fidList.procpar,'image')
  disp(sprintf('(dofmrigru:getVolTrigRatio) !!! Could not find image field in procpar for scan %s - assuming one trigger per volume !!!',fidList.filename));
else
  % get the "image" field, which contains how the triggers work for each volume
  procparImage = fidList.procpar.image;
  % check its length
  if length(procparImage)<fidList.info.dim(4)
    disp(sprintf('(dofmrigru:getVolTrigRatio) !!! Image field for scan %s does not have enough volumes (%i but should be > %i) !!!',fidList.filename,length(procparImage),fidList.info.dim(4)));
  end
  % find unique values other than 0 (which is the steady-state image
  procparImage = unique(procparImage(procparImage~=0));
  % in some cases you can have 1 followed by 2 - this triggers only on the first volume of the run
  % so we handle that case here.
  if isequal(procparImage,[1 2])
    % go back to original image setting
    procparImage = fidList.procpar.image;
    % remove the zeros
    procparImage = procparImage(procparImage>0);
    % find the distance between every 1 (that is when there was a volume acq)
    volPerTrig = unique(diff(find(procparImage == 1)));
    % if there is more than one value here it means that there is not a fixed
    % ratio of triggers to volumes which will be a problem. So provide feedback
    if length(volPerTrig) ~= 1
      disp(sprintf('(dofmrigru:getVolTrigRatio) !!! Scan %s has non-unique number of triggers per volume (volPerTrig: %s)!!!',fidList.filename,num2str(volPerTrig,'%i ')))
      return
    end
    % set the volTrigRatio accordingly
    volTrigRatio = volPerTrig*tsenseAcc;
    % otherwise these should all be the same, and either 1, 3 or 4
  elseif (length(procparImage) > 1) || ~any(procparImage == [1 3 4])
    disp(sprintf('(dofmrigru:getVolTrigRatio) !!! Image field for scan %s has strange values in it (%s) !!!',fidList.filename,num2str(procparImage,'%i ')));
  elseif procparImage == 3
    volTrigRatio = tsenseAcc/numshots;
  elseif any(procparImage == [4 1])
    volTrigRatio = tsenseAcc;
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    checkTsense    %
%%%%%%%%%%%%%%%%%%%%%
function [tf fidList tsense volTrigRatio] = checkTsense(fidList,epiNums,tsense)

tf = true;
volTrigRatio = [];
if ~isempty(tsense)
  for iEPI = 1:length(epiNums)
    % get some things about the scan
    numshots = fidList{epiNums(iEPI)}.procpar.numshots;
    ilts = fidList{epiNums(iEPI)}.procpar.ilts;
    % if set to use tsense, get its tsense acceleration
    if isscalar(tsense{iEPI})
      if isequal(tsense{iEPI},1)
	% make sure it is not SENSE
	if ~(isfield(fidList{epiNums(iEPI)}.info,'accFactor') && (fidList{epiNums(iEPI)}.info.accFactor>1))
	  % calculate optimal tsense acceleration
	  tsense{iEPI} = numshots/ilts;
	end
      elseif ~isequal(tsense{iEPI},0)
	if tsense{iEPI} ~= numshots/ilts;
	  % check to make sure tsense acceleration divides nicely into numshots
	  if rem(numshots,tsense{iEPI}) ~= 0
	    disp(sprintf('(dofmrigru) tSense acceleration set to %i which does not divide evenly into number of shots'));
	    keyboard
	  end
	  % set the number shots
	  tsense{iEPI} = [tsense{iEPI} getTSenseShotOrder(fidList{epiNums(iEPI)},tsense{iEPI})];
	  disp(sprintf('(dofmrigru) tSense acceleration set to %i even though optimal is %i. Setting shot order to %s',tsense{iEPI}(1),numshots/ilts,num2str(tsense{iEPI}(2:end),'%i ')));
	end
      end
    % if set to sepecify num shots and shot order, 
    % then check to make sure that these are
    % set to sensible values
    else
      % first get accerleation
      tSenseAcc = tsense{iEPI}(1);
      % now make sure we have enough arguments to speify the shot order
      if length(tsense{iEPI}) == 1
	tsense{iEPI} = [tsense{iEPI} 1:(numshots/tsense{iEPI})];
	disp(sprintf('(dofmrigru) tSense acceleration setting shot order to %s',num2str(tsense{iEPI}(2:end),'%i ')));
      elseif ~isequal(1:tSenseAcc,sort(tsense{iEPI}(2:end)))
	disp(sprintf('(dofmrigru) !!! tSense acceleration must specify shot order for all shots. Resetting to default acceleration of %i !!!',numshots/ilts))
	tsense{iEPI} = numshots/ilts;
	tf = false;
      else
	% otherwise check whether this is optimal or not
	if tSenseAcc ~= numshots/ilts
	  disp(sprintf('(dofmrigru) Using tSense acceleration factor of %i which is not optimal for numshots: %i ilts: %i',tSenseAcc,numshots,ilts));
	end
      end
    end
    if tsense{iEPI}(1)
      % set the display string
      fidList{epiNums(iEPI)}.dispstr = sprintf('%s tSense: %s tr=%s nShots=%i nVols=%i',fidList{epiNums(iEPI)}.dispstr,num2str(tsense{iEPI},'%i '),mlrnum2str(fidList{epiNums(iEPI)}.info.tr/tsense{iEPI}(1)),numshots/tsense{iEPI}(1),fidList{epiNums(iEPI)}.info.dim(4)*tsense{iEPI}(1));
    end
    % figure out if we had a trigger every shot or not
    volTrigRatio(iEPI) = getVolTrigRatio(fidList{epiNums(iEPI)},tsense{iEPI}(1));
  end
end

%%%%%%%%%%%%%%%%%%%
%    setTsense    %
%%%%%%%%%%%%%%%%%%%
function tsense = setTsense(tsense,nEPI)

if ischar(tsense)
  disp(sprintf('(dofmrigru) tSense input cannot be a string.'));
  keyboard
end

% if this looks like acceleration factor plus shot order
% then make into a cell array so that each epi scan
% will use these parameters
if ~iscell(tsense) && isequal(1:tsense(1),sort(tsense(2:end)))
  tsense = {tsense};
end

if ~iscell(tsense)
  % set values to the ones in the array
  % if there are not enough values, for example
  % there is only one value, set to the last
  % value in the array - thus for example
  % if the value tsense is a scalar will set
  % all values to that scalar
  accFactor = tsense;
  tsense = {};
  for i = 1:nEPI
    tsense{i} = accFactor(min(i,length(accFactor)));
  end
elseif iscell(tsense)
  for i = length(tsense)+1:nEPI
    tsense{i} = tsense{end};
  end
end  

if all(cell2mat(tsense) == 0)
  % no tsense
  tsense = [];
end

%%%%%%%%%%%%%%%%%%%%
%    epibsiArgs    %
%%%%%%%%%%%%%%%%%%%%
function epibsiArgs = setEpibsiArgs(navCorrectMag,navCorrectPhase,dcCorrect,refScan)

% default settings
epibsiArgs = '';

% defaults
defaultNavCorrectMag = 1;
defaultNavCorrectPhase = 1;
defaultDcCorrect = 0;
defaultRefScan = 1;

% default arguments
if isempty(navCorrectMag) && isempty(navCorrectPhase) && isempty(dcCorrect) && isempty(refScan)
  return
end

% check navCorrectMag
if isempty(navCorrectMag) navCorrectMag = defaultNavCorrectMag;end
if ~any(navCorrectMag == [0 1])
  disp(sprintf('(dofmrigru) navCorrectMag should be either 0 or 1. Using default'));
  navCorrectMag = defaultNavCorrectMag;
end

% check navCorrectPhase
if isempty(navCorrectPhase) navCorrectPhase = defaultNavCorrectPhase;end
if ~any(navCorrectPhase == [0 1])
  disp(sprintf('(dofmrigru) navCorrectPhase should be either 0 or 1. Using default'));
  navCorrectPhase = defaultNavCorrectPhase;
end

% check dcCorrect
if isempty(dcCorrect) dcCorrect = defaultDcCorrect;end
if ~any(dcCorrect == [0 1])
  disp(sprintf('(dofmrigru) dcCorrect should be either 0 or 1. Using default'));
  dcCorrect = defaultDcCorrect;
end

% check refScan
if isempty(refScan) refScan = defaultRefScan;end
if refScan == 0, refScan = -1;end
if ~any(refScan == [-1 1])
  disp(sprintf('(dofmrigru) refScan should be either -1 or 1. Using default'));
  refScan = defaultRefScan;
end

epibsiArgs = sprintf('%i %i %i 0 %i -1',navCorrectMag, navCorrectPhase, dcCorrect,refScan);
disp(sprintf('(dofmrigru) epibsiArgs set to: %s (mag correct, phase correct, dc correct, print nav data, ref scan, data type)',epibsiArgs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getTSenseShotOrder    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shotOrderForAccFactor = getTSenseShotOrder(fidInfo,accFactor)

% get pelist
pelist = fidInfo.info.procpar.pelist;
% get num shots
numShots = fidInfo.info.procpar.numshots;

% find the most negative values in the pelist
[pelist shotOrder] = sort(pelist);

% get the ordering of the most negative values
[temp shotOrder] = sort(shotOrder(1:numShots));

% now shotOrder should have the ordering of the shots from most negative
% line of k-space to most positive

% next we group into the number of shots in each accelerated image
% and find the ordering of the shots
shotsPerAcceleratedImage = numShots/accFactor;
for i = 1:accFactor
  shotOrderForAccFactor(i) = min(shotOrder((i-1)*shotsPerAcceleratedImage+1:i*shotsPerAcceleratedImage));
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

