% mlrFidInfo.m
%
%        $Id:$ 
%      usage: mlrFidInfo(v)
%         by: justin gardner
%       date: 12/12/11
%    purpose: Displays fid info
%
function retval = mlrFidInfo(v)

% check arguments
if ~any(nargin == [1])
  help mlrFidInfo
  return
end

% get the fid name
fidFilename = viewGet(v,'auxParam','fidFilename');
if isempty(fidFilename)
  disp(sprintf('(mlrFidInfo) No fid filename associated with this scan.'));
  return
end

if iscell(fidFilename)
  disp(sprintf('(mlrFidInfo) Have not implemented how to handle cell arrays yet'));
else
  fidInfo = viewGet(v,'fidInfo');
  fidInfo = fidInfo{1};

  % get volTrigRatio
  volTrigRatio = viewGet(v,'auxParam','volTrigRatio');

  % get tSense
  tSense = viewGet(v,'auxParam','tSense');

  % get/set fields to skip
  skipFields = {'procpar','startDatevec','endDatevec','elapsedSecs','offset','originOffset','pss','fidname','sliceOrder'};
  fidInfoFields = fieldnames(fidInfo);
  for i = 1:length(fidInfoFields)
    if iscell(fidInfoFields{i})
      skipFields{end+1} = fidInfoFields{i};
    end
  end
  
  % now make paramsInfo
  paramsInfo = {};
  paramsInfo{1} = {'fidFilename',getLastDir(fidInfo.fidname),'editable=0'};
  paramsInfo{end+1} = {'volTrigRatio',volTrigRatio,'editable=0','Number of volumes per each acq trigger'};
  paramsInfo{end+1} = {'tSense',tSense,'editable=0','tSense acceleration factor'};
  % get slice Order/pss as string
  paramsInfo{end+1} = {'sliceOrder',num2str(fidInfo.sliceOrder,'%i '),'editable=0','Order that slices were acquired'};
  paramsInfo{end+1} = {'pss',mlrnum2str(fidInfo.pss,'compact=1'),'editable=0'};
  for i = 1:length(fidInfoFields)
    if ~strcmp(fidInfoFields{i},skipFields)
      paramsInfo{end+1} = {fidInfoFields{i},fidInfo.(fidInfoFields{i}),'editable=0'};
    end
  end

  
  % get some parameters form procpar
  getFromProcpar = {'ss','petable','image','te','fliplist'};
  if isfield(fidInfo,'procpar')
    for iProcparField = 1:length(getFromProcpar)
      if isfield(fidInfo.procpar,getFromProcpar{iProcparField})
	% get the field value from procpar
	fieldVal = fidInfo.procpar.(getFromProcpar{iProcparField});
	% if it is a cell of length one, remove the cell
	if iscell(fieldVal) && (length(fieldVal) == 1)
	  fieldVal = fieldVal{1};
	end
	% for 'image' variable, do some special processing to display 
	if strcmp(getFromProcpar{iProcparField},'image')
	  % find unique values to display
	  valsToDisp = find([1 diff(fieldVal)]);
	  valsToDispNums = diff([valsToDisp length(fieldVal)+1]);
	  dispVal = '[';
	  for iVal = 1:length(valsToDisp)
	    if valsToDispNums(iVal) < 2
	      dispVal = sprintf('%s %i',dispVal,fieldVal(valsToDisp(iVal)));
	    else
	      dispVal = sprintf('%s %ix%i',dispVal,fieldVal(valsToDisp(iVal)),valsToDispNums(iVal));
	    end
	  end
	  fieldVal = sprintf('%s]',dispVal);
	end
	% and add it as a parameter to display
	paramsInfo{end+1} = {getFromProcpar{iProcparField},fieldVal,'editable=0'};
      end
    end
  end

  
  mrParamsDialog(paramsInfo,sprintf('FID Info for %s',fidInfo.fidname));
end  




