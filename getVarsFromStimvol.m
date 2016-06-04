
%getVarsFromStimvol.m
%
%
% author: steeve laquitaine
%   date: 160112
%purpose: retrieve variables associated with instances using stimvols 
%         (from getStimvol)
%
%usage:
%
%e.g., 
%     
%     %set params
%     o.sessPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%     o.myGroup  = 'Concatenation';
%     o.myScan     = 1;
%     o.taskNum    = 2;
%     o.phaseNum   = 1;
%     o.segmentNum = 2;
%     variables = {'myRandomDir','myRandomCoh'};
%       
%     %get stimvols you want variables of
%     cd(o.sessPath)
%     v = mrLoadRet([]);
%     v = viewSet(v,'curGroup',o.myGroup,'curScan',o.myScan);;
%     stimvols = getStimvol(v,'myRandomCoh','taskNum',o.taskNum,'phaseNum',...
%           o.phaseNum,'segmentNum',o.segmentNum);
%       
%     %get variable values associated with stimvols
%     [vars,o,db] = getVarsFromStimvol(stimvols,variables,o)
% 
%
%inputs: 
%     stimvols: input stimvols
%    variables: variable values we want for stimvols
%            o: fmri session parameters
%
%output:           
%    vars: variables associated with stimvols
%      db: sorted variables and stimvols
%       o: parameters

function [vars,o,db] = getVarsFromStimvol(stimvols,variables,o)

%open session
cd(o.sessPath)
v = mrLoadRet([]);

%set group & scan
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',o.myScan);

%get all stimvols and variables
[d,o.variables] = getStimvol(v,variables,'taskNum',o.taskNum,'phaseNum',...
    o.phaseNum,'segmentNum',o.segmentNum);

nCond = length(o.variables);
j = zeros(1,nCond);
for ivar = 1 : length(variables)
    
    %init variable values
    db.(variables{ivar}) = [];
    
    %get variable
    a = regexp(o.variables,variables{ivar});
    
    %find variable values    
    thisVarpos = find(~cellfun('isempty',a));    
    
    %get variable values
    nValueThisVar = length(thisVarpos);
    for k = 1 : nValueThisVar
        %get trial #
        valPos = thisVarpos(k);
        volsThisVarValue = d{valPos};
        nTrialsThisVar = length(volsThisVarValue);
        
        %get this variable value
        thisVarValue = o.variables{valPos};
        posValue = regexp(thisVarValue,variables{ivar},'end')+2;
        valThisVar = str2double(thisVarValue(posValue:end));
        
        %add to previous values
        db.(variables{ivar}) = [db.(variables{ivar}); repmat(valThisVar,nTrialsThisVar,1)];
    end
    
    %get and sort stimvols matched with associated variables
    db.stimvols = cell2mat(d(thisVarpos));
    [db.stimvols,sorted] = sort(db.stimvols);
    db.(variables{ivar}) = db.(variables{ivar})(sorted);    
end

%retrieve input stimvols variables
o.stimvols = stimvols;
for i = 1 :  length(stimvols)
    
    %find stimvols in db
    [~,volsPos] = intersect(db.stimvols,stimvols{i});
    
    %get corresponding variables
    for ivar = 1 : length(variables)
        vars.(variables{ivar}){i} = db.(variables{ivar})(volsPos);
    end        
end
mrQuit






