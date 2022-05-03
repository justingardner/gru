function data = collectSubjTaskParams(taskName, root)
%  collectSubjTaskParams(taskName, root)
%
%     Gets all stimfilles for task taskName and compiles parameters.
%     Usage: data = collectSubjTaskParams(taskName, root);
%
%       e.g. 
%               mglSetSID(350);
%               data = collectSubjTaskParams('mglRetinotopy');
%
if ieNotDefined('root')
    root = '';
end
%% Get the subject ID
sid = mglGetSID;
if isempty(sid)
    disp('*****');
    disp('(mlrReconAll) Please set the subject ID by calling mglSetSID(#)');
    disp('*****');
    return
end

disp('*******************************************************');
disp(sprintf('*** mlrReconAll RUNNING FOR SUBJECT %s TASK %s ***',sid, taskName));
disp('*******************************************************');


%% get the files list
files = dir(fullfile(sprintf('~/data/%s/%s/%s*.mat',taskName, sid, root)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'resp', [], 'nTrials', 0);
for fi = 1:length(files)
  load(fullfile(sprintf('~/data/%s/%s/%s',taskName, sid, files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
  if iscell(e); e = e{1}; end
  if e.nTrials>1
    
    f = fields(e.parameter);
    for i = 1:length(f)
        if ~isfield(data, f{i})
            data.(f{i}) = [];
        end
        data.(f{i}) = [data.(f{i}) e.parameter.(f{i})];
    end
    f = fields(e.randVars);
    for i = 1:length(f)
        if ~isfield(data, f{i})
            data.(f{i}) = [];
        end
        data.(f{i}) = [data.(f{i}) e.randVars.(f{i})];
    end
    
    data.resp = [data.resp e.response];
    data.reaction_time = [data.reaction_time e.reactionTime];
    data.nTrials = data.nTrials + e.nTrials;
        
  end
  count = count + 1;
end

disp(sprintf('SUBJECT %s: Found %i runs with a total of %i trials', sid, count, data.nTrials))
