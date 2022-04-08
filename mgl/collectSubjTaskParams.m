function collectSubjTaskParams(stimfilepath, root)

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
disp(sprintf('*** mlrReconAll RUNNING FOR SUBJECT %s ***',sid));
disp('*******************************************************');


%% get the files list
files = dir(fullfile(sprintf('%s/%s*.mat',stimfilepath, root)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'resp', [], 'nTrials', 0);
for fi = 1:length(files)
  load(fullfile(sprintf('%s/%s',stimfilepath, files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
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

disp(sprintf('SUBJECT %s: Found %i runs with a total of %i trials', sid, count, data.nTrials));
keyboard
