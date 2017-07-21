function prf = pRFController( splits, params )
%PRFCONTROLLER 
%   Launch the controller function and then call the merge functions when
%   necessary.
%
%   Usage
%       Called by pRF
%
%   Author
%       Dan Birman (dbirman@stanford.edu)
%       7/18/17

global controller

%% Setup global variables

controller.prfName = params.saveName;
controller.curPath = pwd;
controller.sherlockSessionPath = ['/share/PI/jlg/' controller.curPath(findstr(controller.curPath, 'data'):end)];
controller.suid = getsuid;
controller.scriptsDir = 'Splits/Scripts';

%% Initial setup and function calls
prf = struct;

prf.sherlock.setup = @setupSherlock;
prf.sherlock.add = @addSherlockJob;
prf.sherlock.check = @checkSherlockJob;
prf.local.setup = @setupLocal;
prf.local.add = @addLocalJob;
prf.local.check = @checkLocalJob;

groups = fields(prf);

%% Startup
pdisp(20);
pdisp('  Spinning up controller...');
pdisp(20);

%% Call setup functions
for gi = 1:length(groups)
    pdisp(sprintf('  Calling %s setup',groups{gi}));
    prf.(groups{gi}).params = prf.(groups{gi}).setup(params);
    prf.(groups{gi}).bins = prf.(groups{gi}).params.bins;
    prf.(groups{gi}).running = cell(1,prf.(groups{gi}).bins);
    if prf.(groups{gi}).bins == 0
        % actually there's no bins here, delete this
        disp(sprintf('Group %i has no bins available',groups{gi}));
        prf = rmfield(prf,groups{gi});
    end
end
groups = fields(prf);
pdisp(20);

% get intial job # (these get added quickly at the start)
sentjobs = 0;
initjobs = 0;
for gi = 1:length(groups)
    initjobs = initjobs + prf.(groups{gi}).bins;
end

%% Start running the splits
nSplits = length(splits);
queue = splits;
finished = {};

%% Send info to the user
pdisp(sprintf('Found %i splits',nSplits));
for gi = 1:length(groups)
    pdisp(sprintf('There are %i bins available on %s',prf.(groups{gi}).bins,groups{gi}));
end
pdisp(20);

%% If no analysis directory, generate it
if ~isdir(fullfile(controller.curPath,'Splits','Analysis'))
    mkdir(fullfile(controller.curPath,'Splits','Analysis'));
end
%% Cycle
tic
ntoc = toc;
while length(finished)<length(splits)
    % ADD
    torun = ~cellfun(@isempty,queue);
    if any(torun)
        next = find(torun,1);
        % we have available jobs
        csplit = queue{next};
        
        % check the bins
        for gi = 1:length(groups)
            emptyslots = cellfun(@isempty,prf.(groups{gi}).running);
            if any(emptyslots)
                slot = find(emptyslots,1);
                pdisp(sprintf('Split %i will be added to %s in slot #%i',csplit.num,groups{gi},slot));
%                 try 
                    prf.(groups{gi}).running{slot} = prf.(groups{gi}).add(csplit);
                    prf.(groups{gi}).running{slot}.timestamp.start = toc;
                    % drop from queue
                    queue{next} = [];
                    sentjobs = sentjobs+1;
                    break;
%                 catch
%                     pdisp(sprintf('Warning: a job on %s failed to be added. Ignoring and attempting to continue automatically.',groups{gi}));
%                 end      
            end 
        end        
    end
    % CHECK (delay until # initjobs have all been sent)
    if sentjobs>=initjobs
        for gi = 1:length(groups)
            for ri = 1:length(prf.(groups{gi}).running)
                csplit = prf.(groups{gi}).running{ri};
                if ~isempty(csplit)
                    if prf.(groups{gi}).check(csplit)
                        pdisp(sprintf('Split %i completed successfully, recycling slot %i on %s',csplit.num,ri,groups{gi}));
                        csplit.timestamp.end = toc;
                        csplit.timestamp.elapsed = csplit.timestamp.end - csplit.timestamp.start;
                        finished{end+1} = csplit;
                        prf.(groups{gi}).running{ri} = [];
                    end
                end
            end
        end
        if (toc-ntoc)>60
            pdisp(sprintf('Elapsed time %s',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')));
            ntoc = 60*floor(toc/60);
        end
        pause(5);
    end
end

%% Cleanup
prf.splits = finished;

%% End
pdisp(20);
pdisp('  Shutting down controller...');
pdisp(20);

function pdisp(str)
if isnumeric(str)
    n = str;
    str = '';
    for i = 1:n
        str = strcat(str,'*');
    end
end
disp(sprintf('(pRFControl) %s',str));

