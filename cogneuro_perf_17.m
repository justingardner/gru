%% Check cogneuro performance

folder = '~/data/cogneuro_workingmemory_corey/s363';

files = dir(fullfile(folder,'*.mat'));

correct = []; taskv = [];
for fi = 5:length(files)
    load(fullfile(folder,files(fi).name));
    
    e = getTaskParameters(myscreen,task);
    e = e{1};
    
    taskv = [taskv e.parameter.task];
    correct = [correct e.randVars.correct];
end

c1 = nanmean(correct(taskv==1));
c2 = nanmean(correct(taskv==2));

%% Check cogneuro performance

folder = '~/data/cogneuro_workingmemory_manasi/s364';

files = dir(fullfile(folder,'*.mat'));

correct = []; attend = [];
for fi = 5:length(files)
    load(fullfile(folder,files(fi).name));
    
    e = getTaskParameters(myscreen,task);
    e = e{1};
    
    attend = [attend e.randVars.length];
    correct = [correct e.randVars.correct];
end

c1 = nanmean(correct(attend==8));
c2 = nanmean(correct(attend==16));
c3 = nanmean(correct(attend==0));

%% Check cogneuro performance

folder = '~/data/cogneuro_attention_isabel/s361';

files = dir(fullfile(folder,'*.mat'));

correct = []; attend = [];
for fi = 5:length(files)
    load(fullfile(folder,files(fi).name));
    
    e = getTaskParameters(myscreen,task);
    e = e{1};
    
    attend = [attend e.parameter.attend];
    correct = [correct e.randVars.correct];
end

c1 = nanmean(correct(attend==1));
c2 = nanmean(correct(attend==2));
c3 = nanmean(correct(attend==0));
%% Check cogneuro performance

folder = '~/data/cogneuro_attention_eshed/s362';

files = dir(fullfile(folder,'*.mat'));

correct = []; attend = [];
for fi = 5:length(files)
    load(fullfile(folder,files(fi).name));
    
    e = getTaskParameters(myscreen,task);
    e = e{1};
    
    attend = [attend e.parameter.attend];
    correct = [correct e.randVars.correct];
end

c1 = nanmean(correct(attend==1));
c2 = nanmean(correct(attend==2));
c3 = nanmean(correct(attend==0));