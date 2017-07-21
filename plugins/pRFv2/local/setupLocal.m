function localParams = setupLocal(~)

%% Setup local (initialize # of parallel processors)

localParams = struct;
localParams.bins = mlrNumWorkers(1);

if localParams.bins<0
    localParams.bins = 0;
end