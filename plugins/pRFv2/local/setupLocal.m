function localParams = setupLocal(~)

%% Setup local (initialize # of parallel processors)

localParams = struct;
localParams.bins = mlrNumWorkers(1);