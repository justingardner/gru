classdef JobMixin < flywheel.mixins.ContainerBase
    properties
        containerType_ = 'Job';
    end
    methods
        function result = changeState(obj, newState)
            result = obj.context_.changeJobState(obj.id, newState);
        end
        function result = getLogs(obj)
            response = obj.context_.getJobLogs(obj.id);
            result = response.logs;
        end
        function obj = printLogs(obj, varargin)
            p = inputParser;
            addParameter(p, 'file', 1);
            
            parse(p, varargin{:});

            logs = obj.getLogs();
            for i = 1:numel(logs)
                fprintf(p.Results.file, '%s', logs{i}.msg);
            end
            fprintf(p.Results.file, '');
        end
    end
end
