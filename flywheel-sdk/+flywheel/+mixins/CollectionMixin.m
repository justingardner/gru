classdef CollectionMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.InfoMethods ...
        & flywheel.mixins.NoteMethods ...
        & flywheel.mixins.FileMethods ...
        & flywheel.mixins.PermissionMethods
    properties
        containerType_ = 'Collection';
        sessions_ = false;
        acquisitions_ = false;
        analyses_ = false;
    end
    properties(Dependent)
        sessions
        acquisitions
    end
    methods
        function sessions = get.sessions(obj)
            sessions = obj.getChildren('Sessions');
        end
        function acquisitions = get.acquisitions(obj)
            acquisitions = obj.getChildren('Acquisitions');
        end
        function [returnData, resp] = addSessions(obj, varargin)
            ids = {};
            for i = 1:numel(varargin)
                arg = varargin{i};
                if ischar(arg)
                    ids = [ids; arg];
                else
                    ids = [ids; arg.id];
                end
            end
            [returnData, resp] = obj.invokeContainerApi('addSessionsTo%s', obj.get('id'), ids);
        end
        function [returnData, resp] = addAcquisitions(obj, varargin)
            ids = {};
            for i = 1:numel(varargin)
                arg = varargin{i};
                if ischar(arg)
                    ids = [ids; arg];
                else
                    ids = [ids; arg.id];
                end
            end
            [returnData, resp] = obj.invokeContainerApi('addAcquisitionsTo%s', obj.get('id'), ids);
        end
        function analyses = get_analyses(obj)
            analyses = obj.getChildren('Analyses');
        end
    end
end