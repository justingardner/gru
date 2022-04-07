classdef ProjectMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.InfoMethods ...
        & flywheel.mixins.NoteMethods ...
        & flywheel.mixins.FileMethods ...
        & flywheel.mixins.PermissionMethods ...
        & flywheel.mixins.DownloadMethods ...
        & flywheel.mixins.AnalysisMethods
    properties
        containerType_ = 'Project';
        subjects_ = false;
        sessions_ = false;
        analyses_ = false;
    end
    properties(Dependent)
        subjects
        sessions
    end
    methods
        function subjects = get.subjects(obj)
            subjects = obj.getChildren('Subjects');
        end
        function sessions = get.sessions(obj)
            sessions = obj.getChildren('Sessions');
        end
        function analyses = get_analyses(obj)
            analyses = obj.getChildren('Analyses');
        end
        function [returnData, resp] = addSubject(obj, varargin)
            [returnData, resp] = obj.addChild('Subject', varargin{:});
        end
        function [returnData, resp] = addSession(obj, varargin)
            [returnData, resp] = obj.addChild('Session', varargin{:});
        end
    end
end