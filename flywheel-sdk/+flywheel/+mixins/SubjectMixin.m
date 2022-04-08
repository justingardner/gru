classdef SubjectMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.InfoMethods ...
        & flywheel.mixins.NoteMethods ...
        & flywheel.mixins.FileMethods ...
        & flywheel.mixins.AnalysisMethods
    properties
        containerType_ = 'Subject';
        sessions_ = false;
        analyses_ = false;
    end
    properties(Dependent)
        sessions
    end
    methods
        function sessions = get.sessions(obj)
            sessions = obj.getChildren('Sessions');
        end
        function analyses = get_analyses(obj)
            analyses = obj.getChildren('Analyses');
        end
        function [returnData, resp] = addSession(obj, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            body.project = obj.project;
            if ~isfield(body, 'subject')
                body.subject = struct();
            end
            body.subject.id = obj.id;
            [sessionId, resp] = addSession(obj.context_, body);
            if ~isempty(sessionId)
                returnData = getSession(obj.context_, sessionId);
            else
                returnData = [];
            end
        end
    end
end