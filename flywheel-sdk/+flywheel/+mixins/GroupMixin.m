classdef GroupMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.PermissionMethods
    properties
        containerType_ = 'Group';
        projects_ = false;
    end
    properties(Dependent)
        projects
    end
    methods
        function projects = get.projects(obj)
            projects = obj.getChildren('Projects');
        end
        function [returnData, resp] = addProject(obj, varargin)
            [returnData, resp] = obj.addChild('Project', varargin{:});
        end
    end
end