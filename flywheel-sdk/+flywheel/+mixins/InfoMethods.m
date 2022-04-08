classdef InfoMethods < handle
    methods
        function [returnData, resp] = replaceInfo(obj, info)
            [returnData, resp] = obj.invokeContainerApi('replace%sInfo', obj.get('id'), info);
        end
        function [returnData, resp] = updateInfo(obj, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [returnData, resp] = obj.invokeContainerApi('set%sInfo', obj.get('id'), body);
        end
        function [returnData, resp] = deleteInfo(obj, keys)
            [returnData, resp] = obj.invokeContainerApi('delete%sInfoFields', obj.get('id'), keys);
        end
    end
end
