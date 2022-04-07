classdef PermissionMethods < handle
    methods
        function [returnData, resp] = addPermission(obj, permission)
            [returnData, resp] = obj.invokeContainerApi('add%sPermission', obj.get('id'), permission);
        end
        function [returnData, resp] = updatePermission(obj, userId, permission)
            [returnData, resp] = obj.invokeContainerApi('modify%sUserPermission', obj.get('id'), userId, permission);
        end
        function [returnData, resp] = deletePermission(obj, userId)
            [returnData, resp] = obj.invokeContainerApi('delete%sUserPermission', obj.get('id'), userId);
        end
    end
end
