classdef NoteMethods < handle
    methods
        function [returnData, resp] = addNote(obj, message)
            [returnData, resp] = obj.invokeContainerApi('add%sNote', obj.get('id'), message);
        end
        function [returnData, resp] = deleteNote(obj, noteId)
            [returnData, resp] = obj.invokeContainerApi('delete%sNote', obj.get('id'), noteId);
        end
    end
end
