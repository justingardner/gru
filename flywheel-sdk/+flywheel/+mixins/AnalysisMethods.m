classdef AnalysisMethods < handle
    methods
        function [returnData, resp] = addAnalysis(obj, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [analysisId, resp] = obj.invokeContainerApi('add%sAnalysis', obj.get('id'), body);
            if ~isempty(analysisId)
                returnData = obj.context_.getAnalysis(analysisId);
            else
                returnData = [];
            end
        end
    end
end
