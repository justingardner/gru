classdef TimestampMethods < handle
    properties(Dependent)
        localTimestamp
        originalTimestamp
    end
    methods
        function result = get.localTimestamp(obj)
            result = obj.localizeDate('timestamp');
        end
        function result = get.originalTimestamp(obj)
            result = obj.get('timestamp');
            if ~isempty(result)
                timezone = obj.get('timezone');
                if isempty(timezone)
                    throw(MException('ApiClient:apiException', 'No original timezone was specified!'));
                end
                result = datetime(result, 'TimeZone', timezone);
            end
        end
    end
end
