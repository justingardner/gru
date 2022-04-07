classdef (Abstract) Util
    methods(Static)
        function result = secondsToYears(seconds)
            result = seconds / 31557600.0;
        end
        function result = yearsToSeconds(years)
            result = int64(years * 31557600.0);
        end
        function result = secondsToMonths(seconds)
            result = seconds / 2592000.0;
        end
        function result = monthsToSeconds(months)
            result = int64(months * 2592000.0);
        end
        function result = secondsToWeeks(seconds)
            result = seconds / 604800.0;
        end
        function result = weeksToSeconds(weeks)
            result = int64(weeks * 604800.0);
        end
        function result = secondsToDays(seconds)
            result = seconds / 86400.0;
        end
        function result = daysToSeconds(days)
            result = int64(days * 86400.0);
        end

        function result = toRef(obj)
            if ismethod(obj, 'ref')
                result = ref(obj);
            else
                result = obj;
            end
        end

        function result = toDestination(obj)
            if ismethod(obj, 'ref')
                result = struct(ref(obj));
            else
                result = obj;
            end
        end

        function applyFnToStructOrCells(fn, value)
            if iscell(value)
                for i = 1:numel(value)
                    item = value{i};
                    if numel(item) ~= 2
                        throw(MException('ApiClient:apiException', 'Expected a key-value pair in cell array!'));
                    end
                    fn(item{1}, item{2});
                end
            else % struct
                fields = fieldnames(value);
                for i = 1:numel(fields)
                    name = fields{i};
                    fn(name, value.(name));
                end
            end
        end       
    end
end