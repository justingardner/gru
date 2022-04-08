classdef ViewBuilder < handle
    % ViewBuilder A builder that assists in constructing a View object.
    properties (Access=private)
        label_ = [];
        public_ = false;
        columns_ = [];
        fileColumns_ = [];
        fileContainer_ = [];
        fileFilter_ = [];
        fileZipFilter_ = [];
        fileFormat_ = [];
        fileFormatOpts_ = [];
        fileMatch_ = [];
        processFiles_ = [];
        analysisFilter_ = [];
        includeIds_ = true;
        includeLabels_ = true;
        missingDataStrategy_ = [];
    end
    methods
        function obj = ViewBuilder(varargin)
            % Create a new ViewBuilder object
            %
            % Parameters:
            %   public: Whether or not to make this data view public when saving it.
            %   match: The file match type, one of: first, last, newest, oldest, all
            %   zipFiles: The zip file filter, see the zip_member_filter function
            %   columns: The columns or column groups to add
            %   processFiles: Whether or not to process files, default is true
            %   includeIds: Whether or not to include id columns, default is true
            %   includeLabels: Whether or not to include label columns, default is true
            %   container: When matching files, the container to match on
            %   filename: When matching files, the filename pattern to match
            %   analysisLabel: When matching analysis files, the label match string
            %   analysisGearName: When matching analysis files, the gear name match string
            %   analysisGearVersion: When matching analysis files, the gear version match string
            p = inputParser;
            addParameter(p, 'label', []);
            addParameter(p, 'public', false);
            addParameter(p, 'match', []);
            addParameter(p, 'zipFiles', []);
            addParameter(p, 'columns', []);
            addParameter(p, 'processFiles', []);
            addParameter(p, 'includeIds', true);
            addParameter(p, 'includeLabels', true);
            addParameter(p, 'container', []);
            addParameter(p, 'filename', []);
            addParameter(p, 'analysisLabel', []);
            addParameter(p, 'analysisGearName', []);
            addParameter(p, 'analysisGearVersion', []);

            parse(p, varargin{:});
            obj.label_ = p.Results.label;
            obj.public_ = p.Results.public;
            obj.fileMatch_ = p.Results.match;
            obj.processFiles_ = p.Results.processFiles;
            obj.includeIds_ = p.Results.includeIds;
            obj.includeLabels_ = p.Results.includeLabels;

            if ~isempty(p.Results.zipFiles)
                obj.zipMemberFilter(p.Results.zipFiles);
            end

            if ~isempty(p.Results.filename)
                obj.files(p.Results.container, p.Results.filename, ...
                    'analysisLabel', p.Results.analysisLabel, ...
                    'analysisGearName', p.Results.analysisGearName, ...
                    'analysisGearVersion', p.Results.analysisGearVersion);
            end

            % Add column/columns
            if ~isempty(p.Results.columns)
                if iscell(p.Results.columns)
                    for i = 1:numel(p.Results.columns)
                        col = p.Results.columns{i};
                        if iscell(col)
                            obj.column(col{:});
                        else
                            obj.column(col);
                        end
                    end
                else
                    obj.column(p.Results.columns);
                end
            end
        end
        function obj = label(obj, varargin)
            % Set the label for this view
            p = inputParser;
            addRequired(p, 'label');
            parse(p, varargin{:});
            obj.label_ = p.Results.label;
        end
        function obj = public(obj, varargin)
            % Set whether or not this data view should be made public.
            p = inputParser;
            addOptional(p, 'public', true);
            parse(p, varargin{:});
            obj.public_ = p.Results.public;
        end
        function obj = column(obj, varargin)
            % Define a column for this data view.
            %
            % Parameters:
            %   src: The source field, or column alias name.
            %   dst: The optional destination field (defaults to source)
            %   type: The optional type for this column, one of: int, float, string bool.
            p = inputParser;
            addRequired(p, 'src');
            addOptional(p, 'dst', [], @ischar);
            addOptional(p, 'type', []);
            parse(p, varargin{:});

            col = p.Results;
            if isstruct(col.src)
                col = col.src;
            end

            col = obj.preprocessColumn_(col);
            obj.columns_{end+1} = flywheel.model.DataViewColumnSpec(col);
        end
        function obj = fileColumn(obj, varargin)
            % Define a column to extract from a file.
            %
            % Parameters:
            %   src: The source field, or column alias name.
            %   dst: The optional destination field (defaults to source)
            %   type: The optional type for this column, one of: int, float, string bool.
            p = inputParser;
            addRequired(p, 'src');
            addOptional(p, 'dst', [], @ischar);
            addOptional(p, 'type', []);
            parse(p, varargin{:});

            if isstruct(p.Results.src)
                obj.fileColumns_{end+1} = flywheel.model.DataViewColumnSpec(p.Results.src);
            else
                obj.fileColumns_{end+1} = flywheel.model.DataViewColumnSpec(p.Results);
            end
        end
        function obj = files(obj, varargin)
            % Set filter for matching files
            %
            % Container is one of project, subject, session, acquisition
            % Filename filters can use the (\*, ?) wildcards
            % Analysis filters also support wildcards
            %
            % Parameters:
            %   container: When matching files, the container to match on: one of project, subject, session, acquisition
            %   filename: When matching files, the filename pattern to match
            %   analysisLabel: When matching analysis files, the label match string
            %   analysisGearName: When matching analysis files, the gear name match string
            %   analysisGearVersion: When matching analysis files, the gear version match string
            p = inputParser;
            addRequired(p, 'container');
            addRequired(p, 'filename');
            addParameter(p, 'analysisLabel', []);
            addParameter(p, 'analysisGearName', []);
            addParameter(p, 'analysisGearVersion', []);
            parse(p, varargin{:});

            obj.fileContainer_ = p.Results.container;
            obj.fileFilter_ = flywheel.model.DataViewNameFilterSpec('value', p.Results.filename);
            if ~isempty(p.Results.analysisLabel) || ~isempty(p.Results.analysisGearName) || ~isempty(p.Results.analysisGearVersion)
                obj.analysisFilter('label', p.Results.analysisLabel, ...
                    'gearName', p.Results.analysisGearName, ...
                    'gearVersion', p.Results. analysisGearVersion);
            end
        end
        function obj = fileContainer(obj, varargin)
            % Set the container where files should be matched.
            p = inputParser;
            addRequired(p, 'container');
            parse(p, varargin{:});
            obj.fileContainer_ = p.Results.container;
        end
        function obj = fileMatch(obj, varargin)
            % Set the resolution strategy if multiple matching files or analyses are encountered.
            %
            % The file match type is one of: first, last, newest, oldest, all
            p = inputParser;
            addRequired(p, 'match');
            parse(p, varargin{:});
            obj.fileMatch_ = p.Results.match;
        end
        function obj = analysisFilter(obj, varargin)
            % Set the filter to use for matching analyses. If this is set, then analyses files will be matched instead of container.
            %
            % Parameters:
            %   label: The label match string, wildcards (\*, ?) are supported.
            %   gearName: The gear name match string, wildcards (\*, ?) are supported.
            %   gearVersion: The gear version match string, wildcards (\*, ?) are supported.
            %   regex: Whether to treat the match string as a regular expression (default is False)
            p = inputParser;
            addParameter(p, 'label', []);
            addParameter(p, 'gearName', []);
            addParameter(p, 'gearVersion', []);
            addParameter(p, 'regex', false);
            parse(p, varargin{:});

            if ~isempty(obj.analysisFilter_)
                obj.analysisFilter_ = flywheel.model.DataViewAnalysisFilterSpec()
            end

            if ~isempty(p.Results.label)
                obj.analysisFilter_.label = flywheel.model.DataViewNameFilterSpec('value', p.Results.label, 'regex', p.Results.regex);
            end
            if ~isempty(p.Results.gearName)
                obj.analysisFilter_.gearName = flywheel.model.DataViewNameFilterSpec('value', p.Results.gearName, 'regex', p.Results.regex);
            end
            if ~isempty(p.Results.gearVersion)
                obj.analysisFilter_.gearVersion = flywheel.model.DataViewNameFilterSpec('value', p.Results.gearVersion, 'regex', p.Results.regex);
            end
        end
        function obj = fileFilter(obj, varargin)
            % Set the filter to use for matching files.
            %
            % Parameters:
            %   value: The filename match string, wildcards (\*, ?) are supported.
            %   regex: Whether to treat the match string as a regular expression (default is False)
            p = inputParser;
            addRequired(p, 'value');
            addParameter(p, 'regex', false);
            parse(p, varargin{:});
            obj.fileFilter_ = flywheel.model.DataViewNameFilterSpec(p.Results);
        end
        function obj = zipMemberFilter(obj, varargin)
            % Set the filter to use for matching members of a zip file.
            %
            % Parameters:
            %   value: The filename match string, wildcards (\*, ?) are supported.
            %   regex: Whether to treat the match string as a regular expression (default is False)
            p = inputParser;
            addRequired(p, 'value');
            addParameter(p, 'regex', false);
            parse(p, varargin{:});
            obj.fileZipFilter_ = flywheel.model.DataViewNameFilterSpec(p.Results);
        end
        function obj = fileFormat(obj, varargin)
            % Set the expected format of files to read.
            %
            % The format is one of: csv, tsv, json.
            % NOTE: This shouldn't be needed very often. If not specified, autodetect will be used for processing files.
            p = inputParser;
            addRequired(p, 'format');
            parse(p, varargin{:});
            obj.fileFormat_ = p.Results.format;
        end
        function obj = fileFormatOptions(obj, varargin)
            % Set additional options for the file format. (e.g. arguments to be passed to csv reader function)
            p = inputParser;
            p.KeepUnmatched = 1;
            parse(p, varargin{:});
            if isempty(obj.fileFormatOpts_)
                obj.fileFormatOpts_ = struct;
            end
            for name = fieldnames(p.Unmatched)'
                obj.fileFormatOpts_.(name{1}) = p.Unmatched.(name{1});
            end
        end
        function obj = processFiles(obj, varargin)
            % Set whether or not to process files (default is True)
            %
            % By default, files will be read and return a row for each row in the file. If you just want file attributes or info
            % instead, you can set this to False.
            p = inputParser;
            addRequired(p, 'value');
            parse(p, varargin{:});
            obj.processFiles_ = p.Results.value;
        end
        function obj = includeLabels(obj, varargin)
            % Set whether or not to include the label columns by default.
            p = inputParser;
            addRequired(p, 'value');
            parse(p, varargin{:});
            obj.includeLabels_ = p.Results.value;
        end
        function obj = includeIds(obj, varargin)
            % Set whether or not to include the id columns by default.
            p = inputParser;
            addRequired(p, 'value');
            parse(p, varargin{:});
            obj.includeIds_ = p.Results.value;
        end
        function obj = missingDataStrategy(obj, varargin)
            % Set the resolution strategy if rows are missing data for a column. The default is to replace the column value with None.
            % The strategy is one of:  none, drop-row
            p = inputParser;
            addRequired(p, 'value');
            parse(p, varargin{:});
            obj.missingDataStrategy_ = p.Results.value;
        end
        function view = build(obj)
            % Build the DataView constructed with this builder.
            fileSpec = [];
            if ~isempty(obj.fileContainer_) && ~isempty(obj.fileFilter_)
                fileSpec = flywheel.model.DataViewFileSpec('container', obj.fileContainer_, ...
                    'analysisFilter', obj.analysisFilter_, ...
                    'filter', obj.fileFilter_, ...
                    'zipMember', obj.fileZipFilter_, ...
                    'match', obj.fileMatch_, ...
                    'format', obj.fileFormat_, ...
                    'formatOptions', obj.fileFormatOpts_, ...
                    'processFiles', obj.processFiles_, ...
                    'columns', obj.fileColumns_);
            end

            view = flywheel.model.DataView('label', obj.label_, ...
                'public', obj.public_, ...
                'columns', obj.columns_, ...
                'fileSpec', fileSpec, ...
                'includeIds', obj.includeIds_, ...
                'includeLabels', obj.includeLabels_, ...
                'missingDataStrategy', obj.missingDataStrategy_);
        end
    end
    methods (Access=private)
        function result = preprocessColumn_(obj, col)
            srcParts = strsplit(col.src, '.');
            fileIdx = find(strcmp(srcParts, 'file'));
            if isempty(fileIdx) || fileIdx < 2
                result = col;
                return
            end

            analysisContainer = fileIdx > 2;
            fileContainer = srcParts{1};

            % If we currently have a file filter, validate
            if ~isempty(obj.fileContainer_)
                if ~strcmp(fileContainer, obj.fileContainer_)
                    error(sprintf('Can only select files on one container (%s already selected)', obj.fileContainer_));
                end
                if analysisContainer && isempty(obj.analysisFilter_)
                    error(sprintf('Can only select files on one container (%s already selected)', obj.fileContainer_));
                elseif ~analysisContainer && ~isempty(obj.analysisFilter_)
                    error(sprintf('Can only select files on one container (%s analyses already selected)', obj.fileContainer_));
                end
            else
                % Setup file matches, matching all files and dropping rows with missing data
                obj.fileContainer_ = fileContainer;
                obj.fileFilter_ = flywheel.model.DataViewNameFilterSpec('value', '*');
                if analysisContainer
                    labelFilter = flywheel.model.DataViewNameFilterSpec('value', '*');
                    obj.analysisFilter_ = flywheel.model.DataViewAnalysisFilterSpec('label', labelFilter);
                end
                obj.fileMatch_ = 'all';
                obj.missingDataStrategy_ = 'drop-row';
                obj.processFiles_ = false;
            end

            if isempty(col.dst)
                col.dst = col.src;
            end
            col.src = strjoin(srcParts(1, fileIdx:end), '.');
            result = col;
        end
    end
end