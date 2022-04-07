% Client - Flywheel SDK Client instance
%
% Client Properties:
%    users - The users finder instance
%    groups - The groups finder instance
%    projects - The projects finder instance
%    subjects - The subjects finder instance
%    sessions - The sessions finder instance
%    acquisitions - The acquisitions finder instance
%    jobs - The jobs finder instance
%    gears - The gears finder instance
%    collections - The collections finder instance
%
% Client Methods:
%    addGroup - Create a new group
%    addGear - Create or update a gear.
%    getGear - Retrieve details about a specific gear
%    addCollection - Create a collection
%    addJob - Add a job
%    getConfig - Return public Scitran configuration information
%    getVersion - Get server and database schema version info
%    getCurrentUser - Get information about the current user
%    getModality - Get a modality's classification specification
%    resolve - Perform a path based lookup of nodes in the Flywheel hierarchy.
%    lookup - Perform a path based lookup of a single node in the Flywheel hierarchy.
%    fileUrl - Perform a path based lookup of a file in the Flywheel hierarchy, and return a single-use download URL.
%    downloadTar - Download the given set of containers as a tarball to dest_file.
%    View - Create a new View object
%    printViewColumns - Print the list of common View columns
%    readViewData - Execute a data view against container, and return the view data
%    saveViewData - Execute a data view against container, and return the view data
%    readViewStruct - Execute a data view against container, and return the view data as a struct array.
%    readViewTable - Execute a data view against container, and return the view data as a table.
classdef Client < handle
    properties(Hidden)
        fw_
    end
    properties(Dependent)
        users
        groups
        projects
        subjects
        sessions
        acquisitions
        jobs
        gears
        collections
    end
    methods
        function obj = Client(varargin)
            apiKey = [];
            % Check if the first argument is an api-key
            if numel(varargin) > 0 && ischar(varargin{1})
                apiKey = varargin{1};
                varargin = varargin(2:end);
            end
            if isempty(apiKey)
                apiKey = flywheel.Client.readApiKeyFromCLI();
            end
            if isempty(apiKey)
                throw(MException('Client:loginException', 'Must login with flywheel command-line interface, or specify an api key'));
            end
            obj.fw_ = flywheel.Flywheel(apiKey, varargin{:});
        end
        function [returnData, resp] = addUser(obj, varargin)
            % Add a user
            % body (User)
            % returns: [CommonObjectCreated, resp]
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [returnData, resp] = obj.fw_.addUser(body);
        end
        function [returnData, resp] = addGroup(obj, varargin)
            % Add a group
            % body (Group)
            % returns: [GroupNewOutput, resp]
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [returnData, resp] = obj.fw_.addGroup(body);
        end
        function [returnData, resp] = addGear(obj, gearName, body, varargin)
            % Create or update a gear.
            % gearName (char):Name of the gear to interact with
            % body (GearDoc)
            % returns: [CollectionNewOutput, resp]
            [returnData, resp] = obj.fw_.addGear(gearName, body, varargin{:});
        end
        function [returnData, resp] = getGear(obj, gearId, varargin)
            % Retrieve details about a specific gear
            % gearId (char):Id of the gear to interact with
            % returns: [GearDoc, resp]
            [returnData, resp] = obj.fw_.getGear(gearId, varargin{:});
        end
        function [returnData, resp] = addCollection(obj, varargin)
            % Create a collection
            % body (Collection)
            % returns: [CollectionNewOutput, resp]
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [returnData, resp] = obj.fw_.addCollection(body);
        end
        function [returnData, resp] = addJob(obj, body, varargin)
            % Add a job
            % body (Job)
            % returns: [CommonObjectCreated, resp]
            [returnData, resp] = obj.fw_.addJob(body, varargin{:});
        end
        function [returnData, resp] = getConfig(obj, varargin)
            % Return public Scitran configuration information
            % returns: [ConfigOutput, resp]
            [returnData, resp] = obj.fw_.getConfig(varargin{:});
        end
        function [returnData, resp] = getVersion(obj, varargin)
            % Get server and database schema version info
            % returns: [VersionOutput, resp]
            [returnData, resp] = obj.fw_.getVersion(varargin{:});
        end
        function [returnData, resp] = getCurrentUser(obj, varargin)
            % Get information about the current user
            % returns: [User, resp]
            [returnData, resp] = obj.fw_.getCurrentUser(varargin{:});
        end
        function [returnData, resp] = getModality(obj, modalityId, varargin)
            % Get a modality's classification specification
            % modalityId (str): The modality id
            % returns: [Modality, resp]
        end
        function [returnData, resp] = get(obj, id, varargin)
            % Retrieve the specified object by id.
            %
            % Objects that can be retrieved in this way are:
            %   group, project, session, subject, acquisition, analysis and collection
            % id (str): Id of the object to retrieve
            % returns: [ContainerOutput, resp]
            [returnData, resp] = obj.fw_.get(id, varargin{:});
        end
        function [returnData, resp] = resolve(obj, path, varargin)
            % Perform a path based lookup of nodes in the Flywheel hierarchy.
            % path (char): The path to resolve
            % returns: [ResolverOutput, resp]
            [returnData, resp] = obj.fw_.resolve(path, varargin{:});
        end
        function [returnData, resp] = lookup(obj, path, varargin)
            % Perform a path based lookup of a single node in the Flywheel hierarchy.
            % path (char): The path to resolve
            % returns: [ResolverOutput, resp]
            [returnData, resp] = obj.fw_.lookup(path, varargin{:});
        end
        function url = fileUrl(obj, path)
            % Perform a path based lookup of a file in the Flywheel hierarchy, and return a single-use download URL.
            % path (char): The path to resolve
            % returns: The file URL if found, otherwise raises an error
            url = obj.fw_.fileUrl(path);
        end
        function summary = downloadTar(obj, varargin)
            % Download the given set of containers as a tarball to dest_file.
            %
            % Supports downloading Projects, Sessions, Acquisitions and/or Analyses.
            %
            % Parameters:
            % containers: (required) The container, or list of containers to download.
            % destFile: (required) The destination file on disk
            % includeTypes: The optional list of types to include in the download (e.g. ['nifti'])
            % excludeTypes: The optional list of types to exclude from the download (e.g. ['dicom'])
            summary = obj.fw_.downloadTar(varargin{:});
        end
        function view = View(obj, varargin)
            % Create a new View object
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
            builder = flywheel.ViewBuilder(varargin{:});
            view = builder.build();
        end
        function printViewColumns(obj)
            % Print the list of common View columns
            obj.fw_.printViewColumns();
        end
        function data = readViewData(obj, view, containerId, varargin)
            % Execute a data view against container, and return the view data
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            data = obj.fw_.readViewData(view, containerId, varargin{:});
        end

        function destFile = saveViewData(obj, view, containerId, destFile, varargin)
            % Execute a data view against container, and return the view data
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            %   destFile: The destination file path
            destFile = obj.fw_.saveViewData(view, containerId, destFile, varargin{:});
        end

        function result = readViewStruct(obj, varargin)
            % Execute a data view against container, and return the view data as a struct array.
            % Requires Matlab 2016 or later
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            %   filter: The filter to apply
            %   skip: The number of rows to skip
            %   limit: The maximum number of rows to return
            result = obj.fw_.readViewStruct(varargin{:});
        end

        function result = readViewTable(obj, varargin)
            % Execute a data view against container, and return the view data as a table.
            % Requires Matlab 2016 or later
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            %   filter: The filter to apply
            %   skip: The number of rows to skip
            %   limit: The maximum number of rows to return
            result = obj.fw_.readViewTable(varargin{:});
        end
        function users = get.users(obj)
            % Returns the users finder object
            users = obj.fw_.users;
        end
        function groups = get.groups(obj)
            % Returns the groups finder object
            groups = obj.fw_.groups;
        end
        function projects = get.projects(obj)
            % Returns the projects finder object
            projects = obj.fw_.projects;
        end
        function subjects = get.subjects(obj)
            % Returns the subjects finder object
            subjects = obj.fw_.subjects;
        end
        function sessions = get.sessions(obj)
            % Returns the sessions finder object
            sessions = obj.fw_.sessions;
        end
        function acquisitions = get.acquisitions(obj)
            % Returns the acquisitions finder object
            acquisitions = obj.fw_.acquisitions;
        end
        function jobs = get.jobs(obj)
            % Returns the jobs finder object
            jobs = obj.fw_.jobs;
        end
        function gears = get.gears(obj)
            % Returns the gears finder object
            gears = obj.fw_.gears;
        end
        function collections = get.collections(obj)
            % Returns the collections finder object
            collections = obj.fw_.collections;
        end
        function varargout = subsref(obj, s)
            if strcmp(s(1).type, '.') && (isprop(obj, s(1).subs) || ismethod(obj, s(1).subs))
                [varargout{1:nargout}] = builtin('subsref', obj, s);
            else
                [varargout{1:nargout}] = subsref(obj.fw_, s);
            end
        end
    end

    methods(Static, Hidden)
        function apiKey = readApiKeyFromCLI()
            apiKey = [];

            try
                homedir = getenv('HOME');
                if isempty(homedir)
                    homedrive = getenv('HomeDrive');
                    homepath = getenv('HomePath');
                    if ~isempty(homedrive) && ~isempty(homepath)
                        homedir = [homedrive, homepath];
                    end
                end
                if isempty(homedir)
                    return
                end

                fid = fopen(sprintf('%s/.config/flywheel/user.json', homedir));
                raw = fread(fid, inf);
                fclose(fid);

                config = jsondecode(char(raw'));

                apiKey = config.key;
            catch ME
            end
        end
    end
end
