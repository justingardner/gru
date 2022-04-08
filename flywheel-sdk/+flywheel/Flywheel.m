% Flywheel - Global namespace for API calls
%
% Flywheel Properties:
%    apiClient - The api client instance
%    acquisitionsApi - Acquisition operations
%    analysesApi - Analysis operations
%    batchApi - Batch job operations
%    bulkApi - Bulk container operations
%    callbacksApi - 
%    collectionsApi - Collection operations
%    containersApi - Abstract container operations
%    dataexplorerApi - Search operations
%    defaultApi - 
%    devicesApi - Device operations
%    dimseApi - DIMSE configuration
%    filesApi - File upload/download operations
%    gearsApi - Gear operations
%    groupsApi - Group operations
%    jobsApi - Job operations
%    modalitiesApi - Modality operations
%    projectsApi - Project operations
%    reportsApi - Site-wide reports
%    rulesApi - Gear rule configuration
%    sessionsApi - Session operations
%    siteApi - Site level operations
%    subjectsApi - 
%    usersApi - User operations
%    viewsApi - Data view operations
%
% Flywheel Methods:
%    addAcquisition - Create a new acquisition
%    addAcquisitionAnalysis - Create an analysis and upload files.
%    addAcquisitionAnalysisNote - Add a note to acquisition analysis.
%    addAcquisitionNote - Add a note to acquisition.
%    addAcquisitionTag - Add a tag to acquisition.
%    deleteAcquisition - Delete a acquisition
%    deleteAcquisitionAnalysis - Delete an anaylsis
%    deleteAcquisitionAnalysisNote - Remove a note from acquisition analysis.
%    deleteAcquisitionFile - Delete a file
%    deleteAcquisitionNote - Remove a note from acquisition
%    deleteAcquisitionTag - Delete a tag
%    downloadAcquisitionAnalysisInputs - Download analysis inputs.
%    downloadAcquisitionAnalysisOutputs - Download analysis outputs.
%    downloadFileFromAcquisition - Download a file.
%    getAcquisitionFileZipInfo - Download a file.
%    getAcquisitionDownloadUrl - Download a file.
%    downloadInputFromAcquisitionAnalysis - Download analysis inputs with filter.
%    getAcquisitionAnalysisInputZipInfo - Download analysis inputs with filter.
%    getAcquisitionAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromAcquisitionAnalysis - Download analysis outputs with filter.
%    getAcquisitionAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getAcquisitionAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    getAcquisition - Get a single acquisition
%    getAcquisitionAnalyses - Get analyses for acquisition.
%    getAcquisitionAnalysis - Get an analysis.
%    getAcquisitionFileInfo - Get info for a particular file.
%    getAcquisitionNote - Get a note on acquisition.
%    getAcquisitionTag - Get the value of a tag, by name.
%    getAllAcquisitions - Get a list of acquisitions
%    modifyAcquisition - Update a acquisition
%    modifyAcquisitionAnalysis - Modify an analysis.
%    modifyAcquisitionFile - Modify a file's attributes
%    modifyAcquisitionFileClassification - Update classification for a particular file.
%    modifyAcquisitionFileInfo - Update info for a particular file.
%    modifyAcquisitionInfo - Update or replace info for a acquisition.
%    modifyAcquisitionNote - Update a note on acquisition.
%    renameAcquisitionTag - Rename a tag.
%    replaceAcquisitionFile - Replace a file
%    uploadFileToAcquisition - Upload a file to acquisition.
%    uploadOutputToAcquisitionAnalysis - Upload an output file to analysis.
%    addAnalysisNote - Add a note to analysis.
%    addAnalysisTag - Add a tag to analysis.
%    deleteAnalysisNote - Remove a note from analysis
%    deleteAnalysisTag - Delete a tag
%    downloadAnalysisInputs - Download analysis inputs.
%    downloadAnalysisOutputs - Download analysis outputs.
%    downloadInputFromAnalysis - Download analysis inputs with filter.
%    getAnalysisInputZipInfo - Download analysis inputs with filter.
%    getAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromAnalysis - Download analysis outputs with filter.
%    getAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    getAnalyses - Get nested analyses for a container
%    getAnalysis - Get an analysis.
%    getAnalysisNote - Get a note on analysis.
%    getAnalysisTag - Get the value of a tag, by name.
%    modifyAnalysis - Modify an analysis.
%    modifyAnalysisInfo - Update or replace info for a analysis.
%    modifyAnalysisNote - Update a note on analysis.
%    renameAnalysisTag - Rename a tag.
%    uploadOutputToAnalysis - Upload an output file to analysis.
%    cancelBatch - Cancel a Job
%    createBatchJobFromJobs - Create a batch job proposal from preconstructed jobs and insert it as 'pending'.
%    getAllBatches - Get a list of batch jobs the user has created.
%    getBatch - Get batch job details.
%    proposeBatch - Create a batch job proposal and insert it as 'pending'.
%    startBatch - Launch a job.
%    bulkCopy - Perform a bulk copy operation
%    bulkDelete - Perform a bulk delete operation
%    bulkMove - Perform a bulk move operation
%    bulkMoveSessions - Perform a bulk move of sessions to either a subject or project
%    callbackVirusScan - Callback url to send the virus scan result of a file.
%    addCollection - Create a collection
%    addCollectionAnalysis - Create an analysis and upload files.
%    addCollectionAnalysisNote - Add a note to collection analysis.
%    addCollectionNote - Add a note to collection.
%    addCollectionPermission - Add a permission
%    addCollectionTag - Add a tag to collection.
%    deleteCollection - Delete a collection
%    deleteCollectionAnalysis - Delete an anaylsis
%    deleteCollectionAnalysisNote - Remove a note from collection analysis.
%    deleteCollectionFile - Delete a file
%    deleteCollectionNote - Remove a note from collection
%    deleteCollectionTag - Delete a tag
%    deleteCollectionUserPermission - Delete a permission
%    downloadCollectionAnalysisInputs - Download analysis inputs.
%    downloadCollectionAnalysisOutputs - Download analysis outputs.
%    downloadFileFromCollection - Download a file.
%    getCollectionFileZipInfo - Download a file.
%    getCollectionDownloadUrl - Download a file.
%    downloadInputFromCollectionAnalysis - Download analysis inputs with filter.
%    getCollectionAnalysisInputZipInfo - Download analysis inputs with filter.
%    getCollectionAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromCollectionAnalysis - Download analysis outputs with filter.
%    getCollectionAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getCollectionAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    getAllCollections - List all collections.
%    getAllCollectionsCurators - List all curators of collections
%    getCollection - Retrieve a single collection
%    getCollectionAcquisitions - List acquisitions in a collection
%    getCollectionAnalyses - Get analyses for collection.
%    getCollectionAnalysis - Get an analysis.
%    getCollectionFileInfo - Get info for a particular file.
%    getCollectionNote - Get a note on collection.
%    getCollectionSessions - List sessions in a collection
%    getCollectionTag - Get the value of a tag, by name.
%    getCollectionUserPermission - List a user's permissions for this collection.
%    modifyCollection - Update a collection and its contents
%    modifyCollectionAnalysis - Modify an analysis.
%    modifyCollectionFile - Modify a file's attributes
%    modifyCollectionFileClassification - Update classification for a particular file.
%    modifyCollectionFileInfo - Update info for a particular file.
%    modifyCollectionInfo - Update or replace info for a collection.
%    modifyCollectionNote - Update a note on collection.
%    modifyCollectionUserPermission - Update a user's permission for this collection.
%    renameCollectionTag - Rename a tag.
%    replaceCollectionFile - Replace a file
%    uploadFileToCollection - Upload a file to collection.
%    uploadOutputToCollectionAnalysis - Upload an output file to analysis.
%    addContainerAnalysis - Create an analysis and upload files.
%    addContainerAnalysisNote - Add a note to container analysis.
%    addContainerNote - Add a note to container.
%    addContainerTag - Add a tag to container.
%    checkUidsExist - Check for existence of UIDs system-wide
%    deleteContainer - Delete a container
%    deleteContainerAnalysis - Delete an anaylsis
%    deleteContainerAnalysisNote - Remove a note from container analysis.
%    deleteContainerFile - Delete a file
%    deleteContainerNote - Remove a note from container
%    deleteContainerTag - Delete a tag
%    downloadContainerAnalysisInputs - Download analysis inputs.
%    downloadContainerAnalysisOutputs - Download analysis outputs.
%    downloadFileFromContainer - Download a file.
%    getContainerFileZipInfo - Download a file.
%    getContainerDownloadUrl - Download a file.
%    downloadInputFromContainerAnalysis - Download analysis inputs with filter.
%    getContainerAnalysisInputZipInfo - Download analysis inputs with filter.
%    getContainerAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromContainerAnalysis - Download analysis outputs with filter.
%    getContainerAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getContainerAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    getContainer - Retrieve a single container
%    getContainerAnalyses - Get analyses for container.
%    getContainerAnalysis - Get an analysis.
%    getContainerFileInfo - Get info for a particular file.
%    getContainerNote - Get a note on container.
%    getContainerTag - Get the value of a tag, by name.
%    modifyContainer - Update a container and its contents
%    modifyContainerAnalysis - Modify an analysis.
%    modifyContainerFile - Modify a file's attributes
%    modifyContainerFileClassification - Update classification for a particular file.
%    modifyContainerFileInfo - Update info for a particular file.
%    modifyContainerInfo - Update or replace info for a container.
%    modifyContainerNote - Update a note on container.
%    renameContainerTag - Rename a tag.
%    replaceContainerFile - Replace a file
%    uploadFileToContainer - Upload a file to container.
%    uploadOutputToContainerAnalysis - Upload an output file to analysis.
%    getSearchQuerySuggestions - Get suggestions for a structured search query
%    getSearchStatus - Get the status of search (Mongo Connector)
%    parseSearchQuery - Parse a structured search query
%    search - Perform a search query
%    cleanPackfiles - Clean up expired upload tokens and invalid token directories.
%    engineUpload - Upload a list of file fields.
%    fetchTree - Query a portion of the flywheel hierarchy, returning only the requested fields.
%    getAuthStatus - Get Login status
%    getConfig - Return public Scitran configuration information
%    getConfigJs - Return public Scitran configuration information in javascript format.
%    getTreeGraph - Get a description of the flywheel hiearchy
%    getVersion - Get server and database schema version info
%    login - Login
%    logout - Log Out
%    lookupPath - Perform path based lookup of a single node in the Flywheel hierarchy
%    resolvePath - Perform path based lookup of nodes in the Flywheel hierarchy
%    createDevice - Create a new device.
%    getAllDevices - List all devices.
%    getAllDevicesStatus - Get status for all known devices.
%    getDevice - Get device details
%    modifyDevice - Update a device
%    regenerateKey - Regenerate device API key
%    updateDevice - Modify a device's type, name, interval, info or set errors.
%    createProjectAet - Create a new DIMSE project AET
%    createServiceAet - Create a new DIMSE service AET
%    deleteProjectAet - Delete a DIMSE project AET
%    deleteServiceAet - Delete a DIMSE service AET
%    getAllProjectAets - List all DIMSE project AETs
%    getAllServiceAets - List all DIMSE services AETs
%    getProjectAet - Get DIMSE project AET
%    getServiceAet - Get DIMSE service by AET or id
%    createDownloadTicket - Create a download ticket
%    downloadTicket - Download files listed in the given ticket.
%    uploadByLabel - Multipart form upload with N file fields, each with their desired filename.
%    uploadByReaper - Bottom-up UID matching of Multipart form upload with N file fields, each with their desired filename.
%    uploadByUid - Multipart form upload with N file fields, each with their desired filename.
%    uploadMatchUid - Multipart form upload with N file fields, each with their desired filename.
%    addGear - Create or update a gear.
%    deleteGear - Delete a gear (not recommended)
%    getAllGears - List all gears
%    getGear - Retrieve details about a specific gear
%    getGearContext - Get context values for the given gear and container.
%    getGearInvocation - Get a schema for invoking a gear.
%    getGearSuggest - Get files with input suggestions, parent containers, and child containers for the given container.
%    getGearTicket - Retrieve a specific gear ticket
%    getMyGearTickets - Retrieve all gear tickets for the current user
%    prepareAddGear - Prepare a gear upload
%    saveGear - Report the result of a gear upload and save the ticket
%    addGroup - Add a group
%    addGroupPermission - Add a permission
%    addGroupTag - Add a tag to group.
%    deleteGroup - Delete a group
%    deleteGroupTag - Delete a tag
%    deleteGroupUserPermission - Delete a permission
%    getAllGroups - List all groups
%    getGroup - Get group info
%    getGroupProjects - Get all projects in a group
%    getGroupTag - Get the value of a tag, by name.
%    getGroupUserPermission - List a user's permissions for this group.
%    modifyGroup - Update group
%    modifyGroupUserPermission - Update a user's permission for this group.
%    renameGroupTag - Rename a tag.
%    acceptFailedOutput - Accept failed job output.
%    addJob - Add a job
%    addJobLogs - Add logs to a job.
%    askJobs - Ask the queue a question
%    completeJob - Complete a job, with information
%    determineProviderForJob - Determine the effective compute provider for a proposed job.
%    getAllJobs - Return all jobs
%    getJob - Get job details
%    getJobConfig - Get a job's config
%    getJobDetail - Get job container details
%    getJobLogs - Get job logs
%    getJobsStats - Get stats about all current jobs
%    getNextJob - Get the next job in the queue
%    modifyJob - Update a job.
%    prepareCompleteJob - Create a ticket for completing a job, with id and status.
%    reapJobs - Reap stale jobs
%    retryJob - Retry a job.
%    updateJobProfile - Update profile information on a job. (e.g. machine type, etc)
%    addModality - Create a new modality.
%    deleteModality - Delete a modality
%    getAllModalities - List all modalities.
%    getModality - Get a modality's classification specification
%    replaceModality - Replace modality
%    addProject - Create a new project
%    addProjectAnalysis - Create an analysis and upload files.
%    addProjectAnalysisNote - Add a note to project analysis.
%    addProjectNote - Add a note to project.
%    addProjectPermission - Add a permission
%    addProjectRule - Create a new rule for a project.
%    addProjectTag - Add a tag to project.
%    deleteProject - Delete a project
%    deleteProjectAnalysis - Delete an anaylsis
%    deleteProjectAnalysisNote - Remove a note from project analysis.
%    deleteProjectFile - Delete a file
%    deleteProjectNote - Remove a note from project
%    deleteProjectTag - Delete a tag
%    deleteProjectUserPermission - Delete a permission
%    downloadFileFromProject - Download a file.
%    getProjectFileZipInfo - Download a file.
%    getProjectDownloadUrl - Download a file.
%    downloadInputFromProjectAnalysis - Download analysis inputs with filter.
%    getProjectAnalysisInputZipInfo - Download analysis inputs with filter.
%    getProjectAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromProjectAnalysis - Download analysis outputs with filter.
%    getProjectAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getProjectAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    downloadProjectAnalysisInputs - Download analysis inputs.
%    downloadProjectAnalysisOutputs - Download analysis outputs.
%    endProjectPackfileUpload - End a packfile upload
%    getAllProjects - Get a list of projects
%    getAllProjectsGroups - List all groups which have a project in them
%    getProject - Get a single project
%    getProjectAcquisitions - List all acquisitions for the given project.
%    getProjectAnalyses - Get analyses for project.
%    getProjectAnalysis - Get an analysis.
%    getProjectFileInfo - Get info for a particular file.
%    getProjectNote - Get a note on project.
%    getProjectRule - Get a project rule.
%    getProjectRules - List all rules for a project.
%    getProjectSessions - List all sessions for the given project.
%    getProjectSubjects - List all subjects for the given project.
%    getProjectTag - Get the value of a tag, by name.
%    getProjectUserPermission - List a user's permissions for this project.
%    modifyProject - Update a project
%    modifyProjectAnalysis - Modify an analysis.
%    modifyProjectFile - Modify a file's attributes
%    modifyProjectFileClassification - Update classification for a particular file.
%    modifyProjectFileInfo - Update info for a particular file.
%    modifyProjectInfo - Update or replace info for a project.
%    modifyProjectNote - Update a note on project.
%    modifyProjectRule - Update a rule on a project.
%    modifyProjectUserPermission - Update a user's permission for this project.
%    projectPackfileUpload - Add files to an in-progress packfile
%    recalcAllProjects - Recalculate all sessions against their project templates.
%    recalcProject - Recalculate if sessions in the project satisfy the template.
%    removeProjectRule - Remove a project rule.
%    removeProjectTemplate - Remove the session template for a project.
%    renameProjectTag - Rename a tag.
%    replaceProjectFile - Replace a file
%    setProjectTemplate - Set the session template for a project.
%    startProjectPackfileUpload - Start a packfile upload to project
%    uploadFileToProject - Upload a file to project.
%    uploadOutputToProjectAnalysis - Upload an output file to analysis.
%    collectUsage - Collect daily usage statistics.
%    getAccessLogReport - Get a report of access log entries for the given parameters
%    getAccessLogTypes - Get the list of types of access log entries
%    getDailyUsageReport - Get a daily usage report for the given month.
%    getLegacyUsageReport - Get a usage report for the site grouped by month or project
%    getProjectReport
%    getSiteReport
%    getUsageAvailability - Get year/month combinations where report data is available.
%    getUsageReport - Get a usage report for the given month.
%    addSiteRule - Create a new site rule.
%    getSiteRule - Get a site rule.
%    getSiteRules - List all site rules.
%    modifySiteRule - Update a site rule.
%    removeSiteRule - Remove a site rule.
%    addSession - Create a new session
%    addSessionAnalysis - Create an analysis and upload files.
%    addSessionAnalysisNote - Add a note to session analysis.
%    addSessionNote - Add a note to session.
%    addSessionTag - Add a tag to session.
%    deleteSession - Delete a session
%    deleteSessionAnalysis - Delete an anaylsis
%    deleteSessionAnalysisNote - Remove a note from session analysis.
%    deleteSessionFile - Delete a file
%    deleteSessionNote - Remove a note from session
%    deleteSessionTag - Delete a tag
%    downloadFileFromSession - Download a file.
%    getSessionFileZipInfo - Download a file.
%    getSessionDownloadUrl - Download a file.
%    downloadInputFromSessionAnalysis - Download analysis inputs with filter.
%    getSessionAnalysisInputZipInfo - Download analysis inputs with filter.
%    getSessionAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromSessionAnalysis - Download analysis outputs with filter.
%    getSessionAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getSessionAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    downloadSessionAnalysisInputs - Download analysis inputs.
%    downloadSessionAnalysisOutputs - Download analysis outputs.
%    getAllSessions - Get a list of sessions
%    getSession - Get a single session
%    getSessionAcquisitions - List acquisitions in a session
%    getSessionAnalyses - Get analyses for session.
%    getSessionAnalysis - Get an analysis.
%    getSessionFileInfo - Get info for a particular file.
%    getSessionJobs - Return any jobs that use inputs from this session
%    getSessionNote - Get a note on session.
%    getSessionTag - Get the value of a tag, by name.
%    modifySession - Update a session
%    modifySessionAnalysis - Modify an analysis.
%    modifySessionFile - Modify a file's attributes
%    modifySessionFileClassification - Update classification for a particular file.
%    modifySessionFileInfo - Update info for a particular file.
%    modifySessionInfo - Update or replace info for a session.
%    modifySessionNote - Update a note on session.
%    renameSessionTag - Rename a tag.
%    replaceSessionFile - Replace a file
%    uploadFileToSession - Upload a file to session.
%    uploadOutputToSessionAnalysis - Upload an output file to analysis.
%    addProvider - Add a new provider
%    getProvider - Return the provider identified by ProviderId
%    getProviderConfig - Return the configuration for provider identified by ProviderId
%    getProviders - Return a list of all providers on the site
%    getSiteSettings - Return administrative site settings
%    modifyProvider - Update the provider identified by ProviderId
%    modifySiteSettings - Update administrative site settings
%    addSubject - Create a new subject
%    addSubjectAnalysis - Create an analysis and upload files.
%    addSubjectAnalysisNote - Add a note to subject analysis.
%    addSubjectNote - Add a note to subject.
%    addSubjectTag - Add a tag to subject.
%    createMasterSubjectCode - Request a master subject code for the given patient
%    deleteSubject - Delete a subject
%    deleteSubjectAnalysis - Delete an anaylsis
%    deleteSubjectAnalysisNote - Remove a note from subject analysis.
%    deleteSubjectFile - Delete a file
%    deleteSubjectNote - Remove a note from subject
%    deleteSubjectTag - Delete a tag
%    downloadFileFromSubject - Download a file.
%    getSubjectFileZipInfo - Download a file.
%    getSubjectDownloadUrl - Download a file.
%    downloadInputFromSubjectAnalysis - Download analysis inputs with filter.
%    getSubjectAnalysisInputZipInfo - Download analysis inputs with filter.
%    getSubjectAnalysisInputDownloadUrl - Download analysis inputs with filter.
%    downloadOutputFromSubjectAnalysis - Download analysis outputs with filter.
%    getSubjectAnalysisOutputZipInfo - Download analysis outputs with filter.
%    getSubjectAnalysisOutputDownloadUrl - Download analysis outputs with filter.
%    downloadSubjectAnalysisInputs - Download analysis inputs.
%    downloadSubjectAnalysisOutputs - Download analysis outputs.
%    getAllSubjects - Get a list of subjects
%    getSubject - Get a single subject
%    getSubjectAnalyses - Get analyses for subject.
%    getSubjectAnalysis - Get an analysis.
%    getSubjectFileInfo - Get info for a particular file.
%    getSubjectNote - Get a note on subject.
%    getSubjectSessions - List sessions of a subject
%    getSubjectTag - Get the value of a tag, by name.
%    modifySubject - Update a subject
%    modifySubjectAnalysis - Modify an analysis.
%    modifySubjectFile - Modify a file's attributes
%    modifySubjectFileClassification - Update classification for a particular file.
%    modifySubjectFileInfo - Update info for a particular file.
%    modifySubjectInfo - Update or replace info for a subject.
%    modifySubjectNote - Update a note on subject.
%    renameSubjectTag - Rename a tag.
%    replaceSubjectFile - Replace a file
%    uploadFileToSubject - Upload a file to subject.
%    uploadOutputToSubjectAnalysis - Upload an output file to analysis.
%    verifyMasterSubjectCode - Verify that the given master subject code exists or not
%    addUser - Add a new user
%    deleteUser - Delete a user
%    getAllUsers - Return a list of all users
%    getCurrentUser - Get information about the current user
%    getCurrentUserAvatar - Get the avatar of the current user
%    getCurrentUserInfo - Get info of the current user
%    getCurrentUserJobs - Return list of jobs created by the current user
%    getUser - Get information about the specified user
%    getUserAcquisitions - Get all acquisitions that belong to the given user.
%    getUserAvatar - Get the avatar of the specified user
%    getUserCollections - Get all collections that belong to the given user.
%    getUserGroups - List all groups the specified user is a member of
%    getUserProjects - Get all projects that belong to the given user.
%    getUserSessions - Get all sessions that belong to the given user.
%    modifyCurrentUserInfo - Update or replace info for the current user.
%    modifyUser - Update the specified user
%    addView - Add a new data view
%    deleteView - Delete a data view
%    evaluateView - Execute a view, returning data in the preferred format.
%    evaluateViewAdhoc - Execute an ad-hoc view, returning data in the preferred format.
%    getView - Return the view identified by ViewId
%    getViewColumns - Return a list of all known column aliases for use in data views
%    getViews - Return a list of all views belonging to container
%    modifyView - Update the view identified by ViewId
%    saveViewDataToContainer - Execute a view, saving data to the target container / file
classdef Flywheel < handle
    % NOTE: This file is auto generated by the swagger code generator program.
    % Do not edit the file manually.
    properties(Constant)
        SDK_VERSION = '11.4.5';
    end
    properties
        apiClient
        acquisitionsApi
        analysesApi
        batchApi
        bulkApi
        callbacksApi
        collectionsApi
        containersApi
        dataexplorerApi
        defaultApi
        devicesApi
        dimseApi
        filesApi
        gearsApi
        groupsApi
        jobsApi
        modalitiesApi
        projectsApi
        reportsApi
        rulesApi
        sessionsApi
        siteApi
        subjectsApi
        usersApi
        viewsApi
        users
        groups
        projects
        sessions
        subjects
        acquisitions
        jobs
        gears
        collections
        checkVersion
    end
    methods
        function obj = Flywheel(apiKey, root, skipVersionCheck, subjectsInResolver)
            obj.apiClient = flywheel.ApiClient(apiKey);

            % Set root mode
            if exist('root', 'var') && root
                fprintf('WARNING: Root mode is deprecated\n');
                obj.apiClient.restClient.addDefaultParameter('root', 'true');
            end

            if exist('skipVersionCheck', 'var') && ~skipVersionCheck
                obj.checkVersion = true;
            else
                skipEnv = getenv('FLYWHEEL_SDK_SKIP_VERSION_CHECK');
                obj.checkVersion = strcmp('0', skipEnv) || strcmpi('false', skipEnv);
            end

            userAgent = sprintf('Flywheel SDK/%s (Matlab %s; %s)', flywheel.Flywheel.SDK_VERSION, version, computer);
            obj.apiClient.restClient.setDefaultHeader('User-Agent', userAgent);

            obj.acquisitionsApi = flywheel.api.AcquisitionsApi(obj.apiClient, obj);
            obj.analysesApi = flywheel.api.AnalysesApi(obj.apiClient, obj);
            obj.batchApi = flywheel.api.BatchApi(obj.apiClient, obj);
            obj.bulkApi = flywheel.api.BulkApi(obj.apiClient, obj);
            obj.callbacksApi = flywheel.api.CallbacksApi(obj.apiClient, obj);
            obj.collectionsApi = flywheel.api.CollectionsApi(obj.apiClient, obj);
            obj.containersApi = flywheel.api.ContainersApi(obj.apiClient, obj);
            obj.dataexplorerApi = flywheel.api.DataexplorerApi(obj.apiClient, obj);
            obj.defaultApi = flywheel.api.DefaultApi(obj.apiClient, obj);
            obj.devicesApi = flywheel.api.DevicesApi(obj.apiClient, obj);
            obj.dimseApi = flywheel.api.DimseApi(obj.apiClient, obj);
            obj.filesApi = flywheel.api.FilesApi(obj.apiClient, obj);
            obj.gearsApi = flywheel.api.GearsApi(obj.apiClient, obj);
            obj.groupsApi = flywheel.api.GroupsApi(obj.apiClient, obj);
            obj.jobsApi = flywheel.api.JobsApi(obj.apiClient, obj);
            obj.modalitiesApi = flywheel.api.ModalitiesApi(obj.apiClient, obj);
            obj.projectsApi = flywheel.api.ProjectsApi(obj.apiClient, obj);
            obj.reportsApi = flywheel.api.ReportsApi(obj.apiClient, obj);
            obj.rulesApi = flywheel.api.RulesApi(obj.apiClient, obj);
            obj.sessionsApi = flywheel.api.SessionsApi(obj.apiClient, obj);
            obj.siteApi = flywheel.api.SiteApi(obj.apiClient, obj);
            obj.subjectsApi = flywheel.api.SubjectsApi(obj.apiClient, obj);
            obj.usersApi = flywheel.api.UsersApi(obj.apiClient, obj);
            obj.viewsApi = flywheel.api.ViewsApi(obj.apiClient, obj);

            % Initialize finders
            obj.users = flywheel.Finder(obj, 'getAllUsers');
            obj.groups = flywheel.Finder(obj, 'getAllGroups');
            obj.projects = flywheel.Finder(obj, 'getAllProjects');
            obj.subjects = flywheel.Finder(obj, 'getAllSubjects');
            obj.sessions = flywheel.Finder(obj, 'getAllSessions');
            obj.acquisitions = flywheel.Finder(obj, 'getAllAcquisitions');
            obj.jobs = flywheel.Finder(obj, 'getAllJobs');
            obj.gears = flywheel.Finder(obj, 'getAllGears');
            obj.collections = flywheel.Finder(obj, 'getAllCollections');

            % Perform version check
            obj.apiClient.setVersionCheckFn(@() obj.performVersionCheck());

            % Enable subjects in resolver by default
            if ~exist('subjectsInResolver', 'var')
                subjectsEnv = getenv('FLYWHEEL_SDK_SUBJECTS_IN_RESOLVER');
                subjectsInResolver = ~(strcmp('0', subjectsEnv) || strcmpi('false', subjectsEnv));
            end

            if subjectsInResolver
                obj.enableFeature('Subject-Container');
            end
            obj.enableFeature('Safe-Redirect');
        end
        function [returnData, resp] = addAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.addAcquisition(varargin{:});
        end
        function [returnData, resp] = addAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.addAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = addAcquisitionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.addAcquisitionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addAcquisitionNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.addAcquisitionNote(varargin{:});
        end
        function [returnData, resp] = addAcquisitionTag(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.addAcquisitionTag(varargin{:});
        end
        function [returnData, resp] = deleteAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisition(varargin{:});
        end
        function [returnData, resp] = deleteAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteAcquisitionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisitionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteAcquisitionFile(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisitionFile(varargin{:});
        end
        function [returnData, resp] = deleteAcquisitionNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisitionNote(varargin{:});
        end
        function [returnData, resp] = deleteAcquisitionTag(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.deleteAcquisitionTag(varargin{:});
        end
        function [returnData, resp] = downloadAcquisitionAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadAcquisitionAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadAcquisitionAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadAcquisitionAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = downloadFileFromAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadFileFromAcquisition(varargin{:});
        end
        function [returnData, resp] = getAcquisitionFileZipInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getAcquisitionDownloadUrl(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadInputFromAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadOutputFromAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = getAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisition(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalyses(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalyses(varargin{:});
        end
        function [returnData, resp] = getAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = getAcquisitionFileInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionFileInfo(varargin{:});
        end
        function [returnData, resp] = getAcquisitionNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionNote(varargin{:});
        end
        function [returnData, resp] = getAcquisitionTag(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAcquisitionTag(varargin{:});
        end
        function [returnData, resp] = getAllAcquisitions(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.getAllAcquisitions(varargin{:});
        end
        function [returnData, resp] = modifyAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisition(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionFile(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFile(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionFileClassification(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileClassification(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionFileInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileInfo(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionInfo(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionInfo(varargin{:});
        end
        function [returnData, resp] = modifyAcquisitionNote(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionNote(varargin{:});
        end
        function [returnData, resp] = renameAcquisitionTag(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.renameAcquisitionTag(varargin{:});
        end
        function [returnData, resp] = replaceAcquisitionFile(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.replaceAcquisitionFile(varargin{:});
        end
        function [returnData, resp] = uploadFileToAcquisition(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.uploadFileToAcquisition(varargin{:});
        end
        function [returnData, resp] = uploadOutputToAcquisitionAnalysis(obj, varargin)
            [returnData, resp] = obj.acquisitionsApi.uploadOutputToAcquisitionAnalysis(varargin{:});
        end
        function [returnData, resp] = addAnalysisNote(obj, varargin)
            [returnData, resp] = obj.analysesApi.addAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addAnalysisTag(obj, varargin)
            [returnData, resp] = obj.analysesApi.addAnalysisTag(varargin{:});
        end
        function [returnData, resp] = deleteAnalysisNote(obj, varargin)
            [returnData, resp] = obj.analysesApi.deleteAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteAnalysisTag(obj, varargin)
            [returnData, resp] = obj.analysesApi.deleteAnalysisTag(varargin{:});
        end
        function [returnData, resp] = downloadAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.analysesApi.downloadAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.analysesApi.downloadAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = downloadInputFromAnalysis(obj, varargin)
            [returnData, resp] = obj.analysesApi.downloadInputFromAnalysis(varargin{:});
        end
        function [returnData, resp] = getAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromAnalysis(obj, varargin)
            [returnData, resp] = obj.analysesApi.downloadOutputFromAnalysis(varargin{:});
        end
        function [returnData, resp] = getAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = getAnalyses(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalyses(varargin{:});
        end
        function [returnData, resp] = getAnalysis(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysis(varargin{:});
        end
        function [returnData, resp] = getAnalysisNote(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisNote(varargin{:});
        end
        function [returnData, resp] = getAnalysisTag(obj, varargin)
            [returnData, resp] = obj.analysesApi.getAnalysisTag(varargin{:});
        end
        function [returnData, resp] = modifyAnalysis(obj, varargin)
            [returnData, resp] = obj.analysesApi.modifyAnalysis(varargin{:});
        end
        function [returnData, resp] = modifyAnalysisInfo(obj, varargin)
            [returnData, resp] = obj.analysesApi.modifyAnalysisInfo(varargin{:});
        end
        function [returnData, resp] = modifyAnalysisNote(obj, varargin)
            [returnData, resp] = obj.analysesApi.modifyAnalysisNote(varargin{:});
        end
        function [returnData, resp] = renameAnalysisTag(obj, varargin)
            [returnData, resp] = obj.analysesApi.renameAnalysisTag(varargin{:});
        end
        function [returnData, resp] = uploadOutputToAnalysis(obj, varargin)
            [returnData, resp] = obj.analysesApi.uploadOutputToAnalysis(varargin{:});
        end
        function [returnData, resp] = cancelBatch(obj, varargin)
            [returnData, resp] = obj.batchApi.cancelBatch(varargin{:});
        end
        function [returnData, resp] = createBatchJobFromJobs(obj, varargin)
            [returnData, resp] = obj.batchApi.createBatchJobFromJobs(varargin{:});
        end
        function [returnData, resp] = getAllBatches(obj, varargin)
            [returnData, resp] = obj.batchApi.getAllBatches(varargin{:});
        end
        function [returnData, resp] = getBatch(obj, varargin)
            [returnData, resp] = obj.batchApi.getBatch(varargin{:});
        end
        function [returnData, resp] = proposeBatch(obj, varargin)
            [returnData, resp] = obj.batchApi.proposeBatch(varargin{:});
        end
        function [returnData, resp] = startBatch(obj, varargin)
            [returnData, resp] = obj.batchApi.startBatch(varargin{:});
        end
        function [returnData, resp] = bulkCopy(obj, varargin)
            [returnData, resp] = obj.bulkApi.bulkCopy(varargin{:});
        end
        function [returnData, resp] = bulkDelete(obj, varargin)
            [returnData, resp] = obj.bulkApi.bulkDelete(varargin{:});
        end
        function [returnData, resp] = bulkMove(obj, varargin)
            [returnData, resp] = obj.bulkApi.bulkMove(varargin{:});
        end
        function [returnData, resp] = bulkMoveSessions(obj, varargin)
            [returnData, resp] = obj.bulkApi.bulkMoveSessions(varargin{:});
        end
        function [returnData, resp] = callbackVirusScan(obj, varargin)
            [returnData, resp] = obj.callbacksApi.callbackVirusScan(varargin{:});
        end
        function [returnData, resp] = addCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollection(varargin{:});
        end
        function [returnData, resp] = addCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = addCollectionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollectionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addCollectionNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollectionNote(varargin{:});
        end
        function [returnData, resp] = addCollectionPermission(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollectionPermission(varargin{:});
        end
        function [returnData, resp] = addCollectionTag(obj, varargin)
            [returnData, resp] = obj.collectionsApi.addCollectionTag(varargin{:});
        end
        function [returnData, resp] = deleteCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollection(varargin{:});
        end
        function [returnData, resp] = deleteCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteCollectionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteCollectionFile(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionFile(varargin{:});
        end
        function [returnData, resp] = deleteCollectionNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionNote(varargin{:});
        end
        function [returnData, resp] = deleteCollectionTag(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionTag(varargin{:});
        end
        function [returnData, resp] = deleteCollectionUserPermission(obj, varargin)
            [returnData, resp] = obj.collectionsApi.deleteCollectionUserPermission(varargin{:});
        end
        function [returnData, resp] = downloadCollectionAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.collectionsApi.downloadCollectionAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadCollectionAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.collectionsApi.downloadCollectionAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = downloadFileFromCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.downloadFileFromCollection(varargin{:});
        end
        function [returnData, resp] = getCollectionFileZipInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getCollectionDownloadUrl(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.downloadInputFromCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.downloadOutputFromCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = getAllCollections(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getAllCollections(varargin{:});
        end
        function [returnData, resp] = getAllCollectionsCurators(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getAllCollectionsCurators(varargin{:});
        end
        function [returnData, resp] = getCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollection(varargin{:});
        end
        function [returnData, resp] = getCollectionAcquisitions(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAcquisitions(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalyses(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalyses(varargin{:});
        end
        function [returnData, resp] = getCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = getCollectionFileInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionFileInfo(varargin{:});
        end
        function [returnData, resp] = getCollectionNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionNote(varargin{:});
        end
        function [returnData, resp] = getCollectionSessions(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionSessions(varargin{:});
        end
        function [returnData, resp] = getCollectionTag(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionTag(varargin{:});
        end
        function [returnData, resp] = getCollectionUserPermission(obj, varargin)
            [returnData, resp] = obj.collectionsApi.getCollectionUserPermission(varargin{:});
        end
        function [returnData, resp] = modifyCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollection(varargin{:});
        end
        function [returnData, resp] = modifyCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = modifyCollectionFile(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionFile(varargin{:});
        end
        function [returnData, resp] = modifyCollectionFileClassification(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileClassification(varargin{:});
        end
        function [returnData, resp] = modifyCollectionFileInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileInfo(varargin{:});
        end
        function [returnData, resp] = modifyCollectionInfo(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionInfo(varargin{:});
        end
        function [returnData, resp] = modifyCollectionNote(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionNote(varargin{:});
        end
        function [returnData, resp] = modifyCollectionUserPermission(obj, varargin)
            [returnData, resp] = obj.collectionsApi.modifyCollectionUserPermission(varargin{:});
        end
        function [returnData, resp] = renameCollectionTag(obj, varargin)
            [returnData, resp] = obj.collectionsApi.renameCollectionTag(varargin{:});
        end
        function [returnData, resp] = replaceCollectionFile(obj, varargin)
            [returnData, resp] = obj.collectionsApi.replaceCollectionFile(varargin{:});
        end
        function [returnData, resp] = uploadFileToCollection(obj, varargin)
            [returnData, resp] = obj.collectionsApi.uploadFileToCollection(varargin{:});
        end
        function [returnData, resp] = uploadOutputToCollectionAnalysis(obj, varargin)
            [returnData, resp] = obj.collectionsApi.uploadOutputToCollectionAnalysis(varargin{:});
        end
        function [returnData, resp] = addContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.addContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = addContainerAnalysisNote(obj, varargin)
            [returnData, resp] = obj.containersApi.addContainerAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addContainerNote(obj, varargin)
            [returnData, resp] = obj.containersApi.addContainerNote(varargin{:});
        end
        function [returnData, resp] = addContainerTag(obj, varargin)
            [returnData, resp] = obj.containersApi.addContainerTag(varargin{:});
        end
        function [returnData, resp] = checkUidsExist(obj, varargin)
            [returnData, resp] = obj.containersApi.checkUidsExist(varargin{:});
        end
        function [returnData, resp] = deleteContainer(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainer(varargin{:});
        end
        function [returnData, resp] = deleteContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteContainerAnalysisNote(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainerAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteContainerFile(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainerFile(varargin{:});
        end
        function [returnData, resp] = deleteContainerNote(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainerNote(varargin{:});
        end
        function [returnData, resp] = deleteContainerTag(obj, varargin)
            [returnData, resp] = obj.containersApi.deleteContainerTag(varargin{:});
        end
        function [returnData, resp] = downloadContainerAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.containersApi.downloadContainerAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadContainerAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.containersApi.downloadContainerAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = downloadFileFromContainer(obj, varargin)
            [returnData, resp] = obj.containersApi.downloadFileFromContainer(varargin{:});
        end
        function [returnData, resp] = getContainerFileZipInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getContainerDownloadUrl(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.downloadInputFromContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = getContainerAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getContainerAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.downloadOutputFromContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = getContainerAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getContainerAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = getContainer(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainer(varargin{:});
        end
        function [returnData, resp] = getContainerAnalyses(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalyses(varargin{:});
        end
        function [returnData, resp] = getContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = getContainerFileInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerFileInfo(varargin{:});
        end
        function [returnData, resp] = getContainerNote(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerNote(varargin{:});
        end
        function [returnData, resp] = getContainerTag(obj, varargin)
            [returnData, resp] = obj.containersApi.getContainerTag(varargin{:});
        end
        function [returnData, resp] = modifyContainer(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainer(varargin{:});
        end
        function [returnData, resp] = modifyContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = modifyContainerFile(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerFile(varargin{:});
        end
        function [returnData, resp] = modifyContainerFileClassification(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerFileClassification(varargin{:});
        end
        function [returnData, resp] = modifyContainerFileInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerFileInfo(varargin{:});
        end
        function [returnData, resp] = modifyContainerInfo(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerInfo(varargin{:});
        end
        function [returnData, resp] = modifyContainerNote(obj, varargin)
            [returnData, resp] = obj.containersApi.modifyContainerNote(varargin{:});
        end
        function [returnData, resp] = renameContainerTag(obj, varargin)
            [returnData, resp] = obj.containersApi.renameContainerTag(varargin{:});
        end
        function [returnData, resp] = replaceContainerFile(obj, varargin)
            [returnData, resp] = obj.containersApi.replaceContainerFile(varargin{:});
        end
        function [returnData, resp] = uploadFileToContainer(obj, varargin)
            [returnData, resp] = obj.containersApi.uploadFileToContainer(varargin{:});
        end
        function [returnData, resp] = uploadOutputToContainerAnalysis(obj, varargin)
            [returnData, resp] = obj.containersApi.uploadOutputToContainerAnalysis(varargin{:});
        end
        function [returnData, resp] = getSearchQuerySuggestions(obj, varargin)
            [returnData, resp] = obj.dataexplorerApi.getSearchQuerySuggestions(varargin{:});
        end
        function [returnData, resp] = getSearchStatus(obj, varargin)
            [returnData, resp] = obj.dataexplorerApi.getSearchStatus(varargin{:});
        end
        function [returnData, resp] = parseSearchQuery(obj, varargin)
            [returnData, resp] = obj.dataexplorerApi.parseSearchQuery(varargin{:});
        end
        function [returnData, resp] = search(obj, varargin)
            [returnData, resp] = obj.dataexplorerApi.search(varargin{:});
        end
        function [returnData, resp] = cleanPackfiles(obj, varargin)
            [returnData, resp] = obj.defaultApi.cleanPackfiles(varargin{:});
        end
        function [returnData, resp] = engineUpload(obj, varargin)
            [returnData, resp] = obj.defaultApi.engineUpload(varargin{:});
        end
        function [returnData, resp] = fetchTree(obj, varargin)
            [returnData, resp] = obj.defaultApi.fetchTree(varargin{:});
        end
        function [returnData, resp] = getAuthStatus(obj, varargin)
            [returnData, resp] = obj.defaultApi.getAuthStatus(varargin{:});
        end
        function [returnData, resp] = getConfig(obj, varargin)
            [returnData, resp] = obj.defaultApi.getConfig(varargin{:});
        end
        function [returnData, resp] = getConfigJs(obj, varargin)
            [returnData, resp] = obj.defaultApi.getConfigJs(varargin{:});
        end
        function [returnData, resp] = getTreeGraph(obj, varargin)
            [returnData, resp] = obj.defaultApi.getTreeGraph(varargin{:});
        end
        function [returnData, resp] = getVersion(obj, varargin)
            [returnData, resp] = obj.defaultApi.getVersion(varargin{:});
        end
        function [returnData, resp] = login(obj, varargin)
            [returnData, resp] = obj.defaultApi.login(varargin{:});
        end
        function [returnData, resp] = logout(obj, varargin)
            [returnData, resp] = obj.defaultApi.logout(varargin{:});
        end
        function [returnData, resp] = lookupPath(obj, varargin)
            [returnData, resp] = obj.defaultApi.lookupPath(varargin{:});
        end
        function [returnData, resp] = resolvePath(obj, varargin)
            [returnData, resp] = obj.defaultApi.resolvePath(varargin{:});
        end
        function [returnData, resp] = createDevice(obj, varargin)
            [returnData, resp] = obj.devicesApi.createDevice(varargin{:});
        end
        function [returnData, resp] = getAllDevices(obj, varargin)
            [returnData, resp] = obj.devicesApi.getAllDevices(varargin{:});
        end
        function [returnData, resp] = getAllDevicesStatus(obj, varargin)
            [returnData, resp] = obj.devicesApi.getAllDevicesStatus(varargin{:});
        end
        function [returnData, resp] = getDevice(obj, varargin)
            [returnData, resp] = obj.devicesApi.getDevice(varargin{:});
        end
        function [returnData, resp] = modifyDevice(obj, varargin)
            [returnData, resp] = obj.devicesApi.modifyDevice(varargin{:});
        end
        function [returnData, resp] = regenerateKey(obj, varargin)
            [returnData, resp] = obj.devicesApi.regenerateKey(varargin{:});
        end
        function [returnData, resp] = updateDevice(obj, varargin)
            [returnData, resp] = obj.devicesApi.updateDevice(varargin{:});
        end
        function [returnData, resp] = createProjectAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.createProjectAet(varargin{:});
        end
        function [returnData, resp] = createServiceAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.createServiceAet(varargin{:});
        end
        function [returnData, resp] = deleteProjectAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.deleteProjectAet(varargin{:});
        end
        function [returnData, resp] = deleteServiceAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.deleteServiceAet(varargin{:});
        end
        function [returnData, resp] = getAllProjectAets(obj, varargin)
            [returnData, resp] = obj.dimseApi.getAllProjectAets(varargin{:});
        end
        function [returnData, resp] = getAllServiceAets(obj, varargin)
            [returnData, resp] = obj.dimseApi.getAllServiceAets(varargin{:});
        end
        function [returnData, resp] = getProjectAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.getProjectAet(varargin{:});
        end
        function [returnData, resp] = getServiceAet(obj, varargin)
            [returnData, resp] = obj.dimseApi.getServiceAet(varargin{:});
        end
        function [returnData, resp] = createDownloadTicket(obj, varargin)
            [returnData, resp] = obj.filesApi.createDownloadTicket(varargin{:});
        end
        function [returnData, resp] = downloadTicket(obj, varargin)
            [returnData, resp] = obj.filesApi.downloadTicket(varargin{:});
        end
        function [returnData, resp] = uploadByLabel(obj, varargin)
            [returnData, resp] = obj.filesApi.uploadByLabel(varargin{:});
        end
        function [returnData, resp] = uploadByReaper(obj, varargin)
            [returnData, resp] = obj.filesApi.uploadByReaper(varargin{:});
        end
        function [returnData, resp] = uploadByUid(obj, varargin)
            [returnData, resp] = obj.filesApi.uploadByUid(varargin{:});
        end
        function [returnData, resp] = uploadMatchUid(obj, varargin)
            [returnData, resp] = obj.filesApi.uploadMatchUid(varargin{:});
        end
        function [returnData, resp] = addGear(obj, varargin)
            [returnData, resp] = obj.gearsApi.addGear(varargin{:});
        end
        function [returnData, resp] = deleteGear(obj, varargin)
            [returnData, resp] = obj.gearsApi.deleteGear(varargin{:});
        end
        function [returnData, resp] = getAllGears(obj, varargin)
            [returnData, resp] = obj.gearsApi.getAllGears(varargin{:});
        end
        function [returnData, resp] = getGear(obj, varargin)
            [returnData, resp] = obj.gearsApi.getGear(varargin{:});
        end
        function [returnData, resp] = getGearContext(obj, varargin)
            [returnData, resp] = obj.gearsApi.getGearContext(varargin{:});
        end
        function [returnData, resp] = getGearInvocation(obj, varargin)
            [returnData, resp] = obj.gearsApi.getGearInvocation(varargin{:});
        end
        function [returnData, resp] = getGearSuggest(obj, varargin)
            [returnData, resp] = obj.gearsApi.getGearSuggest(varargin{:});
        end
        function [returnData, resp] = getGearTicket(obj, varargin)
            [returnData, resp] = obj.gearsApi.getGearTicket(varargin{:});
        end
        function [returnData, resp] = getMyGearTickets(obj, varargin)
            [returnData, resp] = obj.gearsApi.getMyGearTickets(varargin{:});
        end
        function [returnData, resp] = prepareAddGear(obj, varargin)
            [returnData, resp] = obj.gearsApi.prepareAddGear(varargin{:});
        end
        function [returnData, resp] = saveGear(obj, varargin)
            [returnData, resp] = obj.gearsApi.saveGear(varargin{:});
        end
        function [returnData, resp] = addGroup(obj, varargin)
            [returnData, resp] = obj.groupsApi.addGroup(varargin{:});
        end
        function [returnData, resp] = addGroupPermission(obj, varargin)
            [returnData, resp] = obj.groupsApi.addGroupPermission(varargin{:});
        end
        function [returnData, resp] = addGroupTag(obj, varargin)
            [returnData, resp] = obj.groupsApi.addGroupTag(varargin{:});
        end
        function [returnData, resp] = deleteGroup(obj, varargin)
            [returnData, resp] = obj.groupsApi.deleteGroup(varargin{:});
        end
        function [returnData, resp] = deleteGroupTag(obj, varargin)
            [returnData, resp] = obj.groupsApi.deleteGroupTag(varargin{:});
        end
        function [returnData, resp] = deleteGroupUserPermission(obj, varargin)
            [returnData, resp] = obj.groupsApi.deleteGroupUserPermission(varargin{:});
        end
        function [returnData, resp] = getAllGroups(obj, varargin)
            [returnData, resp] = obj.groupsApi.getAllGroups(varargin{:});
        end
        function [returnData, resp] = getGroup(obj, varargin)
            [returnData, resp] = obj.groupsApi.getGroup(varargin{:});
        end
        function [returnData, resp] = getGroupProjects(obj, varargin)
            [returnData, resp] = obj.groupsApi.getGroupProjects(varargin{:});
        end
        function [returnData, resp] = getGroupTag(obj, varargin)
            [returnData, resp] = obj.groupsApi.getGroupTag(varargin{:});
        end
        function [returnData, resp] = getGroupUserPermission(obj, varargin)
            [returnData, resp] = obj.groupsApi.getGroupUserPermission(varargin{:});
        end
        function [returnData, resp] = modifyGroup(obj, varargin)
            [returnData, resp] = obj.groupsApi.modifyGroup(varargin{:});
        end
        function [returnData, resp] = modifyGroupUserPermission(obj, varargin)
            [returnData, resp] = obj.groupsApi.modifyGroupUserPermission(varargin{:});
        end
        function [returnData, resp] = renameGroupTag(obj, varargin)
            [returnData, resp] = obj.groupsApi.renameGroupTag(varargin{:});
        end
        function [returnData, resp] = acceptFailedOutput(obj, varargin)
            [returnData, resp] = obj.jobsApi.acceptFailedOutput(varargin{:});
        end
        function [returnData, resp] = addJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.addJob(varargin{:});
        end
        function [returnData, resp] = addJobLogs(obj, varargin)
            [returnData, resp] = obj.jobsApi.addJobLogs(varargin{:});
        end
        function [returnData, resp] = askJobs(obj, varargin)
            [returnData, resp] = obj.jobsApi.askJobs(varargin{:});
        end
        function [returnData, resp] = completeJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.completeJob(varargin{:});
        end
        function [returnData, resp] = determineProviderForJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.determineProviderForJob(varargin{:});
        end
        function [returnData, resp] = getAllJobs(obj, varargin)
            [returnData, resp] = obj.jobsApi.getAllJobs(varargin{:});
        end
        function [returnData, resp] = getJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.getJob(varargin{:});
        end
        function [returnData, resp] = getJobConfig(obj, varargin)
            [returnData, resp] = obj.jobsApi.getJobConfig(varargin{:});
        end
        function [returnData, resp] = getJobDetail(obj, varargin)
            [returnData, resp] = obj.jobsApi.getJobDetail(varargin{:});
        end
        function [returnData, resp] = getJobLogs(obj, varargin)
            [returnData, resp] = obj.jobsApi.getJobLogs(varargin{:});
        end
        function [returnData, resp] = getJobsStats(obj, varargin)
            [returnData, resp] = obj.jobsApi.getJobsStats(varargin{:});
        end
        function [returnData, resp] = getNextJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.getNextJob(varargin{:});
        end
        function [returnData, resp] = modifyJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.modifyJob(varargin{:});
        end
        function [returnData, resp] = prepareCompleteJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.prepareCompleteJob(varargin{:});
        end
        function [returnData, resp] = reapJobs(obj, varargin)
            [returnData, resp] = obj.jobsApi.reapJobs(varargin{:});
        end
        function [returnData, resp] = retryJob(obj, varargin)
            [returnData, resp] = obj.jobsApi.retryJob(varargin{:});
        end
        function [returnData, resp] = updateJobProfile(obj, varargin)
            [returnData, resp] = obj.jobsApi.updateJobProfile(varargin{:});
        end
        function [returnData, resp] = addModality(obj, varargin)
            [returnData, resp] = obj.modalitiesApi.addModality(varargin{:});
        end
        function [returnData, resp] = deleteModality(obj, varargin)
            [returnData, resp] = obj.modalitiesApi.deleteModality(varargin{:});
        end
        function [returnData, resp] = getAllModalities(obj, varargin)
            [returnData, resp] = obj.modalitiesApi.getAllModalities(varargin{:});
        end
        function [returnData, resp] = getModality(obj, varargin)
            [returnData, resp] = obj.modalitiesApi.getModality(varargin{:});
        end
        function [returnData, resp] = replaceModality(obj, varargin)
            [returnData, resp] = obj.modalitiesApi.replaceModality(varargin{:});
        end
        function [returnData, resp] = addProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProject(varargin{:});
        end
        function [returnData, resp] = addProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = addProjectAnalysisNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addProjectNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectNote(varargin{:});
        end
        function [returnData, resp] = addProjectPermission(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectPermission(varargin{:});
        end
        function [returnData, resp] = addProjectRule(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectRule(varargin{:});
        end
        function [returnData, resp] = addProjectTag(obj, varargin)
            [returnData, resp] = obj.projectsApi.addProjectTag(varargin{:});
        end
        function [returnData, resp] = deleteProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProject(varargin{:});
        end
        function [returnData, resp] = deleteProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteProjectAnalysisNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteProjectFile(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectFile(varargin{:});
        end
        function [returnData, resp] = deleteProjectNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectNote(varargin{:});
        end
        function [returnData, resp] = deleteProjectTag(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectTag(varargin{:});
        end
        function [returnData, resp] = deleteProjectUserPermission(obj, varargin)
            [returnData, resp] = obj.projectsApi.deleteProjectUserPermission(varargin{:});
        end
        function [returnData, resp] = downloadFileFromProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.downloadFileFromProject(varargin{:});
        end
        function [returnData, resp] = getProjectFileZipInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getProjectDownloadUrl(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.downloadInputFromProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getProjectAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getProjectAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.downloadOutputFromProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getProjectAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getProjectAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadProjectAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.projectsApi.downloadProjectAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadProjectAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.projectsApi.downloadProjectAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = endProjectPackfileUpload(obj, varargin)
            [returnData, resp] = obj.projectsApi.endProjectPackfileUpload(varargin{:});
        end
        function [returnData, resp] = getAllProjects(obj, varargin)
            [returnData, resp] = obj.projectsApi.getAllProjects(varargin{:});
        end
        function [returnData, resp] = getAllProjectsGroups(obj, varargin)
            [returnData, resp] = obj.projectsApi.getAllProjectsGroups(varargin{:});
        end
        function [returnData, resp] = getProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProject(varargin{:});
        end
        function [returnData, resp] = getProjectAcquisitions(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAcquisitions(varargin{:});
        end
        function [returnData, resp] = getProjectAnalyses(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalyses(varargin{:});
        end
        function [returnData, resp] = getProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getProjectFileInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectFileInfo(varargin{:});
        end
        function [returnData, resp] = getProjectNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectNote(varargin{:});
        end
        function [returnData, resp] = getProjectRule(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectRule(varargin{:});
        end
        function [returnData, resp] = getProjectRules(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectRules(varargin{:});
        end
        function [returnData, resp] = getProjectSessions(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectSessions(varargin{:});
        end
        function [returnData, resp] = getProjectSubjects(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectSubjects(varargin{:});
        end
        function [returnData, resp] = getProjectTag(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectTag(varargin{:});
        end
        function [returnData, resp] = getProjectUserPermission(obj, varargin)
            [returnData, resp] = obj.projectsApi.getProjectUserPermission(varargin{:});
        end
        function [returnData, resp] = modifyProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProject(varargin{:});
        end
        function [returnData, resp] = modifyProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = modifyProjectFile(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectFile(varargin{:});
        end
        function [returnData, resp] = modifyProjectFileClassification(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectFileClassification(varargin{:});
        end
        function [returnData, resp] = modifyProjectFileInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectFileInfo(varargin{:});
        end
        function [returnData, resp] = modifyProjectInfo(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectInfo(varargin{:});
        end
        function [returnData, resp] = modifyProjectNote(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectNote(varargin{:});
        end
        function [returnData, resp] = modifyProjectRule(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectRule(varargin{:});
        end
        function [returnData, resp] = modifyProjectUserPermission(obj, varargin)
            [returnData, resp] = obj.projectsApi.modifyProjectUserPermission(varargin{:});
        end
        function [returnData, resp] = projectPackfileUpload(obj, varargin)
            [returnData, resp] = obj.projectsApi.projectPackfileUpload(varargin{:});
        end
        function [returnData, resp] = recalcAllProjects(obj, varargin)
            [returnData, resp] = obj.projectsApi.recalcAllProjects(varargin{:});
        end
        function [returnData, resp] = recalcProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.recalcProject(varargin{:});
        end
        function [returnData, resp] = removeProjectRule(obj, varargin)
            [returnData, resp] = obj.projectsApi.removeProjectRule(varargin{:});
        end
        function [returnData, resp] = removeProjectTemplate(obj, varargin)
            [returnData, resp] = obj.projectsApi.removeProjectTemplate(varargin{:});
        end
        function [returnData, resp] = renameProjectTag(obj, varargin)
            [returnData, resp] = obj.projectsApi.renameProjectTag(varargin{:});
        end
        function [returnData, resp] = replaceProjectFile(obj, varargin)
            [returnData, resp] = obj.projectsApi.replaceProjectFile(varargin{:});
        end
        function [returnData, resp] = setProjectTemplate(obj, varargin)
            [returnData, resp] = obj.projectsApi.setProjectTemplate(varargin{:});
        end
        function [returnData, resp] = startProjectPackfileUpload(obj, varargin)
            [returnData, resp] = obj.projectsApi.startProjectPackfileUpload(varargin{:});
        end
        function [returnData, resp] = uploadFileToProject(obj, varargin)
            [returnData, resp] = obj.projectsApi.uploadFileToProject(varargin{:});
        end
        function [returnData, resp] = uploadOutputToProjectAnalysis(obj, varargin)
            [returnData, resp] = obj.projectsApi.uploadOutputToProjectAnalysis(varargin{:});
        end
        function [returnData, resp] = collectUsage(obj, varargin)
            [returnData, resp] = obj.reportsApi.collectUsage(varargin{:});
        end
        function [returnData, resp] = getAccessLogReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getAccessLogReport(varargin{:});
        end
        function [returnData, resp] = getAccessLogTypes(obj, varargin)
            [returnData, resp] = obj.reportsApi.getAccessLogTypes(varargin{:});
        end
        function [returnData, resp] = getDailyUsageReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getDailyUsageReport(varargin{:});
        end
        function [returnData, resp] = getLegacyUsageReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getLegacyUsageReport(varargin{:});
        end
        function [returnData, resp] = getProjectReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getProjectReport(varargin{:});
        end
        function [returnData, resp] = getSiteReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getSiteReport(varargin{:});
        end
        function [returnData, resp] = getUsageAvailability(obj, varargin)
            [returnData, resp] = obj.reportsApi.getUsageAvailability(varargin{:});
        end
        function [returnData, resp] = getUsageReport(obj, varargin)
            [returnData, resp] = obj.reportsApi.getUsageReport(varargin{:});
        end
        function [returnData, resp] = addSiteRule(obj, varargin)
            [returnData, resp] = obj.rulesApi.addSiteRule(varargin{:});
        end
        function [returnData, resp] = getSiteRule(obj, varargin)
            [returnData, resp] = obj.rulesApi.getSiteRule(varargin{:});
        end
        function [returnData, resp] = getSiteRules(obj, varargin)
            [returnData, resp] = obj.rulesApi.getSiteRules(varargin{:});
        end
        function [returnData, resp] = modifySiteRule(obj, varargin)
            [returnData, resp] = obj.rulesApi.modifySiteRule(varargin{:});
        end
        function [returnData, resp] = removeSiteRule(obj, varargin)
            [returnData, resp] = obj.rulesApi.removeSiteRule(varargin{:});
        end
        function [returnData, resp] = addSession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.addSession(varargin{:});
        end
        function [returnData, resp] = addSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.addSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = addSessionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.addSessionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addSessionNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.addSessionNote(varargin{:});
        end
        function [returnData, resp] = addSessionTag(obj, varargin)
            [returnData, resp] = obj.sessionsApi.addSessionTag(varargin{:});
        end
        function [returnData, resp] = deleteSession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSession(varargin{:});
        end
        function [returnData, resp] = deleteSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteSessionAnalysisNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSessionAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteSessionFile(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSessionFile(varargin{:});
        end
        function [returnData, resp] = deleteSessionNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSessionNote(varargin{:});
        end
        function [returnData, resp] = deleteSessionTag(obj, varargin)
            [returnData, resp] = obj.sessionsApi.deleteSessionTag(varargin{:});
        end
        function [returnData, resp] = downloadFileFromSession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.downloadFileFromSession(varargin{:});
        end
        function [returnData, resp] = getSessionFileZipInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getSessionDownloadUrl(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.downloadInputFromSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = getSessionAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getSessionAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.downloadOutputFromSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = getSessionAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getSessionAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadSessionAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.sessionsApi.downloadSessionAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadSessionAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.sessionsApi.downloadSessionAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = getAllSessions(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getAllSessions(varargin{:});
        end
        function [returnData, resp] = getSession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSession(varargin{:});
        end
        function [returnData, resp] = getSessionAcquisitions(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAcquisitions(varargin{:});
        end
        function [returnData, resp] = getSessionAnalyses(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalyses(varargin{:});
        end
        function [returnData, resp] = getSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = getSessionFileInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionFileInfo(varargin{:});
        end
        function [returnData, resp] = getSessionJobs(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionJobs(varargin{:});
        end
        function [returnData, resp] = getSessionNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionNote(varargin{:});
        end
        function [returnData, resp] = getSessionTag(obj, varargin)
            [returnData, resp] = obj.sessionsApi.getSessionTag(varargin{:});
        end
        function [returnData, resp] = modifySession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySession(varargin{:});
        end
        function [returnData, resp] = modifySessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionAnalysis(varargin{:});
        end
        function [returnData, resp] = modifySessionFile(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionFile(varargin{:});
        end
        function [returnData, resp] = modifySessionFileClassification(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionFileClassification(varargin{:});
        end
        function [returnData, resp] = modifySessionFileInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionFileInfo(varargin{:});
        end
        function [returnData, resp] = modifySessionInfo(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionInfo(varargin{:});
        end
        function [returnData, resp] = modifySessionNote(obj, varargin)
            [returnData, resp] = obj.sessionsApi.modifySessionNote(varargin{:});
        end
        function [returnData, resp] = renameSessionTag(obj, varargin)
            [returnData, resp] = obj.sessionsApi.renameSessionTag(varargin{:});
        end
        function [returnData, resp] = replaceSessionFile(obj, varargin)
            [returnData, resp] = obj.sessionsApi.replaceSessionFile(varargin{:});
        end
        function [returnData, resp] = uploadFileToSession(obj, varargin)
            [returnData, resp] = obj.sessionsApi.uploadFileToSession(varargin{:});
        end
        function [returnData, resp] = uploadOutputToSessionAnalysis(obj, varargin)
            [returnData, resp] = obj.sessionsApi.uploadOutputToSessionAnalysis(varargin{:});
        end
        function [returnData, resp] = addProvider(obj, varargin)
            [returnData, resp] = obj.siteApi.addProvider(varargin{:});
        end
        function [returnData, resp] = getProvider(obj, varargin)
            [returnData, resp] = obj.siteApi.getProvider(varargin{:});
        end
        function [returnData, resp] = getProviderConfig(obj, varargin)
            [returnData, resp] = obj.siteApi.getProviderConfig(varargin{:});
        end
        function [returnData, resp] = getProviders(obj, varargin)
            [returnData, resp] = obj.siteApi.getProviders(varargin{:});
        end
        function [returnData, resp] = getSiteSettings(obj, varargin)
            [returnData, resp] = obj.siteApi.getSiteSettings(varargin{:});
        end
        function [returnData, resp] = modifyProvider(obj, varargin)
            [returnData, resp] = obj.siteApi.modifyProvider(varargin{:});
        end
        function [returnData, resp] = modifySiteSettings(obj, varargin)
            [returnData, resp] = obj.siteApi.modifySiteSettings(varargin{:});
        end
        function [returnData, resp] = addSubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.addSubject(varargin{:});
        end
        function [returnData, resp] = addSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.addSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = addSubjectAnalysisNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.addSubjectAnalysisNote(varargin{:});
        end
        function [returnData, resp] = addSubjectNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.addSubjectNote(varargin{:});
        end
        function [returnData, resp] = addSubjectTag(obj, varargin)
            [returnData, resp] = obj.subjectsApi.addSubjectTag(varargin{:});
        end
        function [returnData, resp] = createMasterSubjectCode(obj, varargin)
            [returnData, resp] = obj.subjectsApi.createMasterSubjectCode(varargin{:});
        end
        function [returnData, resp] = deleteSubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubject(varargin{:});
        end
        function [returnData, resp] = deleteSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = deleteSubjectAnalysisNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubjectAnalysisNote(varargin{:});
        end
        function [returnData, resp] = deleteSubjectFile(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubjectFile(varargin{:});
        end
        function [returnData, resp] = deleteSubjectNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubjectNote(varargin{:});
        end
        function [returnData, resp] = deleteSubjectTag(obj, varargin)
            [returnData, resp] = obj.subjectsApi.deleteSubjectTag(varargin{:});
        end
        function [returnData, resp] = downloadFileFromSubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.downloadFileFromSubject(varargin{:});
        end
        function [returnData, resp] = getSubjectFileZipInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectFileZipInfo(varargin{:});
        end
        function [returnData, resp] = getSubjectDownloadUrl(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadInputFromSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.downloadInputFromSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalysisInputZipInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalysisInputZipInfo(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalysisInputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalysisInputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadOutputFromSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.downloadOutputFromSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalysisOutputZipInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalysisOutputZipInfo(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalysisOutputDownloadUrl(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalysisOutputDownloadTicket(varargin{:}, 'ticket', true);
            if ~isempty(returnData)
                reqUrl = resp.getRequestUrl().toCharArray';
                returnData = strcat(reqUrl, '=', returnData.ticket);
            end
        end
        function [returnData, resp] = downloadSubjectAnalysisInputs(obj, varargin)
            [returnData, resp] = obj.subjectsApi.downloadSubjectAnalysisInputs(varargin{:});
        end
        function [returnData, resp] = downloadSubjectAnalysisOutputs(obj, varargin)
            [returnData, resp] = obj.subjectsApi.downloadSubjectAnalysisOutputs(varargin{:});
        end
        function [returnData, resp] = getAllSubjects(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getAllSubjects(varargin{:});
        end
        function [returnData, resp] = getSubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubject(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalyses(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalyses(varargin{:});
        end
        function [returnData, resp] = getSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = getSubjectFileInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectFileInfo(varargin{:});
        end
        function [returnData, resp] = getSubjectNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectNote(varargin{:});
        end
        function [returnData, resp] = getSubjectSessions(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectSessions(varargin{:});
        end
        function [returnData, resp] = getSubjectTag(obj, varargin)
            [returnData, resp] = obj.subjectsApi.getSubjectTag(varargin{:});
        end
        function [returnData, resp] = modifySubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubject(varargin{:});
        end
        function [returnData, resp] = modifySubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = modifySubjectFile(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectFile(varargin{:});
        end
        function [returnData, resp] = modifySubjectFileClassification(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectFileClassification(varargin{:});
        end
        function [returnData, resp] = modifySubjectFileInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectFileInfo(varargin{:});
        end
        function [returnData, resp] = modifySubjectInfo(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectInfo(varargin{:});
        end
        function [returnData, resp] = modifySubjectNote(obj, varargin)
            [returnData, resp] = obj.subjectsApi.modifySubjectNote(varargin{:});
        end
        function [returnData, resp] = renameSubjectTag(obj, varargin)
            [returnData, resp] = obj.subjectsApi.renameSubjectTag(varargin{:});
        end
        function [returnData, resp] = replaceSubjectFile(obj, varargin)
            [returnData, resp] = obj.subjectsApi.replaceSubjectFile(varargin{:});
        end
        function [returnData, resp] = uploadFileToSubject(obj, varargin)
            [returnData, resp] = obj.subjectsApi.uploadFileToSubject(varargin{:});
        end
        function [returnData, resp] = uploadOutputToSubjectAnalysis(obj, varargin)
            [returnData, resp] = obj.subjectsApi.uploadOutputToSubjectAnalysis(varargin{:});
        end
        function [returnData, resp] = verifyMasterSubjectCode(obj, varargin)
            [returnData, resp] = obj.subjectsApi.verifyMasterSubjectCode(varargin{:});
        end
        function [returnData, resp] = addUser(obj, varargin)
            [returnData, resp] = obj.usersApi.addUser(varargin{:});
        end
        function [returnData, resp] = deleteUser(obj, varargin)
            [returnData, resp] = obj.usersApi.deleteUser(varargin{:});
        end
        function [returnData, resp] = getAllUsers(obj, varargin)
            [returnData, resp] = obj.usersApi.getAllUsers(varargin{:});
        end
        function [returnData, resp] = getCurrentUser(obj, varargin)
            [returnData, resp] = obj.usersApi.getCurrentUser(varargin{:});
        end
        function [returnData, resp] = getCurrentUserAvatar(obj, varargin)
            [returnData, resp] = obj.usersApi.getCurrentUserAvatar(varargin{:});
        end
        function [returnData, resp] = getCurrentUserInfo(obj, varargin)
            [returnData, resp] = obj.usersApi.getCurrentUserInfo(varargin{:});
        end
        function [returnData, resp] = getCurrentUserJobs(obj, varargin)
            [returnData, resp] = obj.usersApi.getCurrentUserJobs(varargin{:});
        end
        function [returnData, resp] = getUser(obj, varargin)
            [returnData, resp] = obj.usersApi.getUser(varargin{:});
        end
        function [returnData, resp] = getUserAcquisitions(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserAcquisitions(varargin{:});
        end
        function [returnData, resp] = getUserAvatar(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserAvatar(varargin{:});
        end
        function [returnData, resp] = getUserCollections(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserCollections(varargin{:});
        end
        function [returnData, resp] = getUserGroups(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserGroups(varargin{:});
        end
        function [returnData, resp] = getUserProjects(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserProjects(varargin{:});
        end
        function [returnData, resp] = getUserSessions(obj, varargin)
            [returnData, resp] = obj.usersApi.getUserSessions(varargin{:});
        end
        function [returnData, resp] = modifyCurrentUserInfo(obj, varargin)
            [returnData, resp] = obj.usersApi.modifyCurrentUserInfo(varargin{:});
        end
        function [returnData, resp] = modifyUser(obj, varargin)
            [returnData, resp] = obj.usersApi.modifyUser(varargin{:});
        end
        function [returnData, resp] = addView(obj, varargin)
            [returnData, resp] = obj.viewsApi.addView(varargin{:});
        end
        function [returnData, resp] = deleteView(obj, varargin)
            [returnData, resp] = obj.viewsApi.deleteView(varargin{:});
        end
        function [returnData, resp] = evaluateView(obj, varargin)
            [returnData, resp] = obj.viewsApi.evaluateView(varargin{:});
        end
        function [returnData, resp] = evaluateViewAdhoc(obj, varargin)
            [returnData, resp] = obj.viewsApi.evaluateViewAdhoc(varargin{:});
        end
        function [returnData, resp] = getView(obj, varargin)
            [returnData, resp] = obj.viewsApi.getView(varargin{:});
        end
        function [returnData, resp] = getViewColumns(obj, varargin)
            [returnData, resp] = obj.viewsApi.getViewColumns(varargin{:});
        end
        function [returnData, resp] = getViews(obj, varargin)
            [returnData, resp] = obj.viewsApi.getViews(varargin{:});
        end
        function [returnData, resp] = modifyView(obj, varargin)
            [returnData, resp] = obj.viewsApi.modifyView(varargin{:});
        end
        function [returnData, resp] = saveViewDataToContainer(obj, varargin)
            [returnData, resp] = obj.viewsApi.saveViewDataToContainer(varargin{:});
        end

        function [returnData, resp] = addNodesToCollection(obj, collectionId, level, nodeIds, varargin)
            nodes = cellfun(@(id) flywheel.model.CollectionNode('id', id, 'level', level), ...
                nodeIds, 'UniformOutput', false);
            contents = flywheel.model.CollectionOperation('operation', 'add', 'nodes', nodes);
            update = flywheel.model.Collection('contents', contents);
            [returnData, resp] = obj.collectionsApi.modifyCollection(collectionId, update, varargin{:});
        end

        function [returnData, resp] = addSessionsToCollection(obj, collectionId, sessionIds, varargin)
            [returnData, resp] = obj.addNodesToCollection(collectionId, 'session', sessionIds, varargin{:});
        end

        function [returnData, resp] = addAcquisitionsToCollection(obj, collectionId, acquisitionIds, varargin)
            [returnData, resp] = obj.addNodesToCollection(collectionId, 'acquisition', acquisitionIds, varargin{:});
        end

        function [returnData, resp] = changeJobState(obj, jobId, state, varargin)
            [returnData, resp] = obj.modifyJob(jobId, struct('state', state), varargin{:});
        end

        function [returnData, resp] = get(obj, id, varargin)
            [returnData, resp] = obj.getContainer(id, varargin{:});
        end

        function [returnData, resp] = resolve(obj, path, varargin)
            [returnData, resp] = obj.resolvePath(flywheel.model.ResolverInput('path', strsplit(path, '/')), varargin{:});
        end

        function [returnData, resp] = lookup(obj, path, varargin)
            [returnData, resp] = obj.lookupPath(flywheel.model.ResolverInput('path', strsplit(path, '/')), varargin{:});
        end

        function url = fileUrl(obj, path)
            result = obj.resolve(path);
            tail = result.path{end};
            if ~isprop(tail, 'containerType') || ~strcmp(tail.containerType, 'file')
                throw(MException('ApiClient:apiException', 'Resolved path is not a file!'));
            end
            url = tail.url();
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

            p = inputParser;
            p.StructExpand = false;
            addRequired(p, 'containers');
            addRequired(p, 'destFile');
            addParameter(p, 'includeTypes', []);
            addParameter(p, 'excludeTypes', []);
            parse(p, varargin{:});

            if iscell(p.Results.containers)
                containers = p.Results.containers;
            else
                containers = {p.Results.containers};
            end

            % Extract the list of nodes
            nodes = [];
            for i = 1:numel(containers)
                container = containers{i};
                if ~isprop(container, 'containerType')
                    throw(MException('ApiClient:apiException', 'Unknown container specified!'));
                end
                nodes = horzcat(nodes, {flywheel.model.DownloadNode('level', container.containerType, 'id', container.id)});
            end

            % Setup filters
            typeFilter = [];
            if ~isempty(p.Results.includeTypes) || ~isempty(p.Results.excludeTypes)
                typeFilter = flywheel.model.DownloadFilterDefinition('plus', p.Results.includeTypes, ...
                    'minus', p.Results.excludeTypes);
            end

            downloadFilters = [];
            if ~isempty(typeFilter)
                downloadFilters = {flywheel.model.DownloadFilter('types', typeFilter)};
            end

            % Create download request
            request = flywheel.model.Download('nodes', nodes, 'filters', downloadFilters, 'optional', true);
            summary = obj.createDownloadTicket(request);

            % Perform download
            obj.downloadTicket(summary.ticket, p.Results.destFile);
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
            columns = obj.getViewColumns();
            for i = 1:numel(columns)
                column = columns{i};
                if ~isempty(column.group)
                    coltype = 'group';
                else
                    coltype = column.type;
                end
                fprintf('%s (%s): %s\n', column.name, coltype, column.description);
            end
        end

        function data = readViewData(obj, view, containerId, varargin)
            % Execute a data view against container, and return the view data
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            data = obj.evalView(view, containerId, [], varargin{:});
        end

        function destFile = saveViewData(obj, view, containerId, destFile, varargin)
            % Execute a data view against container, and return the view data
            %
            % Parameters:
            %   view: The view id or instance
            %   containerId: The id of the container to execute the view against
            %   destFile: The destination file path
            destFile = obj.evalView(view, containerId, destFile, varargin{:});
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
            obj.requireVersion(9, 0);
            % webread/webwrite did not support Headers before 2016

            inp = inputParser;
            inp.StructExpand = false;
            addRequired(inp, 'view');
            addRequired(inp, 'containerId');

            addParameter(inp, 'filter', []);
            addParameter(inp, 'skip', []);
            addParameter(inp, 'limit', []);

            parse(inp, varargin{:});

            % Query parameters
            queryParams = { 'containerId', inp.Results.containerId, 'format', 'json-flat' };
            if ~isempty(inp.Results.filter)
                queryParams = [queryParams, 'filter', inp.Results.filter];
            end
            if ~isempty(inp.Results.skip)
                queryParams = [queryParams, 'skip', inp.Results.skip];
            end
            if ~isempty(inp.Results.limit)
                queryParams = [queryParams, 'limit', inp.Results.limit];
            end

            % Get default web options
            opts = obj.getWebOptions();
            opts.ContentType = 'json';

            if ischar(inp.Results.view)
                % Saved view execution
                path = '/views/{ViewId}/data';
                pathParams = { 'ViewId', inp.Results.view };

                opts.RequestMethod = 'GET';

                url = obj.buildUrl(path, pathParams, queryParams);
                result = webread(url, opts);
            else
                % Ad-hoc view execution
                path = '/views/data';
                pathParams = {};

                opts.RequestMethod = 'POST';
                opts.MediaType = 'application/json';

                url = obj.buildUrl(path, pathParams, queryParams);
                body = flywheel.model.DataView.ensureIsInstance(inp.Results.view);
                body = flywheel.ApiClient.encodeJson(body.toJson());
                result = webwrite(url, body, opts);
            end
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
            s = obj.readViewStruct(varargin{:});
            result = struct2table(s);
        end

        function [returnData, resp] = downloadFileFromAcquisitionAsData(obj, acquisitionId, fileName, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadFileFromAcquisition(acquisitionId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromAcquisitionAnalysisAsData(obj, acquisitionId, analysisId, filename, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadInputFromAcquisitionAnalysis(acquisitionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromAcquisitionAnalysisAsData(obj, acquisitionId, analysisId, filename, varargin)
            [returnData, resp] = obj.acquisitionsApi.downloadOutputFromAcquisitionAnalysis(acquisitionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setAcquisitionFileClassification(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileClassification(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceAcquisitionFileClassification(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileClassification(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteAcquisitionFileClassificationFields(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileClassification(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setAcquisitionFileInfo(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileInfo(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceAcquisitionFileInfo(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileInfo(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteAcquisitionFileInfoFields(obj, acquisitionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionFileInfo(acquisitionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setAcquisitionInfo(obj, acquisitionId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionInfo(acquisitionId, body, varargin{:});
        end

        function [returnData, resp] = replaceAcquisitionInfo(obj, acquisitionId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionInfo(acquisitionId, body, varargin{:});
        end

        function [returnData, resp] = deleteAcquisitionInfoFields(obj, acquisitionId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.acquisitionsApi.modifyAcquisitionInfo(acquisitionId, body, varargin{:});
        end

        function [returnData, resp] = downloadInputFromAnalysisAsData(obj, analysisId, filename, varargin)
            [returnData, resp] = obj.analysesApi.downloadInputFromAnalysis(analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromAnalysisAsData(obj, analysisId, filename, varargin)
            [returnData, resp] = obj.analysesApi.downloadOutputFromAnalysis(analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setAnalysisInfo(obj, analysisId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.analysesApi.modifyAnalysisInfo(analysisId, body, varargin{:});
        end

        function [returnData, resp] = replaceAnalysisInfo(obj, analysisId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.analysesApi.modifyAnalysisInfo(analysisId, body, varargin{:});
        end

        function [returnData, resp] = deleteAnalysisInfoFields(obj, analysisId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.analysesApi.modifyAnalysisInfo(analysisId, body, varargin{:});
        end

        function [returnData, resp] = downloadFileFromCollectionAsData(obj, collectionId, fileName, varargin)
            [returnData, resp] = obj.collectionsApi.downloadFileFromCollection(collectionId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromCollectionAnalysisAsData(obj, collectionId, analysisId, filename, varargin)
            [returnData, resp] = obj.collectionsApi.downloadInputFromCollectionAnalysis(collectionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromCollectionAnalysisAsData(obj, collectionId, analysisId, filename, varargin)
            [returnData, resp] = obj.collectionsApi.downloadOutputFromCollectionAnalysis(collectionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setCollectionFileClassification(obj, collectionId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileClassification(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceCollectionFileClassification(obj, collectionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileClassification(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteCollectionFileClassificationFields(obj, collectionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileClassification(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setCollectionFileInfo(obj, collectionId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileInfo(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceCollectionFileInfo(obj, collectionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileInfo(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteCollectionFileInfoFields(obj, collectionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionFileInfo(collectionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setCollectionInfo(obj, collectionId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionInfo(collectionId, body, varargin{:});
        end

        function [returnData, resp] = replaceCollectionInfo(obj, collectionId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionInfo(collectionId, body, varargin{:});
        end

        function [returnData, resp] = deleteCollectionInfoFields(obj, collectionId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.collectionsApi.modifyCollectionInfo(collectionId, body, varargin{:});
        end

        function [returnData, resp] = downloadFileFromContainerAsData(obj, containerId, fileName, varargin)
            [returnData, resp] = obj.containersApi.downloadFileFromContainer(containerId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromContainerAnalysisAsData(obj, containerId, analysisId, filename, varargin)
            [returnData, resp] = obj.containersApi.downloadInputFromContainerAnalysis(containerId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromContainerAnalysisAsData(obj, containerId, analysisId, filename, varargin)
            [returnData, resp] = obj.containersApi.downloadOutputFromContainerAnalysis(containerId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setContainerFileClassification(obj, containerId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileClassification(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceContainerFileClassification(obj, containerId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileClassification(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteContainerFileClassificationFields(obj, containerId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileClassification(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setContainerFileInfo(obj, containerId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileInfo(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceContainerFileInfo(obj, containerId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileInfo(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteContainerFileInfoFields(obj, containerId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.containersApi.modifyContainerFileInfo(containerId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setContainerInfo(obj, containerId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.containersApi.modifyContainerInfo(containerId, body, varargin{:});
        end

        function [returnData, resp] = replaceContainerInfo(obj, containerId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.containersApi.modifyContainerInfo(containerId, body, varargin{:});
        end

        function [returnData, resp] = deleteContainerInfoFields(obj, containerId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.containersApi.modifyContainerInfo(containerId, body, varargin{:});
        end

        function [returnData, resp] = downloadTicketAsData(obj, ticket, varargin)
            [returnData, resp] = obj.filesApi.downloadTicket(ticket, [], varargin{:});
        end

        function [returnData, resp] = downloadFileFromProjectAsData(obj, projectId, fileName, varargin)
            [returnData, resp] = obj.projectsApi.downloadFileFromProject(projectId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromProjectAnalysisAsData(obj, projectId, analysisId, filename, varargin)
            [returnData, resp] = obj.projectsApi.downloadInputFromProjectAnalysis(projectId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromProjectAnalysisAsData(obj, projectId, analysisId, filename, varargin)
            [returnData, resp] = obj.projectsApi.downloadOutputFromProjectAnalysis(projectId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setProjectFileClassification(obj, projectId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileClassification(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceProjectFileClassification(obj, projectId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileClassification(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteProjectFileClassificationFields(obj, projectId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileClassification(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setProjectFileInfo(obj, projectId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileInfo(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceProjectFileInfo(obj, projectId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileInfo(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteProjectFileInfoFields(obj, projectId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.projectsApi.modifyProjectFileInfo(projectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setProjectInfo(obj, projectId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.projectsApi.modifyProjectInfo(projectId, body, varargin{:});
        end

        function [returnData, resp] = replaceProjectInfo(obj, projectId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.projectsApi.modifyProjectInfo(projectId, body, varargin{:});
        end

        function [returnData, resp] = deleteProjectInfoFields(obj, projectId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.projectsApi.modifyProjectInfo(projectId, body, varargin{:});
        end

        function [returnData, resp] = downloadFileFromSessionAsData(obj, sessionId, fileName, varargin)
            [returnData, resp] = obj.sessionsApi.downloadFileFromSession(sessionId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromSessionAnalysisAsData(obj, sessionId, analysisId, filename, varargin)
            [returnData, resp] = obj.sessionsApi.downloadInputFromSessionAnalysis(sessionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromSessionAnalysisAsData(obj, sessionId, analysisId, filename, varargin)
            [returnData, resp] = obj.sessionsApi.downloadOutputFromSessionAnalysis(sessionId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setSessionFileClassification(obj, sessionId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileClassification(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceSessionFileClassification(obj, sessionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileClassification(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteSessionFileClassificationFields(obj, sessionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileClassification(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setSessionFileInfo(obj, sessionId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileInfo(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceSessionFileInfo(obj, sessionId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileInfo(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteSessionFileInfoFields(obj, sessionId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.sessionsApi.modifySessionFileInfo(sessionId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setSessionInfo(obj, sessionId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.sessionsApi.modifySessionInfo(sessionId, body, varargin{:});
        end

        function [returnData, resp] = replaceSessionInfo(obj, sessionId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.sessionsApi.modifySessionInfo(sessionId, body, varargin{:});
        end

        function [returnData, resp] = deleteSessionInfoFields(obj, sessionId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.sessionsApi.modifySessionInfo(sessionId, body, varargin{:});
        end

        function [returnData, resp] = downloadFileFromSubjectAsData(obj, subjectId, fileName, varargin)
            [returnData, resp] = obj.subjectsApi.downloadFileFromSubject(subjectId, fileName, [], varargin{:});
        end

        function [returnData, resp] = downloadInputFromSubjectAnalysisAsData(obj, subjectId, analysisId, filename, varargin)
            [returnData, resp] = obj.subjectsApi.downloadInputFromSubjectAnalysis(subjectId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = downloadOutputFromSubjectAnalysisAsData(obj, subjectId, analysisId, filename, varargin)
            [returnData, resp] = obj.subjectsApi.downloadOutputFromSubjectAnalysis(subjectId, analysisId, filename, [], varargin{:});
        end

        function [returnData, resp] = setSubjectFileClassification(obj, subjectId, fileName, body, varargin)
            body = struct( 'add', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileClassification(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceSubjectFileClassification(obj, subjectId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileClassification(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteSubjectFileClassificationFields(obj, subjectId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileClassification(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setSubjectFileInfo(obj, subjectId, fileName, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileInfo(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = replaceSubjectFileInfo(obj, subjectId, fileName, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileInfo(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = deleteSubjectFileInfoFields(obj, subjectId, fileName, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectFileInfo(subjectId, fileName, body, varargin{:});
        end

        function [returnData, resp] = setSubjectInfo(obj, subjectId, body, varargin)
            body = struct( 'set', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectInfo(subjectId, body, varargin{:});
        end

        function [returnData, resp] = replaceSubjectInfo(obj, subjectId, body, varargin)
            body = struct( 'replace', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectInfo(subjectId, body, varargin{:});
        end

        function [returnData, resp] = deleteSubjectInfoFields(obj, subjectId, body, varargin)
            body = struct( 'delete', body );
            [returnData, resp] = obj.subjectsApi.modifySubjectInfo(subjectId, body, varargin{:});
        end



    end
    methods (Access=private)
        function [returnData, resp] = evalView(obj, varargin)
            inp = inputParser;
            inp.StructExpand = false;
            inp.KeepUnmatched = true;
            addRequired(inp, 'view');
            addRequired(inp, 'containerId');
            addRequired(inp, 'destFile');
            addParameter(inp, 'OutputType', 'double');
            parse(inp, varargin{:});

            fwdArgs = {};
            fwdArgNames = fieldnames(inp.Unmatched);
            for i = 1:numel(fwdArgNames)
                name = fwdArgNames{i};
                fwdArgs = [fwdArgs, {name, inp.Unmatched.(name)}];
            end

            if ischar(inp.Results.view)
                [returnData, resp] = obj.viewsApi.evaluateView(inp.Results.view, inp.Results.containerId, fwdArgs{:});
            else
                [returnData, resp] = obj.viewsApi.evaluateViewAdhoc(inp.Results.containerId, inp.Results.view, fwdArgs{:});
            end

            status = resp.getStatusCode();
            switch num2str(status)
                case '200'
                    destFile = inp.Results.destFile;
                    if ~isempty(destFile)
                        resp.saveResponseBodyToFile(destFile);
                        returnData = destFile;
                    else
                        returnData = resp.getBodyData(inp.Results.OutputType);
                    end
                otherwise
                    returnData = [];
            end
        end
        function requireVersion(obj, major, minor)
            versionParts = strsplit(version, '.');
            majorVersion = str2num(versionParts{1});
            minorVersion = str2num(versionParts{2});
            ok = true;

            if majorVersion == major && minorVersion < minor
                ok = false;
            elseif majorVersion < major
                ok = false;
            end
            if ~ok
                error(sprintf('This feature is not supported on Matlab version %s. (Minimum required version: %d.%d)', ...
                    version, major, minor));
            end
        end
        function opts = getWebOptions(obj)
            headers = {};
            userAgent = [];

            defaultHeaders = obj.apiClient.restClient.getDefaultHeaders();
            for i = 1:numel(defaultHeaders)
                header = cell(defaultHeaders(i));
                if strcmpi('user-agent', header{1})
                    userAgent = header{2};
                else
                    headers = [headers; { header{1} header{2} }];
                end
            end

            opts = weboptions('HeaderFields', headers);
            if ~isempty(userAgent)
                opts.UserAgent = userAgent;
            end
        end
        function url = buildUrl(obj, path, pathParams, queryParams)
            url = char(obj.apiClient.restClient.buildUrl(path, pathParams, queryParams));
        end
        function performVersionCheck(obj)
            sdkVersion = flywheel.Flywheel.SDK_VERSION;
            [sdkMajor, sdkMinor] = obj.parseVersion(sdkVersion);
            releaseVersion = '';
            try
                serverVersion = obj.defaultApi.getVersion();
                releaseVersion = serverVersion.release;
            catch ME
            end

            [releaseMajor, releaseMinor] = obj.parseVersion(releaseVersion);

            % Log conditionals:
            % 1. Client or server version not set
            % 2. Major version mismatch
            % 3. SDK Minor version > Server Minor version (client features not available on server)
            % 4. SDK Minor version < Server Minor version (new features on server)
            showPackageVersion = false;

            if releaseMajor > 0 && sdkMajor > 0
                if sdkMajor ~= releaseMajor
                    if isempty(strfind(sdkVersion, 'dev'))
                        warning('Flywheel:versionMismatch', sprintf('Client version %s does not match server version %s. Please update your client version!', sdkVersion, releaseVersion));
                        showPackageVersion = true;
                    end
                elseif obj.checkVersion
                    if sdkMinor > releaseMinor
                        fprintf('Client version %s is ahead of server version %s. Not all client functionality will be supported by the server.\n', sdkVersion, releaseVersion);
                        showPackageVersion = true;
                    elseif sdkMinor < releaseMinor
                        fprintf('Client version %s is behind of server version %s. Please consider upgrading your client to access all available functionality.\n', sdkVersion, releaseVersion);
                        showPackageVersion = true;
                    end
                end
            elseif obj.checkVersion
                warning('Flywheel:versionUnavailable', 'Client or server version not available! This is an unsupported configuration!');
            end

            if showPackageVersion
                fprintf('Go to https://github.com/flywheel-io/core/releases to find a matching package for version: %s\n', releaseVersion);
            end
        end
        function [major, minor] = parseVersion(obj, s)
            if isempty(s)
                major = 0;
                minor = 0;
            else
                parts = strsplit(s, '[^\d\.]+', 'DelimiterType', 'RegularExpression');
                semver = parts{1};
                parts = strsplit(semver, '.');
                major = str2num(parts{1});
                minor = str2num(parts{2});
            end
        end
        function enableFeature(obj, value)
            features = [];

            defaultHeaders = obj.apiClient.restClient.getDefaultHeaders();
            for i = 1:numel(defaultHeaders)
                header = cell(defaultHeaders(i));
                if strcmpi('X-Accept-Feature', header{1})
                    features = header{2};
                    break
                end
            end

            if ~isempty(features)
                features = sprintf('%s %s', features, value);
            else
                features = value;
            end

            obj.apiClient.restClient.setDefaultHeader('X-Accept-Feature', features);
        end


    end
end
