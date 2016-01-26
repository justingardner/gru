% makeSearchROI.m
%
%        $Id: makeSearchROI.m,v 1.1 2008/12/28 15:16:55 justin Exp $ 
%      usage: makeSearchROI(v,surfFilename,anatFilename)
%         by: justin gardner & taosheng liu
%       date: 11/25/08
%    purpose: 
%
function rois = makeSearchROI(v,surfFilename,anatFilename,roiName,varargin)

% check arguments
if nargin < 4
  help makeSearchROI
  return
end

% radius of the search light in mm
searchRadius = [];
getArgs(varargin,{'searchRadius=30'});

% if this debug flag is set, will display the ROIs in the current mrLoadRet view
debug = 1;

% load the anat
anatPath = fileparts(anatFilename);
anatName = getLastDir(anatFilename);
v = loadAnat(v,anatName,anatPath);
anatHdr = viewGet(v,'baseHdr');

% get the roi voxel coordinates
roiVoxelCoords=getROICoordinates(v,roiName,0); 

% load the surface
surf=loadSurfOFF(surfFilename);
surf = xformSurfaceWorld2Array(surf,anatHdr);

% create distanceMatrix
mesh.uniqueVertices = surf.vtcs;
mesh.uniqueFaceIndexList = surf.tris;
mesh.connectionMatrix = findConnectionMatrix(mesh);
surf.distanceMatrix = find3DNeighbourDists(mesh);

for voxnum = 1:size(roiVoxelCoords,2)
  % get the closest vertex in the surface to this coordinate of the roi
  startVertexNum = assignToNearest(surf.vtcs,roiVoxelCoords(:,voxnum)');

  % now compute the distance from this vertex to every other vertex in the surface
  distance = dijkstra(surf.distanceMatrix,startVertexNum);

  % find all the vertices that are with in the searchRadius 
  verticesWithinSearchRadius= find(distance<searchRadius);
  
  % and put the coordinates of those vertices into our search ROI
  searchROI.coords = surf.vtcs(verticesWithinSearchRadius,:)';
  searchROI.coords(4,:) = 1;

  % set other parts of the search ROI
  searchROI.name = sprintf('searchROI_%i_radius_%i',voxnum,searchRadius);
  searchROI.voxelSize=[1 1 1];
  searchROI.xform=viewGet(v,'basesform');
  [tf searchROI] = isroi(searchROI);
  
  % display in current view
  if debug
    viewSet(getMLRView,'newROI',searchROI,1);
    refreshMLRDisplay(viewGet(getMLRView,'viewNum'));
  end

  % save the ROI
%  v = viewSet(v,'newROI',searchROI,1);
%  saveROI(v,searchROI.name);

  rois{voxnum} = searchROI;
end

deleteView(v);


