function G = skel2graph(IMskel)

r = regionprops3(IMskel,'VoxelList','VoxelIdxList');

% MAKE ALL VOXELS IN SKELETON NODES AND PUT EDGES BETWEEN THE ONES THAT ARE NEIGBORS (26 - neigborhood)

coor = r.VoxelList{1};% Get coordinates of all voxels in skeleton
IDs  = r.VoxelIdxList{1}; % get their liniar index in image

% put this information in a Node table (each voxel in skel is a node)
Nodes = table((1:size(coor,1))',coor,IDs);
Nodes.Properties.VariableNames = {'Index','PixelCoor','LiniarIndexInImage'};

Edges = table; % make empty table for edges

for ii=1:ceil(size(coor,1)/2) % Do circshift length(coor)-1 times to go through all possible node combinations
    clear TEMP % clear temp table 

    dist = sqrt(sum((circshift(coor,ii,1)-coor).^2,2)); % calculate distance between coor and circshift(coor,ii,1)

    TEMP = table([(1:size(coor,1))' , circshift((1:size(coor,1))',ii,1)],dist); % Make temp table with endnode index and dist

    TEMP.Properties.VariableNames = {'EndNodes','Dist'}; % give variables meaningfull names
    
    Edges = [Edges; TEMP(TEMP.Dist<=sqrt(3),:)]; % concatenate to build Edge table but add only node pairs where dist is less than sqrt(3) (distance to diagonal neigbor in 26 connected neigborhood)
end


% Make graph
G = graph(Edges,Nodes,"omitselfloops");



end