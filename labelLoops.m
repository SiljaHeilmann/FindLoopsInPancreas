function [IMloops] = labelLoops(IMskel,G,cycles)

% this function takes a BW skeletonizes image (IMskel), a graph (G) made
% from that skel and a list of cycles (example found with MATLAB function
% cyclebasis) that you want to label in the returned image IMloops
% Since loops may overlap we may need to use one of the 26 neiboring voxels of the voxel in the skel
% structure once that one has been labelled. (Approx. 26-2=24 neiboring voxel should be available max sqrt(3)=1.7321 pixel lengths away for multible labels)
% so this will only fail when more than 24 loops overlap (which should be
% very rare and only occur for a small fraction of a loop)... meaning
% extraction of example long axis will still give a good enough result...

IMloops = uint16(IMskel);% make uint16 from start
IMloops(IMloops>0)=0;% set to zero

% make array containing all the possible combination of -1 0 1 that will
% give the 26 neigbors of point (0 0 0) in a 3x3x3 cube
addToCoor = zeros(3*3*3,3); % size 27x3 - first row is [0 0 0]
possibilities = [0, -1,1];
count = 0;
for ii=1:3
    for jj=1:3
        for kk=1:3
            count = count+1;
            addToCoor(count,1:3)=[possibilities(ii),possibilities(jj),possibilities(kk)];
        end
    end
end

[maxY,maxX,maxZ]=size(IMloops);

tic
for ii=1:length(cycles)
    curLoopNodeTable = G.Nodes(cycles{ii},:);% make a table with only the Nodes in current cycle
    for jj = 1:length(cycles{ii})% go through all nodes in loop/cycle
        flag = 0; % set flag to 1 when a node has been labelled (either focal node or a neigbor node)
        count = 0; % in case we do have more than 24 loops overlapping we still want to leave while loop - we know we failed when count exceeds 26
        while flag==0 & count<=26
            count = count+1;
            curNodeCoor = curLoopNodeTable.PixelCoor(jj,:) + addToCoor(count,:);% first row of addToCoor is 0 0 0
            curNodeCoor = [max(1,min(curNodeCoor(2),maxY)) max(1,min(curNodeCoor(1),maxX)) max(1,min(curNodeCoor(3),maxZ))];% Make sure coor are not outside image
            if IMloops(curNodeCoor(1),curNodeCoor(2),curNodeCoor(3))==0 & flag==0
                IMloops(curNodeCoor(1),curNodeCoor(2),curNodeCoor(3))=ii;
                flag=1;
            end
        end
     end
end  
end