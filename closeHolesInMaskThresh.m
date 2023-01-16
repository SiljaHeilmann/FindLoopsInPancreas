function [BWout] = closeHolesInMaskThresh(BWin,MASK,T)
% Close holes smaller than T inside MASK

BW_RemoveSmallHoles = logical(neg(BWin));

BW_RemoveSmallHoles = imopen(BW_RemoveSmallHoles,strel('disk',1));

for zz=1:size(BWin,3)% remove small elements 
    BW_RemoveSmallHoles(:,:,zz) = bwareaopen(BW_RemoveSmallHoles(:,:,zz),T,4);    
end

BW_JustSmallHoles = logical(neg(BWin)) - BW_RemoveSmallHoles;% some of these may have to be saved

MASK_holesToSave = logical(neg(MASK));% - imclearborder(neg(loclapBW));### I NEED TO CLEAR BORDER ONE SLICE AT A TIME!!!

for zz=1:size(BWin,3)% remove elements touching border slice by slice
    MASK_holesToSave(:,:,zz) = imclearborder(MASK_holesToSave(:,:,zz),4);   
    BW_JustSmallHoles(:,:,zz) = imclearborder(BW_JustSmallHoles(:,:,zz),4);   
%     if sum(sum(MASK_noBorderObjects(:,:,zz)))==0% if everything got removed - put it back
%    MASK_holesToSave(:,:,zz) = neg(MASK(:,:,zz));        
%     end
end

MASK_holesToSave = imdilate(MASK_holesToSave,strel('sphere',2));
MASK_holesToSave = imdilate(MASK_holesToSave,strel('disk',1));
MASK_holesToSave = imclose(MASK_holesToSave,strel('sphere',1));
MASK_holesToSave = imdilate(MASK_holesToSave,strel('cuboid',[1 1 7]));

%implaypair(BW_JustSmallHoles,MASK_holesToSave)


BW_JustSmallHoles(MASK_holesToSave==1)=0;

BW_JustSmallHoles = imdilate(BW_JustSmallHoles,strel('disk',1));

BWout = BWin;
BWout(BW_JustSmallHoles==1)=1;

end

