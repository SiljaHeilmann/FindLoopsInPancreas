function [out] = holefill3D(in)
% remove all 2D black holes in a 3D BW image

out = in;

stackSizeZ = size(in,3);
if stackSizeZ==1
    stackSizeZ = size(in,4);    
    for zz = 1:stackSizeZ    
       out(:,:,:,zz) = imfill(in(:,:,:,zz),'holes');    
    end
else
    for zz = 1:stackSizeZ    
       out(:,:,zz) = imfill(in(:,:,zz),'holes');    
    end
    
end


end
