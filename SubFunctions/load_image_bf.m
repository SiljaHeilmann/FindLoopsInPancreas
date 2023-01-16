function [RAW] = load_image_bf(fileName)
% IMPORT tiff

readOME = bfGetReader([fileName]);

omeMeta = readOME.getMetadataStore();
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % number of time points
stackSizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of Channels

display(['stackSizeX is ' num2str(stackSizeX) ', stackSizeY is ' num2str(stackSizeY) ', stackSizeZ is ' ...
    num2str(stackSizeZ) ', stackSizeT is ' num2str(stackSizeT) ', stackSizeC is ' num2str(stackSizeC)])

%ALLOCATE SPACE FOR STACK
RAW = uint16(zeros(stackSizeY,stackSizeX,stackSizeZ,stackSizeC,stackSizeT));

for cc=1:stackSizeC % readOME in all channels and save as mat files (channel order is 1 membrane (red), 2 DIC (gray), 3 lumen/ZO1 (green) )
%cc
    for tt=1:stackSizeT
%        tt
        for zz=1:stackSizeZ
%            zz
            iPlane = readOME.getIndex(zz - 1, cc - 1, tt - 1) + 1;
            RAW(:,:,zz,cc,tt) = bfGetPlane(readOME,iPlane);
        end
    end   
end


end

