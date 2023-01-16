%% get image info
fileName = 'explant1_loop1.lsm';
path = '/Users/xqk313/Documents/TemporaryBigFiles';

read = bfGetReader([path filesep fileName]);

omeMeta = read.getMetadataStore();
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % number of time points
stackSizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of Channels

display(['stackSizeX is ' num2str(stackSizeX) ', stackSizeY is ' num2str(stackSizeY) ', stackSizeZ is ' ...
    num2str(stackSizeZ) ', stackSizeT is ' num2str(stackSizeT) ', stackSizeC is ' num2str(stackSizeC)])

%% load in image
%ALLOCATE SPACE FOR STACK
RAW = uint16(zeros(stackSizeY,stackSizeX,stackSizeZ,4,stackSizeT));

for cc=1:stackSizeC % read in all channels and save as mat files (channel order is 1 membrane (red), 2 DIC (gray), 3 lumen/ZO1 (green) )
cc
    for tt=1:stackSizeT
        tt
        for zz=1:stackSizeZ
            zz
            iPlane = read.getIndex(zz - 1, cc - 1, tt - 1) + 1;
            RAW(:,:,zz,cc,tt) = bfGetPlane(read,iPlane);   
        end
    end   
end


size(RAW)
max(RAW(:))

%% Get first look
figure

for cc=1:stackSizeC
    nexttile
    imshow(max(RAW(:,:,:,cc),[],3),[])
    nexttile
    imshow(max(flipxz(RAW(:,:,:,cc)),[],3),[])
end

%% enhance edges in all channels using loclap3 (prepare to make super voxel image) 
C = zeros(size(RAW));
sigma = 0.5;
alpha = 1;
beta = 1;
tic
for cc=1:stackSizeC
    C(:,:,:,cc) = locLap3(imgaussfilt3(RAW(:,:,:,cc)),sigma,alpha,beta);
end
toc

close all
figure
for cc=1:stackSizeC
    nexttile
    imshow(max(C(:,:,:,cc),[],3),[])
    nexttile
    imshow(max(flipxz(C(:,:,:,cc)),[],3),[])
end
%%
G = zeros(size(RAW));
sigma = 0.5;
alpha = 1;
beta = 1;

tic
for cc=1:stackSizeC
    G(:,:,:,cc) = imgradient3(imgaussfilt3(C(:,:,:,cc),1));%,'method','sobel'
end
toc

close all
figure
for cc=1:stackSizeC
    nexttile
    imshow(max(G(:,:,:,cc),[],3),[])
    nexttile
    imshow(max(flipxz(G(:,:,:,cc)),[],3),[])
end
%%
tic
RG = imregionalmax(imgaussfilt3(C(:,:,:,1),[2 2 1]),26);
toc

implaypair(RG,C(:,:,:,1).*15)


%%

Edge1 = imclose(imgaussfilt3(double(otsu3D(C(:,:,:,2),0.1)),1)>0.5,strel('sphere',3));

%RG(Edge1==1)=0; %remove regional maxes situated on clear edges...

RG=imclose(RG,strel('sphere',5));
RG=imclose(RG,strel('disk',10));

RG(Edge1==1)=0; %remove regional maxes situated on clear edges...

implaypair(RG,G(:,:,:,2).*15)

%% watershed C1

B = imimposemin(imgaussfilt3(G(:,:,:,2),1),RG);% catch basin image (impose minima at RG==1 on edge inhanced image)

tic
W = watershed(B,6);% connectivity 6 is lowest possible for 3D im, default is maximal which is 26 for 3D
toc

implaypair(W==0,uint16(G(:,:,:,2)).*15)


%%

implaypair(flipxz(W==0),flipxz(uint16(G(:,:,:,2)).*15))

%%

vtkwrite('watershedlines.vtk', 'structured_points', 'wa', uint8(W==0));
vtkwrite('C1.vtk', 'structured_points', 'C1', uint16(C(:,:,:,1)));
vtkwrite('C2.vtk', 'structured_points', 'C2', uint16(C(:,:,:,2)));

%%
headerInfo = nhdr_nrrd_read('/Volumes/GoogleDrive/My Drive/Segmentation.nrrd', 1);

AN = uint8(headerInfo.data);

%% extend handdrawn annotation to edge of supervoxel



numClasses = max(AN(:))

%%
implaypair(AN,C(:,:,:,1).*15)

ANex = AN;

for nn = 1:numClasses
    TEMP1 = AN==nn;% isolate annotations for just one class
    TEMP2 = W;
    TEMP2(TEMP1==0) = 0;
    superVoxelIndW = nonzeros(unique(TEMP2));
    for ii = 1:length(superVoxelIndW)
        ANex(W==superVoxelIndW(ii))=nn;
    end
end

%%
implaypair(flipxz(ANex),flipxz(C(:,:,:,1).*55))

%%
implaypair(ANex==1,C(:,:,:,2).*15)
