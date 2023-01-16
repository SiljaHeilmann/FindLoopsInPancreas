close all

path =['/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles'];
%    '/Volumes/GoogleDrive-105399969456537266325/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles'];%'/Volumes/GoogleDrive-105399969456537266325/My Drive/Piljabi/QuantifyLoopInWM/CroppedFiles';%'/Users/xqk313/Google Drive/Piljabi/WholeMountsForLoopQuanitification';% '/Users/xqk313/Documents/Piljabi/WholeMountsForLoopQuantification';
%'/Volumes/GoogleDrive/My Drive/Curvature2021';

filelist = dir([path  filesep  '*.ims'])
%                 ff  1       2       3      4      5      6      7
YVoxelSizeInMicron = [0.568   0.568   0.757  0.461  0.692  0.568   0.598]; 
XVoxelSizeInMicron = [0.568   0.568   0.757  0.461  0.692  0.568   0.843];
ZVoxelSizeInMicron = [2       2       2      1.4    1      2       2    ];

%%


for ff = 1:length(filelist)
    ff
    filelist(ff).name
    read = bfGetReader([path filesep filelist(ff).name ]);
    
    omeMeta = read.getMetadataStore();
    stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % number of time points
    stackSizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of Channels
    
    display(['stackSizeX is ' num2str(stackSizeX) ', stackSizeY is ' num2str(stackSizeY) ', stackSizeZ is ' ...
        num2str(stackSizeZ) ', stackSizeT is ' num2str(stackSizeT) ', stackSizeC is ' num2str(stackSizeC)])
    
    filelist(ff).stackSize       = [stackSizeY , stackSizeX, stackSizeZ];
    filelist(ff).voxelSizeMicron = [YVoxelSizeInMicron(ff) , XVoxelSizeInMicron(ff), ZVoxelSizeInMicron(ff)];
    filelist(ff).reader = read;
    
end

filelist
%%
save('filelist.mat','filelist')
%%
%ALLOCATE SPACE FOR STACK

for ff = 1:length(filelist)
    ff
    %read = filelist(ff).reader;
    read = bfGetReader([path filesep filelist(ff).name ]);
    
    omeMeta = read.getMetadataStore();
    stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    
    % ff          1 2 3 4 5 6 7
    stackSizeC = [3 3 3 2 2 3 3];
    RAW1 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ));
    
    if stackSizeC(ff)==3 % ins channel exist and we want to subtract from Muc1 because of bleed through

        RAW3 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ));
        
        for cc=1:stackSizeC(ff)
            cc
            for zz=1:stackSizeZ
                if cc==1
                    iPlane = read.getIndex(zz - 1, cc - 1, 1 - 1) + 1;
                    RAW1(:,:,zz) = bfGetPlane(read,iPlane);
                elseif cc==3
                    iPlane = read.getIndex(zz - 1, cc - 1, 1 - 1) + 1;
                    RAW3(:,:,zz) = bfGetPlane(read,iPlane);
                end
            end
        end
        RAW1 = RAW1-RAW3;
        save(['RAW1_ff' num2str(ff)   '.mat'],'RAW1');
        save(['RAW3_ff' num2str(ff)   '.mat'],'RAW3');

        clear RAW3
        RAW2 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ));
        for cc=1:stackSizeC(ff)
            cc
            for zz=1:stackSizeZ
                if cc==2
                    iPlane = read.getIndex(zz - 1, cc - 1, 1 - 1) + 1;
                    RAW2(:,:,zz) = bfGetPlane(read,iPlane);
                end
            end
        end
        save(['RAW2_ff' num2str(ff)   '.mat'],'RAW2');

    elseif stackSizeC(ff)==2 % no insulin channel just import channel 1 and 2
        clear RAW3
        RAW2 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ));
        for cc=1:stackSizeC(ff)
            cc
            for zz=1:stackSizeZ
                iPlane = read.getIndex(zz - 1, cc - 1, 1 - 1) + 1;
                if cc==1
                    RAW2(:,:,zz) = bfGetPlane(read,iPlane);% channel order is different for 4 and 5 who has only two channels
                elseif cc==2
                    RAW1(:,:,zz) = bfGetPlane(read,iPlane);
                end
            end
        end
        save(['RAW1_ff' num2str(ff) '.mat'],'RAW1');
        save(['RAW2_ff' num2str(ff) '.mat'],'RAW2');
        
    end
       
    clear RAW2

end

%% Look at raw data

%ff=1

clear RAW1 RAW2 RAW3

%implaypair(RAW1(:,:,100:110),RAW2(:,:,100:110))

for ff = 1:length(filelist)
    ff

    load(['RAW2_ff' num2str(ff)   '.mat']);
    figure
    imshow(max(RAW2(:,:,round(size(RAW2,3)/2)-10:round(size(RAW2,3)/2)+10),[],3))
end


%% DOWNSCALE SIZE TO MAKE ISOTROPIC - one micron per pixel
% and normalize with respect to 99.9999 percentile

for ff = 1:length(filelist)
    ff
    yscale = filelist(ff).voxelSizeMicron(1);
    xscale = filelist(ff).voxelSizeMicron(2);
    zscale = filelist(ff).voxelSizeMicron(3);

    if (ff==4 || ff==5)
        endcc=2;
    else
        endcc=3;
    end

    for cc = 1:endcc

        load(['RAW' num2str(cc) '_ff' num2str(ff)   '.mat']);

        if cc==1
            C1 = imresize3(RAW1(:,:,:),[round(size(RAW1,1).*yscale) round(size(RAW1,2).*xscale) round(size(RAW1,3).*zscale)]);
            C1 = uint8(255.*(double(C1)./double(prctile(C1(:),99.9999))));
            save(['ISO' num2str(cc) '_ff' num2str(ff)   '.mat'],'C1');
            clear C1
        elseif cc==2
            C2 = imresize3(RAW2(:,:,:),[round(size(RAW2,1).*yscale) round(size(RAW2,2).*xscale) round(size(RAW2,3).*zscale)]);
            C2 = uint8(255.*(double(C2)./double(prctile(C2(:),99.9999))));
            save(['ISO' num2str(cc) '_ff' num2str(ff)   '.mat'],'C2');
            clear C2
        elseif cc==3
            C3 = imresize3(RAW3(:,:,:),[round(size(RAW3,1).*yscale) round(size(RAW3,2).*xscale) round(size(RAW3,3).*zscale)]);
            C3 = uint8(255.*(double(C3)./double(prctile(C3(:),99.9999))));
            save(['ISO' num2str(cc) '_ff' num2str(ff)   '.mat'],'C3');
            clear C3
        end

    end
end


%% LOOK AT ISOTROPIC NORMALIZED DATA
close all

for ff = 1:length(filelist)
    ff
    load(['ISO' num2str(3) '_ff' num2str(ff)   '.mat']);
    figure
    imshow(max(C3(:,:,round(size(C3,3)/2)-10:round(size(C3,3)/2)+10),[],3))
end

%% save as h5 files


for ff = 1:length(filelist)
    ff
    if (ff==4 || ff==5)
        load(['ISO1_ff' num2str(ff)   '.mat']);
        load(['ISO2_ff' num2str(ff)   '.mat']);

        IM = cat(4,C1,C2);
        filename = ['ff=' num2str(ff) '_DownScaled.h5'];

        group = [ '/' num2str(ff) ];

        h5create(filename,group ,size(IM),'ChunkSize',[10 10 10 2],'DataType','uint8')
        h5write(filename,group,IM)

        clear IM C1 C2
    else
        load(['ISO1_ff' num2str(ff)   '.mat']);
        load(['ISO2_ff' num2str(ff)   '.mat']);
        load(['ISO3_ff' num2str(ff)   '.mat']);

        IM = cat(4,C1,C2,C3);
        filename = ['ff=' num2str(ff) '_DownScaled.h5'];

        group = [ '/' num2str(ff) ];

        h5create(filename,group ,size(IM),'ChunkSize',[10 10 10 3],'DataType','uint8')
        h5write(filename,group,IM)

        clear IM C1 C2 C3
    end
end

%% save ff=1 as smaller h5 files

for ff = 1:1%2:length(filelist)
    
    load(['ISO1_ff' num2str(ff)   '.mat']);
    load(['ISO2_ff' num2str(ff)   '.mat']);
    load(['ISO3_ff' num2str(ff)   '.mat']);
        
    IM = cat(4,C1,C2,C3);
    size(IM)
    
    IM = IM(round(size(IM,1)/2)-round(size(IM,1)/4):round(size(IM,1)/2)+round(size(IM,1)/4),round(size(IM,2)/2)-round(size(IM,2)/4):round(size(IM,2)/2)+round(size(IM,2)/4),round(size(IM,3)/2)-round(size(IM,3)/4):round(size(IM,3)/2)+round(size(IM,3)/4),:);
    size(IM)
    
    filename = ['ff=' num2str(ff) '_DownScaled_Small.h5'];

    group = [ '/' num2str(ff) ];

    h5create(filename,group ,size(IM),'ChunkSize',[10 10 10 3],'DataType','uint8')
    h5write(filename,group,IM)
    
    clear IM C1 C2 C3
end
%% save as smaller h5 files

for ff = 2:length(filelist)
    
    load(['ISO1_ff' num2str(ff)   '.mat']);
    load(['ISO2_ff' num2str(ff)   '.mat']);
        
    IM = cat(4,C1,C2);
    size(IM)
    
    IM = IM(round(size(IM,1)/2)-round(size(IM,1)/4):round(size(IM,1)/2)+round(size(IM,1)/4),round(size(IM,2)/2)-round(size(IM,2)/4):round(size(IM,2)/2)+round(size(IM,2)/4),round(size(IM,3)/2)-round(size(IM,3)/4):round(size(IM,3)/2)+round(size(IM,3)/4),:);
    size(IM)
    
    filename = [filelist(1).name(1:end-4) '_ff=' num2str(ff) '_DownScaled_Small.h5'];

    group = [ '/' num2str(ff) ];

    h5create(filename,group ,size(IM),'ChunkSize',[10 10 10 2],'DataType','uint8')
    h5write(filename,group,IM)
    
    clear IM C1 C2
end

%% save as smaller vtk files (for annotation in itksnap)

for ff = 1:1%2:length(filelist)
    
    load(['/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles/1h5WholeFiles/ISO1_ff' num2str(ff)   '.mat']);
    load(['/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles/1h5WholeFiles/ISO2_ff' num2str(ff)   '.mat']);
    load(['/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles/1h5WholeFiles/ISO3_ff' num2str(ff)   '.mat']);
        
  
    vtkwrite(['C1_ff' num2str(ff) '.vtk'], 'structured_points', 'C1', uint8(C1(round(size(C1,1)/2)-round(size(C1,1)/4):round(size(C1,1)/2)+round(size(C1,1)/4),round(size(C1,2)/2)-round(size(C1,2)/4):round(size(C1,2)/2)+round(size(C1,2)/4),round(size(C1,3)/2)-round(size(C1,3)/4):round(size(C1,3)/2)+round(size(C1,3)/4))));
    vtkwrite(['C2_ff' num2str(ff) '.vtk'], 'structured_points', 'C2', uint8(C2(round(size(C1,1)/2)-round(size(C1,1)/4):round(size(C1,1)/2)+round(size(C1,1)/4),round(size(C1,2)/2)-round(size(C1,2)/4):round(size(C1,2)/2)+round(size(C1,2)/4),round(size(C1,3)/2)-round(size(C1,3)/4):round(size(C1,3)/2)+round(size(C1,3)/4))));
    vtkwrite(['C3_ff' num2str(ff) '.vtk'], 'structured_points', 'C3', uint8(C3(round(size(C3,1)/2)-round(size(C3,1)/4):round(size(C3,1)/2)+round(size(C3,1)/4),round(size(C3,2)/2)-round(size(C3,2)/4):round(size(C3,2)/2)+round(size(C3,2)/4),round(size(C3,3)/2)-round(size(C3,3)/4):round(size(C3,3)/2)+round(size(C3,3)/4))));
   
    clear C1 C2 C3
end

%% BEEN IN ITKSNAP READ IN ANNOTATION AND SAVE AS TIFF

headerInfo = nhdr_nrrd_read('/Volumes/GoogleDrive-105399969456537266325/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles/ff1_anno.nrrd', 1);

AN = uint8(headerInfo.data);

ff=1
bfsave(uint8(AN),['ANNO_ff' num2str(ff) '.tiff'])


%% Been in ILASTIK - now look at output

path = '/Volumes/GoogleDrive-105399969456537266325/My Drive/Piljabi/QuantifyLoopInWM/WholeRawFiles/1h5WholeFiles';%'/Volumes/GoogleDrive-105399969456537266325/My Drive/Piljabi/QuantifyLoopInWM/CroppedFiles';
probfilename ='ff=3_DownScaled_Probabilities.h5';% 'ff=1_probabilities.h5';
INFO = h5info([path filesep probfilename])


cc=1;

start = [1 1 1 cc];
count = [Inf Inf Inf 1];
P = h5read([path filesep probfilename],'/probabilities',start, count);
class(P)
%%
size(P)

implayS(P,1)
%%

N1 = imclose(C1,strel('sphere',2))>imgaussfilt3(C1,[5 5 5]);
Ot1 = otsu3D(P(:,:,:,2),.2);

%implaypair(C1.*5,N1)

%implaypair(N,otsu3D(P(:,:,:,2),.2))

N1filt = N1 & Ot1;
implaypair(C1.*5,bwperim(N1filt))
%%

N2 = C2>imgaussfilt3(C2,[5 5 5]);
Ot2 = otsu3D(P(:,:,:,1),.5);

%implaypair(C2.*5,N2)

implaypair(N2,otsu3D(P(:,:,:,2),.2))
%%

Ot3 = otsu3D(P(:,:,:,3),.5);


N1filtBAO = N1filt;
N1filtBAO(Ot3==1)=0;

N1filtBAO = bwareaopen(N1filtBAO,30);
%%
implaypair(Ot3,N1filtBAO)

%%

N2filt = N2 & Ot2;

%%

implaypair(N1filtBAO,C1.*5)


%%
L = uint8(Ot2);
L(N1filtBAO)=2;

implaypair(C1,L)

%%

skel = bwskel(L==2);
%%
implaypair(C1,skel)


%%
vtkwrite('ff1Labels.vtk', 'structured_points', 'ff1Labels', uint8(L));
%%
vtkwrite('C1.vtk', 'structured_points', 'C1', uint8(C1));
vtkwrite('C2.vtk', 'structured_points', 'C2', uint8(C2));

%%
vtkwrite('skelff1.vtk', 'structured_points', 'skel', uint8(skel));

%%

P1 = P(:,:,:,2);
P1 = uint8(255.*P1./max(P1(:)));

implaypair(P1,C1)
%%

vtkwrite('P1.vtk', 'structured_points', 'P1', uint8(P1));

%%
volumeViewer(C1)
%%
volumeViewer(N1filt)

%%
yscale = 1;
xscale = filelist(ff).voxelSizeMicron(2)/filelist(ff).voxelSizeMicron(1);
zscale = filelist(ff).voxelSizeMicron(3)/filelist(ff).voxelSizeMicron(1);


 C1 = imresize3(RAW(:,:,:,1),[round(size(RAW,1).*yscale) round(size(RAW,2).*xscale) round(size(RAW,3).*zscale)]);
 C2 = imresize3(RAW(:,:,:,2),[round(size(RAW,1).*yscale) round(size(RAW,2).*xscale) round(size(RAW,3).*zscale)]);
 C3 = imresize3(RAW(:,:,:,3),[round(size(RAW,1).*yscale) round(size(RAW,2).*xscale) round(size(RAW,3).*zscale)]);
%%
tic
 C1m = medfilt3(C1,[3 3 3]);
 C2m = medfilt3(C2,[3 3 3]);
 C3m = medfilt3(C3,[3 3 3]);
toc
 %%

implaypair(flipxz(C1m),flipxz(C2m))
implaypair((C1m),(C2m))

%% EQUALIZE TEST

C2eq = double(C2m(:,:,50:100))-double(imgaussfilt3(C2m(:,:,50:100),15));


implayS(C2eq,1)
%%
figure
histogram(C2eq(:))


%%
OT1 = imgaussfilt3(double(otsu3D(C1m,0.1)),1)>0.5;
BW1 = double(C1m>imgaussfilt3(C1m,3));
OT2 = imgaussfilt3(double(otsu3D(C2m,0.1)),1)>0.5;
BW2 = double(C2m>imgaussfilt3(C2m,3));
OT3 = imgaussfilt3(double(otsu3D(C3m,0.1)),1)>0.5;
BW3 = double(C3m>imgaussfilt3(C3m,3));

BW1(OT1==0)=0;
BW2(OT2==0)=0;
BW3(OT3==0)=0;

%%

implaypair(BW1,C1m.*5)
implaypair(BW2,C2m.*5)
implaypair(BW3,C3m.*5)

%%
EPME = imgaussfilt3(double(imclose(OT2+imdilate(OT1,strel('sphere',7)),strel('sphere',7))),2)>0.5;

%%
implaypair(EPME(:,:,50:150),C1m(:,:,50:150).*5)

%%
G = otsu3D(imgradient3(EPME),0.5);
%%
G1 = otsu3D(imgradient3(C1m),0.2);
%%
G3 = otsu3D(imgradient3(C3m),0.2);
%%
implaypair(flipxz(bwperim(BW1(:,:,50:150))),flipxz(G3(:,:,50:150)))
%%
implaypair(G(:,:,50:150),C2m(:,:,50:150).*5)
%%
BW = bwperim(BW1) | BW2 | bwperim(BW3) | G;

implaypair(BW(:,:,50:150),C2m(:,:,50:150).*5)


%%
D = -1.*bwdist(BW);
W = watershed(D,6);

implaypair(C1m(:,:,50:250).*5,W(:,:,50:250)==0)


%%
implaypair(flipxz(W26(:,:,50:250)==0),flipxz(W(:,:,50:250)==0))
%%
implaypair((W26(:,:,50:250)==0),(C1m(:,:,50:250).*5))

%% SAVE ITK FILES FOR HAND ANNOTATION IN ITK SNAP
tic
vtkwrite('C1.vtk', 'structured_points', 'C1', C1m(:,:,50:250));
vtkwrite('C2.vtk', 'structured_points', 'C2', C2m(:,:,50:250));
vtkwrite('C3.vtk', 'structured_points', 'C3', C3m(:,:,50:250));
%%
vtkwrite('watershedlines.vtk', 'structured_points', 'W', W(:,:,50:250)==0);
%%
toc
%%
%IM = uint8(zeros([size(C1m(:,:,50:250)) 3]));
IM = uint8(zeros([size(C1m) 3]));

%%
% IM(:,:,:,1) = C1m(:,:,50:250);
% IM(:,:,:,2) = C2m(:,:,50:250);
% IM(:,:,:,3) = C3m(:,:,50:250);
%%
IM(:,:,:,1) = C1m;
IM(:,:,:,2) = C2m;
IM(:,:,:,3) = C3m;


%%
h5create('WHtestBIGcropped.h5','/PiasData',[size(IM,1) size(IM,2) size(IM,3) size(IM,4)],'ChunkSize',[10 10 10 3])
h5write('WHtestBIGcropped.h5','/PiasData',IM)

%%
headerInfo = nhdr_nrrd_read('/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/LazyAnnotation.nrrd', 1);

AN = uint8(headerInfo.data);


numClasses = max(AN(:))

implaypair(C2m(:,:,50:250),AN)

%% extend handdrawn annotation to edge of supervoxel
ANex = AN;

for nn = 1:numClasses
    nn
    TEMP1 = AN==nn;% isolate annotations for just one class
    TEMP2 = W(:,:,50:250);
    TEMP2(TEMP1==0) = 0;
    superVoxelIndW = nonzeros(unique(TEMP2));
    for ii = 1:length(superVoxelIndW)
        ANex(W(:,:,50:250)==superVoxelIndW(ii))=nn;
    end
end


%implaypair(imdilate(OT1(:,:,50:250),strel('sphere',5)),ANex==3)
%%
implaypair(255.*flipxz(uint8(AN==1) + uint8(ANex==1)),flipxz(255.*uint8(AN==1)+C2m(:,:,50:250)))

%%

%vtkwrite('annotation1.vtk', 'structured_points', 'ANex', ANex);
vtkwrite('annotation2.vtk', 'structured_points', 'ANex', ANex);

%%
%bfsave(ANex,'annotation1.ome.tiff')
bfsave(ANex,'annotation2.ome.tiff')
%%

path = "/Users/xqk313/Google Drive/Piljabi/QuantifyLoopInWM";


prob = h5read([path + filesep + "probabilitiesff2.h5"],"/exported_data");

size(prob)
%%
save('Workspaceff2.mat','-v7.3')
%%


implaypair(prob(:,:,:,1),prob(:,:,:,2))
%%
implaypair(prob(:,:,:,1),C2m)
%%
implaypair(flipxz(prob(:,:,:,2)),flipxz(C1m(:,:,50:250)))
%%
implaypair(flipxz(prob(:,:,:,3)),flipxz(C1m(:,:,50:250)))

%%

volumeViewer(C1m(:,:,50:250))

volumeViewer(prob(:,:,:,3)+prob(:,:,:,2))

%%
% %%
% sigma = 1;
% alpha = 0.1;
% beta = 1;
% 
% 
% tic
% LLR = locLap3(imgaussfilt3(C2(:,:,1:5)),sigma,alpha,beta);
% toc
% 
% implaypair(C2(:,:,1:5),LLR)
% 
% %%
% 
% LC = zeros(size(R_cor));
% 
% for zz=1:size(LC,3)
%     zz
%     LC(:,:,zz) = localcontrast(R_cor(:,:,zz),1);
% end
% 
% implaypair(R_cor,LC)
% 
% %%
% 
% 
% 
% %%
% G = imgradient3(imgaussfilt3(R_cor,2));
% 
% implaypair(R_cor,G)
% %%
% 
% AR1 = double(imgaussfilt3(LC,.1)>imgaussfilt3(LC,0.8));
% AR2 = imgaussfilt3(AR1,0.7)>0.5;
% AR3 = imclose(AR2,strel('sphere',1));
% AR3(AR1==0)=0;
% AR3 = imclose(AR3,strel('sphere',1));
% %%
% implaypair(AR3,R_cor)
% %%
% AG = imgaussfilt3(double(imgaussfilt3(G_cor,.5)>imgaussfilt3(G_cor,4)),1)>0.5;
% %AG = imclose(AG,strel('sphere',2));
% %%
% 
% D = zeros(size(AR));
% 
% for zz=1:size(AR,3)
%    D(:,:,zz) = imregionalmin( -1.*bwdist(AR3(:,:,zz))+1); 
% end
% 
% D = imclose(D,strel('sphere',2));
% %implaypair(flipxz(R_cor.*5),flipxz(D))
% 
% implaypair(AR3,(D))
% 
% 
% %%
% BOTH = AG & AR;
% 
% CL  = imclose(BOTH,strel('disk',10));
% INS = imopen(CL-BOTH,strel('sphere',2));
% 
% implaypair(INS,BOTH);
% 
% %%
% im = uint16(G_cor);
% 
% LM = imgaussfilt3(im,0.001);
% LM = imregionalmin(LM);
% LM1 = imopen(LM,strel('sphere',2));
% LM2 = bwareaopen(LM1,5000);
% LM = LM-LM2;
% LM = imclose(LM,strel('sphere',2));
% 
% implaypair(5.*im,LM);
% %%
% 
% W2D = zeros(size(AR));
% 
% for zz=1:size(AR,3)
% %W(:,:,zz) = watershed(imimposemin(imgaussfilt3(G_cor(:,:,zz),1),D(:,:,zz)),6);
%     W2D(:,:,zz) = watershed(imgaussfilt3(R_cor(:,:,zz),1),6);
%     %W2D(:,:,zz) = imerode(W2D(:,:,zz)>0,strel('disk',2));
% end
% 
% %%
% W2Dxz = flipxz(zeros(size(AR)));
% R_corxz = flipxz(R_cor);
% for zz=1:size(W2Dxz,3)
%     zz
%     W2Dxz(:,:,zz) = watershed(imgaussfilt3(R_corxz(:,:,zz),1),6);
% end
% 
% 
% implaypair(W2Dxz==0,R_corxz)
% %%
% implaypair(permute(W2Dxz==0,[3 2 1]),W2D)
% %%
% W2Dyz = flipyz(zeros(size(AR)));
% R_coryz = flipyz(R_cor);
% for zz=1:size(W2Dyz,3)
%     W2Dyz(:,:,zz) = watershed(imgaussfilt3(R_coryz(:,:,zz),1),6);
% end
% 
% 
% %%
% implaypair(flipxz(5.*R_cor),flipxz(W2D>0));
% implaypair((5.*R_cor),(W2D>0));
% %%
% W3D = watershed(imgaussfilt3(R_cor,1),6);
% implaypair(flipxz(5.*R_cor),flipxz(W3D==0));
% %%
% implaypair(flipxz(W2D==0),flipxz(W3D==0));
% 
% %%
% implaypair(flipxz(5.*R_cor),flipxz(W==0));
% %implaypair((GRID),(W==0));
% 
% 
% 
