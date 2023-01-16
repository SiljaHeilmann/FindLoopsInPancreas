

filename ='LI 2015-07-12 MIPGFP_Muc1_pos2.ims';%'P42 D13 pdx ezrin bcatenin-gfp dapi.lif';%% %'LI 2015-07-12 MIPGFP_Muc1_40x_2015_07_14__pos1.ims';%
path ='/Volumes/danstem/Semb/Silja/From Pia/Muc1mCherry X MIPGFP/LI 2015-07-12 MIPGFP_Muc1/';%'/Volumes/DANSTEM/Semb/Christy/Confocal/SP8 Danstem/PP diff betalog';%

read = bfGetReader([path '/' filename]);

omeMeta = read.getMetadataStore();
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % number of time points
stackSizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of Channels

totalNumberOfImages = read.getImageCount()

display(['stackSizeX is ' num2str(stackSizeX) ', stackSizeY is ' num2str(stackSizeY) ', stackSizeZ is ' ...
    num2str(stackSizeZ) ', stackSizeT is ' num2str(stackSizeT) ', stackSizeC is ' num2str(stackSizeC)])

%% IMPORT PART OF THE STACK
clear IM1 IM2 IM3 IM_col_MIP C1ent C1med5

TT_start = 5;
delta_TT = 1;%stackSizeT; 

IM1 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ,delta_TT+1));
IM2 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ,delta_TT+1));
IM3 = uint8(zeros(stackSizeY,stackSizeX,stackSizeZ,delta_TT+1));

IM_col_MIP = uint8(zeros(stackSizeY,stackSizeX,3,delta_TT+1));


% READ IN DATA
count_tt = 0;    

tic
for tt = TT_start:TT_start+delta_TT-1
    tt

    count_tt = count_tt+1
    
    for zz=1:stackSizeZ-1
        iPlane = read.getIndex(zz - 1, 1 - 1, tt - 1) + 1;
        IM1(:,:,zz,count_tt) = bfGetPlane(read,iPlane);
        
         iPlane = read.getIndex(zz - 1, 2 - 1, tt - 1) + 1;
         IM2(:,:,zz,count_tt) = bfGetPlane(read,iPlane);
         
         iPlane = read.getIndex(zz - 1, 3 - 1, tt - 1) + 1;
         IM3(:,:,zz,count_tt) = bfGetPlane(read,iPlane);
         
    end
end
toc

%%

IM_col_MIP(:,:,1,:) = max(IM2,[],3);
IM_col_MIP(:,:,2,:) = max(IM1,[],3);
IM_col_MIP(:,:,3,:) = max(IM3,[],3)-max(IM3,[],3);


implayS(IM_col_MIP,1)
%%
implayS(medfilt3(IM2(:,:,:,5)),2)
%%
volumeViewer(IM2(:,:,:,5))

%%
close all
test=medfilt3(IM2(:,:,:,5));
imhist(test(:))
set(gca,'yscale','log')

TT = otsuthresh(test(:)).*max(test(:))

%%
close all
figure
imshow(IM2(:,:,20,5)>TT.*0.3)

%%

stack1 = medfilt3(IM1(:,:,:,100),[3 3 3]);

stack2 = entropyfilt(stack1);

implayS(stack1,1)

%%
implayS(otsu3D(stack2,0.5),2)

%%
tt= randi(delta_TT+1);

implayS(permute(IM1(:,:,:,tt),[1 2 4 3]),1) 


