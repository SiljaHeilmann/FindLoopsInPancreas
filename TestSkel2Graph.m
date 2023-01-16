% %% test skel to graph
% MAKE FAKE DATA TO TEST
% 
% IM = uint8(rand(40,40,8).*255)>246;
% 
% IM = imclose(IM,strel('sphere',3));
% 
% IM = imfill(bwareaopen(IM,1000,6),'holes');
% 
% IM = bwareaopen(IM,10);
% 
% IMskel = bwskel(IM);
% IMskelDil = imdilate(IMskel,strel('sphere',1));
% 
% %implaypair(IM,IMskelDil)
% 
% 
% 
% volumeViewer(IMskelDil)
%% IMPORT H5 prediction file (here only part of H5 file)

ff = 7;
path = '/Users/xqk313/Documents/Piljabi-LOCAL/H5files-org-prob-seg';

filename = ['ff=' num2str(ff) '_DownScaled_Simple Segmentation.h5'];

INFO = h5info([path filesep filename])

filename2 = ['ff=' num2str(ff) '_DownScaled.h5'];
INFO2 = h5info([path filesep filename2])

%%
imsizeP = INFO.Datasets.Dataspace.Size% x y z c
imsizeIM = INFO2.Datasets.Dataspace.Size% x y z c

%start = [round(imsize(1)/3) 20 round(imsize(3)/3) 1]
start = [1 1 1 1];

countP = imsizeP;%[500 500 100 3];
countIM = imsizeIM;%[500 500 100 3];

P = h5read([path filesep filename],['/' INFO.Datasets.Name],start,countP);
IM = h5read([path filesep filename2],['/' INFO2.Datasets.Name],start,countIM);

%% mesen and ins?
implaypair(P==3,P==4)

%% lumen and epi
implaypair(P==1,P==2)

%%

implaypair(IM(:,:,:,2).*5,IM(:,:,:,1))


%%
implaypair(uint8(bwperim(BW)).*255 +IM(:,:,:,1).*3,uint8(bwperim(EPI)).*255 +IM(:,:,:,2).*1.5)
%% MESEN (and nothing) MASK

MESENorg = P==3;

MESEN = imclose(MESENorg,strel('sphere',4));
%MESEN = bwmorph3(MESEN,'majority');

implaypair(IM(:,:,:,2).*5,MESEN)
%% PROCESS LUMEN Simple segmentation OUTPUT TO MAKE BETTER SEGMENTATION 

BWorg = P==1;%otsu3D(IM,1);

%BW = bwareaopen(BWorg,400);
BW = imclose(BWorg,strel('sphere',1));
BW = imfill(BW,"holes");
BW = bwmorph3(BW,'clean');
BW = bwmorph3(BW,'fill');
BW = bwmorph3(BW,'majority');
BW = imfill(BW,"holes");

tic
BWconnMask = imclose(BW,strel('cube',20));
toc

BWconnMask = bwareaopen(BWconnMask,10000);

implaypair(BWconnMask,BW)

BW(BWconnMask==0)=0;% use BWconnMask to get rid of small objects far aways from network (we will still need more clean up... in itksnap?)

implaypair(flipxz(BWorg),flipxz(BW))

%% SKELETONIZE EPI SEGMENTATION 


EPIorg = P==2 | P==4 | BW==1;%otsu3D(IM,1);

%EPI = bwareaopen(EPIorg,400);
EPI = imopen(EPIorg,strel('cube',2));
EPI = imfill(EPI,"holes");
EPI = bwmorph3(EPI,'clean');
EPI = bwmorph3(EPI,'fill');
EPI = bwmorph3(EPI,'majority');
%EPI = imfill(EPI,"holes");

EPIconnMask = imclose(EPI,strel('cube',15));

EPIconnMask = bwareaopen(EPIconnMask,200000.*10); % EPI should be really large - thereshold for getting rid of small object is more than 5 times that of lumen mask 

EPI(EPIconnMask==0)=0;

EPIskel = bwskel(EPI);
%EPIskel = imclose(EPIskel,strel('sphere',3));

EPIskelDil = imdilate(EPIskel,strel('cube',1));

% EPIskel = bwskel(EPIskelDil,'MinBranchLength',1);
% EPIskelDil = imdilate(EPIskel,strel('sphere',1));

EPIskelDil(EPI==0) = 0;
EPIskel(EPI==0) = 0;

implaypair(EPI,EPIorg)
implaypair(EPI,MESEN)

%volumeViewer([EPIskelDil, EPIorg])

volumeViewer(uint8(EPIskelDil).*3+uint8(EPI))



%%
IMskel = bwskel(BW);
IMskel = imclose(IMskel,strel('sphere',3));

IMskelDil = imdilate(IMskel,strel('sphere',1));

% IMskel = bwskel(IMskelDil,'MinBranchLength',1);
% IMskelDil = imdilate(IMskel,strel('sphere',1));

IMskelDil(BW==0) = 0;
IMskel(BW==0) = 0;

%implaypair(IMskelDil,BWorg)


%volumeViewer([IMskelDil, BWorg])
%% MAKE LABELLED IMAGE FOR CLEAN UP IN ITK SNAP

LBW = bwlabeln(BW,26);
max(LBW(:));
LBWshow = uint8(relabelForVolView(LBW));


%implaypair(LBWshow>0,LBW>0)

LBWshow = LBWshow + uint8(LBWshow>1);% increase all labels by 1 (thus freeing up label 1 for EPI) 

TEMP = uint8(EPI);% MAKE EPI MASK HAVE LABEL 1 (EVERYWHERE THATS NOT LUMEN...)
TEMP(LBWshow>1)=0;

LBWshow = TEMP+LBWshow;
LBWshow = uint8(relabelForVolView(LBWshow));

%volumeViewer(uint8(IMskelDil).*3+uint8(BW))

vtkwrite(['Labels_ff' num2str(ff) '.vtk'], 'structured_points', 'LBWshow', LBWshow);
%vtkwrite(['Labels_ff' num2str(ff) '2ndRound.vtk'], 'structured_points', 'LBWshow', LBWshow);


 vtkwrite(['C1_ff' num2str(ff) '.vtk'], 'structured_points', 'C1', IM(:,:,:,1));
 vtkwrite(['C2_ff' num2str(ff) '.vtk'], 'structured_points', 'C2', IM(:,:,:,2));

%% Been in ITK snap... connected lumen network and broke connections to EPI that needs to be removed

headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/Labels_ff' num2str(ff) '.nrrd'], 1);
%headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/Labels_ff' num2str(ff) '2ndRound.nrrd'], 1);

AN = uint8(headerInfo.data);


EPI = AN==1;% annotatation with 1 is EPITHELIUM

%EPI(P==4)=1;% I forgot ins cells for ff = 1 DOH!

%BW = AN==2;


openingNumPixels = 4;

EPI = imopen(EPI,strel('cube',openingNumPixels));% they are &#ˆ%$##&* still connected!!! grr


LABELS = bwlabeln(EPI,6);

% find largest conn comp % that should be the pancreas

r = regionprops3(LABELS,'Volume');
[sorted,index] = sort(r.Volume);
biggestIndex = index(end);

LABELSshow = relabelForVolView(LABELS);

volumeViewer(LABELS==biggestIndex)
%%
EPI = LABELS==biggestIndex | AN>1;

EPI = returnLargestConnComp(EPI,26);

% EPI = imdilate(BW,strel('sphere',1)) | EPI;
% 
 EPI = imfill(EPI,"holes");
 EPI = bwmorph3(EPI,'clean');
 EPI = bwmorph3(EPI,'fill');
 EPI = bwmorph3(EPI,'majority');

implaypair(EPI,AN>1)

implaypair(EPI,LUMENbiggest)

%%
%EPI = imdilate(BW,strel('sphere',1)) | EPI;

implaypair(BW,EPI)
%volumeViewer(EPI)

%%
vtkwrite(['Labels_ff' num2str(ff) 'Cleaned1st.vtk'], 'structured_points', 'JustEpi', LABELS==biggestIndex);
%%
implaypair(LABELS,LABELS==biggestIndex)
%%

% EPI MASK AFTER CORRECTION!

EPImask = LABELS==biggestIndex;
EPImask = imclose(EPImask,strel('cube',openingNumPixels)); % to compensate for the imopen earlier/before...

%EPI = AN==1;% Start with fresh annotation
EPI(EPImask==0)=0;% remove things outside biggest conn comp

tic
% fill and close to overlap with lumen inside
EPIcl = imclose(EPI,strel('cube',3));
EPIcl = imfill(EPIcl,'holes');
toc

% find conn comp in LUMEN that overlap with EPI mask
LUMEN = bwlabeln(AN>1); % start with all conn comp in annotation
%LUMEN = bwlabeln(BW>0); % If youve been twice in itk snap it might be BW...

TEMP = LUMEN;
TEMP(EPIcl==0)=0;

indexToKeep = unique(nonzeros(TEMP(:)));% index that overlap with imclosed EPI biggest conn comp (the ones we want to keep!)

tic
LUMEN = ismember(LUMEN,indexToKeep);
toc

LABELLUMEN = bwlabeln(LUMEN,26);

howManyConnComp = max(LABELLUMEN(:))

LABELLUMENshow = relabelForVolView(LABELLUMEN);

LUMENbiggest = returnLargestConnComp(LUMEN,26);

volumeViewer(LUMENbiggest)


implaypair(LUMENbiggest,AN>1)
%% Use lumen segmentation as a mask to filter regional thresh image (only for ff 4 and 5 ?)

BWg = imgaussfilt3(IM(:,:,:,1),[15,15,15]);

BWmedorg = IM(:,:,:,1)>BWg;
%BWmedorg(LUMEN==0)=0;
BWmedorg(LUMENbiggest==0)=0;

BWmedorg(AN==2)=1;

BWmed = imclose(BWmedorg,strel('sphere',2));
BWmed = imfill(BWmed,"holes");
BWmed = bwmorph3(BWmed,'clean');
BWmed = bwmorph3(BWmed,'fill');
BWmed = bwmorph3(BWmed,'majority');
BWmed = imfill(BWmed,"holes");


imForShow = IM(:,:,:,1).*5;
imForShow(bwperim(LUMENbiggest)==1)=255;
%imForShow(bwperim(BW)==1)=255;

implaypair(flipxz(imForShow),flipxz(bwperim(BWmed)))
%%
implaypair((imForShow),(bwperim(BWmed)))
%%
implaypair((imForShow),(bwperim(BWmed)))

% check if you can get rid of small conn comp
LBWmed = bwlabeln(BWmed,26);
max(LBWmed(:))
LBWmedSHOW = relabelForVolView(LBWmed);
volumeViewer(BWmed)

LUMEN = returnLargestConnComp(BWmed,26);
%%  EPI mask NEEDS TO BE BOTH EPI AND LUMEN


%EPI = imdilate(LUMEN,strel('sphere',1)) | EPIcl;
EPI = imdilate(LUMEN,strel('sphere',1)) | EPI;
%%
implaypair(P==4,EPI)
%%
vtkwrite(['Lumen_ff' num2str(ff) '.vtk'], 'structured_points', 'JustLUMEN', BWmed);


%%


implaypair(flipxz(IM(:,:,:,1)),flipxz(bwperim(LUMEN)))
%%

LUMENbiggest = LUMEN;



%%
implaypair(LUMENbiggest,EPI)

%% SKELETONIZE LUMEN SEGMENTATION (PERFORM A FEW TRICKS TO MINIMIZE 'FAKE' LOOPS)

BW = LUMENbiggest;

BW = imfill(BW,"holes");
BW = bwmorph3(BW,'clean');
BW = bwmorph3(BW,'fill');
BW = bwmorph3(BW,'majority');
BW = imfill(BW,"holes");

IMskel = bwskel(BW,'MinBranchLength',15);
%IMskel = imclose(IMskel,strel('sphere',3));
%IMskel = bwskel(IMskel,'MinBranchLength',10);
%IMskel(BW==0)=0;

IMskelDil = imdilate(IMskel,strel('sphere',1));

IMskelDil(BW==0) = 0;
%IMskel(BW==0) = 0;
%%
implaypair(P==4,EPI)
%% SAVE LUMEN AND EPI FOR RECORD 

vtkwrite(['LUMEN_FINAL_ff' num2str(ff) '.vtk'], 'structured_points', 'LUMEN', BW);
vtkwrite(['IMskelDilFINAL_ff' num2str(ff) '.vtk'], 'structured_points', 'IMskelDil', IMskelDil);
vtkwrite(['EPI_FINAL_ff' num2str(ff) '.vtk'], 'structured_points', 'EPI', EPI);
%% IF YOU HAVE TO START FROM RECORD - go in itk snap and make nrrd files and go from here...

headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/LUMEN_FINAL_ff' num2str(ff) '.nrrd'], 1);
LUMEN = uint8(headerInfo.data);

headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/EPI_FINAL_ff' num2str(ff) '.nrrd'], 1);
EPI = uint8(headerInfo.data);


implaypair(LUMEN,EPI)
%%

implaypair(BW,EPI)
%%
% 
% % check if you can get rid of small conn comp
% LBW = bwlabeln(BW,26);
% max(LBW(:))
% LBWSHOW = relabelForVolView(LBW);
% volumeViewer(BW)
% LUMEN = BW;
% LUMENbiggest = returnLargestConnComp(BW,26);
%% ONLY FOR IMAGES WITH INS channel (not 4 and 5)

LUMENplusINS = LUMEN;
INS = P==4;
INS(EPI==0)=0;
INS = imfill(INS,"holes");
INS = bwmorph3(INS,'clean');
INS = bwmorph3(INS,'fill');
INS = bwmorph3(INS,'majority');
LUMENplusINS(INS==1)=2;
 %%
implaypair(EPI,LUMENplusINS)
%%
vtkwrite(['LUMENplusINS_FINAL_ff' num2str(ff) '.vtk'], 'structured_points', 'LUMENplusINS', LUMENplusINS);

%%

% headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/LUMEN_FINAL_ff' num2str(ff) '.nrrd'], 1);
% BW2 = uint8(headerInfo.data);
% IMskel = bwskel(logical(BW2),'MinBranchLength',15);
% IMskelDil = imdilate(IMskel,strel('sphere',1));

LIMskel = bwlabeln(IMskel,26); % how many conn comp in skeleton?
NUMccompINskel = max(LIMskel(:))
%implaypair(BW,IMskel)
%%
IMskel = LIMskel==1;
%%
IMskelDil = imdilate(IMskel,strel('sphere',1));
%%
implaypair(LIMskel==1,LIMskel==2)

%% 
close all
figure
imshowpair(max(LIMskel==1,[],3),max(LIMskel>1,[],3))
%%

r = regionprops3(LIMskel,'Volume');

[sorted, index]=sort(r.Volume)

index(end)
%%
max(LIMskel(:))
implaypair(LIMskel>1,LIMskel==1)
implaypair(BW,LIMskel==1)

%%

LABELskel = bwlabeln(IMskelDil,26);
LABELskelSHOW = relabelForVolView(LABELskel);

volumeViewer([uint8(IMskelDil)])
%% IDENTIFY LIST OF VOXELS IN SKELETON AND THEIR VOXEL ID

tic
G_BIG = skel2graph(IMskel);% approx 2 min for whole pancreas
toc

%% how many connected components are in graph? Hopefully just one...
bins = conncomp(G_BIG);

unique(bins)
%% FIND END POINTS AND BRANCHPOINTS AND LABEL THEM
tic
IMep = bwmorph3(IMskel,'endpoints');
IMbp = bwmorph3(IMskel,'branchpoints');
Lep = bwlabeln(IMep,6); % label endpoints
Lbp = bwlabeln(IMbp,26); % a lot of bp sits in clusters and should not be counted individually (all endpoints sits free)
tic

%% FIND BRANCH MIDPOINTS - first get geodesic distance map inside IMskelDil from branch points - them locate its regional max (these are found on the mid point between branch points)
tic
DfromBP = bwdistgeodesic(IMskelDil==1,Lbp>0);
toc

DfromBPNoNANs = DfromBP;
DfromBPNoNANs(isnan(DfromBP)) = min(DfromBP(:));

BranchMidPoints = imregionalmax(DfromBPNoNANs,26);
% clean up regional max that overlap with end or branch points
BranchMidPoints(imdilate(Lep>0,strel("sphere",2))==1)=0;%remove regional max that overlap with an end point!
BranchMidPoints(imdilate(Lbp>0,strel("sphere",2))==1)=0;%remove regional max that overlap with a branch point!

%%

implaypair(single(IMskelDil==1)+2.*single(Lbp>0)+3.*single(imdilate(Lep>0,strel("sphere",2))),BranchMidPoints)

%%
implaypair(flipxz(single(IMskelDil==1)+2.*single(Lbp>0)+2.*single(Lep>0)),flipxz(BranchMidPoints))

%%
implaypair(IMskelDil==1,DfromBP+50.*single(Lbp>0))

%%

figure
imshow(max(DfromBP,[],3),[1 20])
%%

BranchMidPointsSHOW = relabelForVolView(bwlabeln(BranchMidPoints));
BranchMidPointsSHOW = imdilate(BranchMidPointsSHOW,strel("sphere",1));

volumeViewer(uint8(IMskelDil));
%%
implaypair(bwdist(imcomplement(logical(BW))),BranchMidPoints)
%% extract max intensity (distance from nearest branch point) in each branch midpoint - this times 2 gives length of branches between BPs

r1 = regionprops3(bwlabeln(BranchMidPoints),DfromBPNoNANs,"MaxIntensity")
r2 = regionprops3(bwlabeln(BranchMidPoints),bwdist(imcomplement(logical(BW))),"MaxIntensity")
%%
close all

TEMP1 = r1.MaxIntensity.*2;
TEMP2 = r2.MaxIntensity;

data1 = TEMP1(isfinite(TEMP1) & isfinite(TEMP2));
data2 = TEMP2(isfinite(TEMP1) & isfinite(TEMP2));

figure
boxchart(discretize(data1,20),data2,'Notch','on')

[rho,p]=corr(data1,data2,'type','Spearman') % slight anti correlation between length of lumen and width at center (NB only extracting lumen width at central point!)

%plot(data1,data2,'o')
xlabel('lumen length')
ylabel('lumen radius at midpoint')
%%

tic
LbpSHOW = relabelForVolView(Lbp);
toc

LbpSHOW = imdilate(LbpSHOW,strel("sphere",1));

volumeViewer(uint8(IMskelDil))




%% REMOVE SMALL CYCLES FROM GRAPH 


% find cycles
tic
[cycles, edgecycles ] = cyclebasis(G_BIG);
toc
size(cycles)


thresh = 15;
tic
Gnew = removeSmallCycles(G_BIG,thresh);
toc
tic
Gnew = removeSmallCycles(Gnew,thresh);% another pass removed 3 more edges?! Then no more...
toc

% find cycles in new graph
tic
[cycles, edgecycles ] = cyclebasis(Gnew);
toc
size(cycles)

%%
% MATLAB cyclebasis tend to find a basis where some loops overlap more than
% they have to. I wish to find the loop basis with minimal overlap (this function does that and returns number of times cycles where rearranged to reduce overlap)
tic
[cyclesNEW,numChanges] = minimizeLOOPoverlap(cycles)% takes a looooooong time (5 hours or more?)
toc

size(cyclesNEW)

% are there any small lops left? If yes filter out - still nessesary!
count=0;
indexToRemove =[];
for ii = 1:length(cyclesNEW)
    if length(cyclesNEW{ii})<=thresh
        count = count + 1;
        indexToRemove(count)=ii;
    end
end

count
cyclesFilt = cyclesNEW;
cyclesFilt(indexToRemove)=[];
%%
size(cyclesFilt)

save(['/Users/xqk313/Documents/DocumentsSUN1008821/Piljabi-LOCAL/PackagesWM/cyclesFilt_ff' num2str(ff) 'extra.mat'],'cyclesFilt');

%%
save(['/Users/xqk313/Documents/Piljabi-LOCAL/cyclesFilt/cyclesFilt_ff' num2str(ff) '.mat'],'cyclesFilt');

%%
% cycles now contain a list of all midline voxels (nodes) that participate in
% a loop... 

tic
IMloops = labelLoops(IMskel,Gnew,cyclesFilt);
toc


%%

tic
IMloopsclean = cleanUpLabelledIm(IMloops,2);
toc
%%
save(['/Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/IMloopsclean_ff' num2str(ff) '.mat'],'IMloopsclean');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare packages for making loop tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all

ff = 5;

headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/LUMEN_FINAL_ff' num2str(ff) '.nrrd'], 1);
BW = logical(headerInfo.data);

% headerInfo = nhdr_nrrd_read(['/Users/xqk313/Documents/Piljabi-LOCAL/vtkFilesANDnrrdFiles/EPI_FINAL_ff' num2str(ff) '.nrrd'], 1);
% EPI = uint8(headerInfo.data);


implaypair(BW,EPI)
%%

load(['/Users/xqk313/Documents/Piljabi-LOCAL/cyclesFilt/IMloopsclean_ff' num2str(ff) '.mat']); % gets IMloopsclean
load(['/Users/xqk313/Documents/Piljabi-LOCAL/cyclesFilt/cyclesFilt_ff' num2str(ff) '.mat']); % gets cyclesFilt
%%

IMskel = bwskel(logical(BW),'MinBranchLength',15);

IMskelDil = imdilate(IMskel,strel('sphere',1));
IMskelDil(BW==0) = 0;

LIMskel = bwlabeln(IMskel,26); % how many conn comp in skeleton?
NUMccompINskel = max(LIMskel(:))

implaypair(LIMskel==1,IMskel)

%%
tic
G_BIG = skel2graph(LIMskel==1);% approx 2 min for whole pancreas
toc

thresh = 15;
tic
Gnew = removeSmallCycles(G_BIG,thresh);
toc
tic
Gnew = removeSmallCycles(Gnew,thresh);% another pass removed 3 more edges?! Then no more...
toc
%%
%'/Users/xqk313/Documents/Piljabi-LOCAL/Packages'

path = '/Users/xqk313/Documents/DocumentsSUN1008821/Piljabi-LOCAL/PackagesWM';
filename = ['Package_ff' num2str(ff) 'extra.mat'];

save([path filesep filename],'cyclesFilt','IMloopsclean','IMskelDil','IMskel','Gnew','BW','EPI')

%%
ff=7
path = '/Users/xqk313/Documents/Piljabi-LOCAL/Packages';
filename = ['Package_ff' num2str(ff) '.mat'];

load([path filesep filename],'cyclesFilt','IMloopsclean','IMskelDil','IMskel','Gnew','BW','EPI')




%%
IMloopscleanSHOW = relabelForVolView(imdilate(IMloopsclean,strel('sphere',1)));

volumeViewer(uint8(IMskelDil).*3+uint8(BW)) 
%% CALCULATE PRINCIPLES AXIS LENGTHS BEFORE MAKING CONVEX HULL REPRESENTATION

close all
figure
imshow(max(IMloopsclean,[],3),[]); colormap('colorcube')

r = regionprops3(IMloopsclean,"PrincipalAxisLength");
%%
edges = logspace(0,3,50);%0:5:300;
close all
nexttile
histogram(r.PrincipalAxisLength(:,1),edges,'Normalization','pdf');hold on
histogram(r.PrincipalAxisLength(:,2),edges,'Normalization','pdf');hold on
ylabel('Fraction of loops')
xlabel('Length of loop principal axis (measured on midline rep.), in microns')
legend('Long axis','Short axis')
set(gca,"XScale",'log');

nexttile
histogram(r.PrincipalAxisLength(:,1)./r.PrincipalAxisLength(:,2),'Normalization','pdf');hold on
ylabel('Fraction of loops')
xlabel('Ratio of long and short principal axis (L/S) (measured on midline rep.)')

nexttile
plot(r.PrincipalAxisLength(:,1),r.PrincipalAxisLength(:,1)./r.PrincipalAxisLength(:,2),'.')
ylabel('Ratio of long and short principal axis (L/S) (measured on midline rep.)')
xlabel('Length of loop long principal axis (measured on midline rep.), in microns')

nexttile
boxchart(discretize(r.PrincipalAxisLength(:,1),5)',(r.PrincipalAxisLength(:,1)./r.PrincipalAxisLength(:,2))','Notch','on')
ylabel('Ratio of long and short principal axis (L/S) (measured on midline rep.)')
xlabel('Length of loop long principal axis (measured on midline rep.), in microns')





%%  MAKE LOOP TABLE!
% index/row number in cyclesFilt, IMloops, IMloopsclean and loopTable is
% the same!!! (IMloops got it from cyclesFilt and so forth...)

close all
tic
loopTable = regionprops3(IMloopsclean,'Volume','ConvexImage','BoundingBox','ConvexVolume','Centroid','PrincipalAxisLength','VoxelIdxList');%
toc

LumenWidthImage = bwdist(imcomplement(BW));
EPIWidthImage   = bwdist(imcomplement(EPI));

IMholeOrganCVH = logical(zeros(size(BW)));

tableHoleOrganCVH = regionprops3(uint8(BW).*2,"ConvexImage","BoundingBox");
BB = tableHoleOrganCVH(2,:).BoundingBox;
IMholeOrganCVH(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1)=tableHoleOrganCVH.ConvexImage{2}; 

DistIMholeOrganCVH = bwdist(imcomplement(IMholeOrganCVH));%distance to convex hull of whole organ

% FIND END POINTS AND BRANCHPOINTS AND LABEL THEM
tic
IMep = bwmorph3(IMskel,'endpoints');
IMbp = bwmorph3(IMskel,'branchpoints');
Lep = bwlabeln(IMep,6); % label endpoints
Lbp = bwlabeln(IMbp,26); % a lot of bp sits in clusters and should not be counted individually (all endpoints sits free)
tic

DfromBP = bwdistgeodesic(IMskelDil==1,Lbp>0);

DfromBPNoNANs = DfromBP;
DfromBPNoNANs(isnan(DfromBP)) = min(DfromBP(:));

BranchMidPoints = imregionalmax(DfromBPNoNANs,26);
BranchMidPoints(imdilate(Lep>0,strel("sphere",2))==1)=0;%remove regional max that overlap with an end point!
BranchMidPoints(imdilate(Lbp>0,strel("sphere",2))==1)=0;%remove regional max that overlap with a branch point!

LBranchMidPoints = bwlabeln(BranchMidPoints,26);

tic
loopTable.DonutHoleImage     = cell(size(loopTable.ConvexImage));
loopTable.DonutHoleVolume    = NaN(size(loopTable.ConvexImage));
loopTable.DonutHolePrinAxis  = NaN([size(loopTable.ConvexImage,1) 3]);
loopTable.DonutHoleCentroid  = NaN([size(loopTable.ConvexImage,1) 3]);
loopTable.NumMesenGoThru     = NaN(size(loopTable.ConvexImage));
loopTable.LiniarIndexInImage = cell(size(loopTable.ConvexImage));
loopTable.BranchPointsIndex  = cell(size(loopTable.ConvexImage));
loopTable.NumBranchPoints    = NaN(size(loopTable.ConvexImage));
loopTable.LumenRadiusVoxelValues = cell(size(loopTable.ConvexImage));
loopTable.OrganCVHVoxelValues    = cell(size(loopTable.ConvexImage));
loopTable.EPIRadiusVoxelValues   = cell(size(loopTable.ConvexImage));

% Empty images for convex hull labelled image and mask (Mask showing where regions overlap with 2, 3 ect )
IMloopsCVH      = uint16(zeros(size(IMloopsclean)));
IMloopsCVHmask  = uint16(zeros(size(IMloopsclean)));
IMdonutHole     = uint16(zeros(size(IMloopsclean)));
IMdonutHoleMask = uint16(zeros(size(IMloopsclean)));

TEMP = uint16(zeros(size(IMloopsclean)));

tic
for ii=1:size(loopTable,1)
    ii
    TEMP(TEMP>0)=0;
    currentConvexImage = uint16(cell2mat(loopTable(ii,:).ConvexImage));
    BB = ceil(loopTable(ii,:).BoundingBox);
    if isempty(currentConvexImage)==0 & sum(currentConvexImage(:))>10% filter out tiny loops
     %   nexttile
     %   imshow(sum(currentConvexImages,3),[])
        TEMP(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1) = currentConvexImage;

        currentConvexImageMINUSEPI = currentConvexImage;
        currentConvexImageMINUSEPI(EPI(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1)==0) = 0;% remove mesencyme

%        volumeViewer([currentConvexImageMINUSEPI,bwskel(logical(currentConvexImageMINUSEPI)),currentConvexImage])

        currentSkel = bwskel(logical(currentConvexImageMINUSEPI));

        if sum(currentSkel(:))>2 % more than one pixel left...

            Gcurrent = skel2graph(currentSkel);
            thresh = 15;
            Gcurrent = removeSmallCycles(Gcurrent,thresh);
            Gcurrent = removeSmallCycles(Gcurrent,thresh);

            currentcycles = cyclebasis(Gcurrent);

            numMesenGoThru = size(currentcycles,1);
        else
            numMesenGoThru = NaN;
        end

        currentConvexImageMINUSBW = currentConvexImage;
        currentConvexImageMINUSBW(BW(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1)==1) = 0;
        currentConvexImageMINUSBW = imopen(currentConvexImageMINUSBW,strel('cube',2));% remove very thin structures
        currentConvexImageMINUSBW = cleanUpLabelledIm(currentConvexImageMINUSBW,0); %keep only largest conn comp
        loopTable.DonutHoleImage(ii) = {logical(currentConvexImageMINUSBW)};% add DonutHole Image on to table

        tempR = regionprops3(currentConvexImageMINUSBW.*2,'Volume','PrincipalAxisLength','Centroid');% multiply by two so that every pixel is like region number 2 ()
       if isempty(tempR)==0
        loopTable.DonutHoleVolume(ii) = tempR.Volume(2);
        loopTable.DonutHolePrinAxis(ii,:) = tempR.PrincipalAxisLength(2,:);
        loopTable.DonutHoleCentroid(ii,:) = tempR.Centroid(2,:) + [BB(1) BB(2) BB(3)];
        loopTable.NumMesenGoThru(ii) = numMesenGoThru;
       end
        % volumeViewer([BW(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1),currentConvexImage,currentConvexImageMINUSBW,IMloopsclean(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1)])
        
        % I make these for visualization purposes but due to overlap they
        % should not be used for quantification!!!!
        IMloopsCVHmask = IMloopsCVHmask + TEMP;% will have >1 when regions overlap
        IMloopsCVH = IMloopsCVH + TEMP.*ii;
        
        TEMP(TEMP>0)=0;% reuse TEMP for donut hole labelled image

        TEMP(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1) = uint16(currentConvexImageMINUSBW);
        IMdonutHoleMask = IMdonutHoleMask + TEMP;% will have >1 when regions overlap

        IMdonutHole = IMdonutHole + TEMP.*ii;
    end
    loopTable.LiniarIndexInImage(ii) = {Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii})};
    loopTable.BranchPointsIndex(ii) = {nonzeros(unique(Lbp(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))))};% Index of branch points in loop
    loopTable.NumBranchPoints(ii) = size(cell2mat(loopTable.BranchPointsIndex(ii)),1);% number of branch points in loop

    loopTable.MidBranchPointsIndex(ii) = {nonzeros(unique(LBranchMidPoints(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))))}; % index of midbranch points in loop
    
    loopTable.LumenRadiusVoxelValues(ii) = {LumenWidthImage(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))};
    loopTable.OrganCVHVoxelValues(ii)    = {DistIMholeOrganCVH(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))};

    loopTable.EPIRadiusVoxelValues(ii) = {EPIWidthImage(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))};

    %Lep(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii})) % extract pixel
    %values in loop

end
toc

% NB these two are only used for visualization!
IMloopsCVH(IMloopsCVHmask>1)=0;% set index in overlapping regions to zero (adding will introduce new indexes...)
IMdonutHole(IMdonutHoleMask>1)=0;

loopTable.CycleNodeID = cyclesFilt;

loopTable.DeviationFromFlatCirclePlane = loopTable.DonutHoleVolume ./(pi.*(loopTable.Volume./(2.*pi)).^2); 

loopTable.MeanLumenRadius = cellfun(@mean,loopTable.LumenRadiusVoxelValues);
loopTable.MedianLumenRadius = cellfun(@median,loopTable.LumenRadiusVoxelValues);
loopTable.ModeLumenRadius = cellfun(@mode,loopTable.LumenRadiusVoxelValues);
loopTable.STDLumenRadius = cellfun(@std,loopTable.LumenRadiusVoxelValues);
loopTable.ExcessKurtosisLumenRadius = cellfun(@kurtosis,loopTable.LumenRadiusVoxelValues)-3;
loopTable.LumenRadiusDistSkewness = 3.*(loopTable.MeanLumenRadius - loopTable.MedianLumenRadius)./loopTable.STDLumenRadius;

loopTable.MeanOrganCVHdist = cellfun(@mean,loopTable.OrganCVHVoxelValues);
loopTable.MedianOrganCVHdist = cellfun(@median,loopTable.OrganCVHVoxelValues);
loopTable.ModeOrganCVHdist = cellfun(@mode,loopTable.OrganCVHVoxelValues);
loopTable.STDOrganCVHdist = cellfun(@std,loopTable.OrganCVHVoxelValues);
loopTable.ExcessKurtosisOrganCVHdist = cellfun(@kurtosis,loopTable.OrganCVHVoxelValues)-3;
loopTable.OrganCVHdistSkewness = 3.*(loopTable.MeanOrganCVHdist - loopTable.MedianOrganCVHdist)./loopTable.STDOrganCVHdist;

loopTable.approxRho = loopTable.DonutHolePrinAxis(:,1)./(loopTable.MeanLumenRadius.*2);% ratio between diameter of donut hole and diameter of surrounding lumen (higher thin, smaller fat loop)
loopTable.approxRho_min = loopTable.DonutHolePrinAxis(:,2)./(loopTable.MeanLumenRadius.*2);% I assume here that third prin axis will be perpendicular to the plane of the loop 

save(['/Users/xqk313/Documents/Piljabi-LOCAL/cyclesFilt/loopTable_ff' num2str(ff) '.mat'],'loopTable');

loopTable

%%
load(['/Users/xqk313/Documents/Piljabi-LOCAL/cyclesFilt/loopTable_ff' num2str(ff) '.mat']);


%%
NumLoops = size(unique(nonzeros(IMloopsclean(:))),1)
NumLoopsMesen = size(nonzeros(loopTable.NumMesenGoThru(~isnan(loopTable.NumMesenGoThru))),1)

frac = NumLoopsMesen/NumLoops

%% VISUALIZE LOOP DENSITY?

tic
skelOfLoops = bwskel(imfill(IMdonutHoleMask>0,'holes'));% reduce donut holes to their skeleton (one line or one point per hole)
toc
%%
tic
loopGaussBlur = imgaussfilt3(double(imfill(IMdonutHoleMask>0,'holes')),[25 25 25]);
toc
%%
implaypair(imgaussfilt3(double(INS),[10 10 10]),loopGaussBlur)
%%
implaypair(IMskelDil,loopGaussBlur)

%%

volumeViewer(double(IMskelDil)+loopGaussBlur.*100)

%%
close all
    
nexttile
histogram(loopTable.NumBranchPoints,50)

%%

C = centrality(Gnew,"degree")



close all
    
nexttile
histogram(C)

Gnew.Nodes.Degree = C;
Gnew.Nodes.EndNode = Gnew.Nodes.Degree==1;% will have a one if its an end node
Gnew.Nodes.BranchNode = Gnew.Nodes.Degree>2;% will have a one if its an branch nodes

subNodeIDs = Gnew.Nodes.Index(Gnew.Nodes.BranchNode==1 | Gnew.Nodes.EndNode==1,:);% the nodes we want to use in simplified network bn's and en's

vec1 = single(Gnew.Nodes.Degree(Gnew.Edges.EndNodes(:,1),:)==2);% vec with 1 everywhere where first endnode has degree 2
vec2 = single(Gnew.Nodes.Degree(Gnew.Edges.EndNodes(:,2),:)==2);% vec with 1 everywhere where second endnode has degree 2

edgeWeights = single(vec1+vec2 ~= 2);% 1 if above 2 or below 2, 0 if 2... ( since we want edge weigths to be 0 between two nodes that both have degree 2)

Gnew.Edges.Weight = edgeWeights;

%%
tic
d1 = distances(Gnew,subNodeIDs,subNodeIDs,'Method','positive');% distances along paths of nodes with degree 2 will have zero weigth
d2 = distances(Gnew,subNodeIDs,subNodeIDs,"Method",'unweighted');
toc
%%
close all
nexttile
histogram(d1); hold on
nexttile
histogram(d2); hold on

%%

close all
nexttile
imshow(d1,[])

A = d1==2 | d1==1 | d1==0;

Aupper = triu(A);
Aupper = Aupper - diag(ones(1,size(A,1)));


nexttile
imshow(Aupper,[0 1])

sum(A(:))

Edges = table(zeros(sum(Aupper(:)),2),'VariableNames',{'EndNodes'});

count = 1;
for cc = 1:size(Aupper,1)
    rr = find(Aupper(:,cc));

    if isempty(rr)==0

        Edges.EndNodes(count:count+size(rr,1)-1,:) = [ones(size(rr,1),1).*cc,rr];
        count = count+size(rr,1);
    end
end

Nodes = Gnew.Nodes(subNodeIDs,:)



nexttile
histogram(d1(:))

Gsimple = graph(Edges,Nodes,'omitselfloops')%,string(subNodeIDs)

figure
plot(Gsimple,'XData',Gsimple.Nodes.PixelCoor(:,1),'YData',Gsimple.Nodes.PixelCoor(:,2),'ZData',Gsimple.Nodes.PixelCoor(:,3))

figure
plot(Gsimple)


[cc, sizes] =conncomp(Gsimple);
sort(sizes)

GsimpleCycleBasis = cyclebasis(Gsimple);
%%
loopSizes = cellfun('size',GsimpleCycleBasis,2)

close all
figure
histogram(loopSizes,0:1:20)
%%

IMloopsCVH(IMloopsCVHmask>1)=0;% take out all overlapping regions - it can introduce new indexes!!!!

IMdonutHole(IMdonutHoleMask>1)=0;% take out all overlapping regions - it can introduce new indexes!!!!
toc

IMloopsCVHshow   = relabelForVolView(IMloopsCVH);
IMdonutHoleShow  = relabelForVolView(IMdonutHole);

volumeViewer(IMskelDil)

% ####### TO DO!
% save BW, IMskel, IMloopsclean, Gnew, Lep, Lbp, LBranchMidPoints,
% cyclesFilt
% save in loop table: extract # of branch points per loop, end points per
% loop? mean Geodesic distance to nearest end point? Local concentration of endpoints? 
% INS cells!_!_!_??
% - loop network? Usefull? (every loop a node and put in edges when loops
% share an edge?)... measure loop node centrality? Other node properties?
% (there is one movie where this could be done over time for a small loop
% network...) - are you less likely to close of you have a loop neiborg?
% Probably!?


%% WHOLE ORGAN CONVEX HULL
IMholeOrganCVH = logical(zeros(size(BW)));

tableHoleOrganCVH = regionprops3(uint8(BW).*2,"ConvexImage","BoundingBox");
BB = tableHoleOrganCVH(2,:).BoundingBox;
IMholeOrganCVH(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1)=tableHoleOrganCVH.ConvexImage{2}; 

DIMholeOrganCVH = bwdist(imcomplement(IMholeOrganCVH));

implaypair(DIMholeOrganCVH,BW)

%volumeViewer([IMholeOrganCVH,BW])



%%

TEMP = IMloopsclean;
TEMP(IMloopsCVHmask<6)=0;

implaypair(IMloops10biggest,IMloopsCVHmask>1)

%ismember(biggest,nonzeros(unique(TEMP)))

%%
implaypair(IMloops,IMloopsCVHmask>1)
%%
 
L = IMloopsCVH;
L(BW>0)=0;% remove labels inside lumen
%L(IMloopsCVHmask>1)=0;
L2 = imopen(L,strel('cube',2));
L2 = cleanUpLabelledIm(L2,0);


volumeViewer([IMskelDil,IMskelDil,BW])
%Ltot = [L,IMloopsCVHfilt];

Ltot = [IMloopsCVH,L2,L2];

%%
LL = int8(L2);
%%

volumeViewer(uint8(IMskelDil) + imdilate(3.*uint8(IMloops==index(end-20) | IMloops==index(end)),strel('sphere',5)))



%%


%%
r2 = regionprops3(L2,'Volume','PrincipalAxisLength');
r3 = regionprops3(IMloopsCVH,'Volume');

D = bwdist(imcomplement(BW));
r1 = regionprops3(IMloopsClean,D,'Volume','PrincipalAxisLength','VoxelValues');

LoopsTable = table(); 

LoopsTable.DeviationFromFlatCirclePlane = r2.Volume ./(pi.*(r1.Volume./(2.*pi)).^2); 
LoopsTable.Perim = r1.Volume; 
LoopsTable.VolConvexHull = r3.Volume; 
LoopsTable.VolDonutHole = r2.Volume; % if deviation from flat circle plane is large then VolDonutHole may be including parts of outside of loop... If so DO NOT use measurments extracted from dunut hole

LoopsTable.PrinAxis1SkelLoop = r1.PrincipalAxisLength(:,1);  
LoopsTable.PrinAxis2SkelLoop = r1.PrincipalAxisLength(:,2);  
LoopsTable.PrinAxis3SkelLoop = r1.PrincipalAxisLength(:,3);  

LoopsTable.PrinAxis1DonutHole = r2.PrincipalAxisLength(:,1);  
LoopsTable.PrinAxis2DonutHole = r2.PrincipalAxisLength(:,2);  
LoopsTable.PrinAxis3DonutHole = r2.PrincipalAxisLength(:,3);  

LoopsTable.MeanLumenRadius = cellfun(@mean,r1.VoxelValues);
LoopsTable.MedianLumenRadius = cellfun(@median,r1.VoxelValues);
LoopsTable.ModeLumenRadius = cellfun(@mode,r1.VoxelValues);
LoopsTable.STDLumenRadius = cellfun(@std,r1.VoxelValues);
LoopsTable.ExcessKurtosisLumenRadius = cellfun(@kurtosis,r1.VoxelValues)-3;

%LoopsTable.LumenRadiusDistSkewness = (LoopsTable.MeanLumenRadius - LoopsTable.ModeLumenRadius)./LoopsTable.STDLumenRadius;
% More robust formula for skewness for small data sets: (cite https://towardsdatascience.com/skewness-kurtosis-simplified-1338e094fc85)
LoopsTable.LumenRadiusDistSkewness = 3.*(LoopsTable.MeanLumenRadius - LoopsTable.MedianLumenRadius)./LoopsTable.STDLumenRadius;


LoopsTable.approxRho = LoopsTable.PrinAxis1DonutHole./(LoopsTable.MeanLumenRadius.*2);% ratio between diameter of donut hole and diameter of surrounding lumen (higher thin, smaller fat loop)
LoopsTable.approxRho_min = LoopsTable.PrinAxis2DonutHole./(LoopsTable.MeanLumenRadius.*2);% I assume here that third prin axis will be perpendicular to the plane of the loop 

LoopsTable

%%
implaypair(IMloopsCVHmask==2,IMloops)

%% close all
figure
nexttile
imshowpair(max(IMskel,[],3),max(IMloops>0,[],3))
nexttile
imshow(max(IMloops,[],3),[1 30]);colormap("jet")

%%
close all

rloop = regionprops3(IMloops,'PrincipalAxisLength','Volume');

PALmax = rloop.PrincipalAxisLength(:,1);

PALmax(PALmax==0)=[];

perimeter = rloop.Volume;% voxel count is a form of pseudo perimeter....
perimeter(perimeter==0)=0;


figure
histogram(PALmax,30)
figure
histogram(perimeter,30)


%%
Gnew




