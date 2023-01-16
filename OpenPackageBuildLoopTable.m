
%% IMPORT iso images and original prediction file (H5 file)

ff = 4;
path = '/Users/xqk313/Documents/Piljabi-LOCAL/H5files-org-prob-seg';

filename = ['ff=' num2str(ff) '_DownScaled_Simple Segmentation.h5'];

INFO = h5info([path filesep filename])

filename2 = ['ff=' num2str(ff) '_DownScaled.h5'];
INFO2 = h5info([path filesep filename2])

imsizeP = INFO.Datasets.Dataspace.Size% x y z c
imsizeIM = INFO2.Datasets.Dataspace.Size% x y z c

start = [1 1 1 1];

countP = imsizeP;%[500 500 100 3];
countIM = imsizeIM;%[500 500 100 3];

P = h5read([path filesep filename],['/' INFO.Datasets.Name],start,countP);
IM = h5read([path filesep filename2],['/' INFO2.Datasets.Name],start,countIM);


%% LOAD PACKAGE FOR MAKING LOOP TABLE


BW = logical(BW);
% make BW logical!!!


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
    loopTable.BranchPointsIndex(ii)  = {nonzeros(unique(Lbp(Gnew.Nodes.LiniarIndexInImage(cyclesFilt{ii}))))};% Index of branch points in loop
    loopTable.NumBranchPoints(ii)    = size(cell2mat(loopTable.BranchPointsIndex(ii)),1);% number of branch points in loop

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

ff=7

save(['/Users/xqk313/Library/CloudStorage/GoogleDrive-heilmann.silja@gmail.com/My Drive/Piljabi/QuantifyLoopInWM/loopTables/loopTable_ff' num2str(ff) '.mat'],'loopTable');

%save(['/
% 
% Volumes/GoogleDrive/My Drive/Piljabi/QuantifyLoopInWM/loopTable_ff' num2str(ff) '.mat'],'loopTable');

loopTable

%%
size(nonzeros(loopTable.NumMesenGoThru(~isnan(loopTable.NumMesenGoThru))))

%%

