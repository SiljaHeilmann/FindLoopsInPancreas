function [Lclean] = cleanUpLabelledIm(L,imcloseThres)


% clean up L - make sure only to include biggest connected component for each index. Maybe not nesseary - but in case rearrange for loops overlap left a few pixels far from rest of loop!
L = uint16(L);
Lclean = uint16(L);
Lclean(Lclean>0)=0;% make empty labelled image ready for output
TEMPbig = Lclean;% Initialize empty image used in loop

s = regionprops3(L,'Image','BoundingBox');


for ii = 1:uint16(max(L(:)))
    TEMPbig(TEMPbig>0)=0;
    
    if iscell(s(ii,:).Image)
        TEMPsmall = uint16(cell2mat(s(ii,:).Image));
    else
        TEMPsmall = uint16(s(ii,:).Image);
    end

    BB = s(ii,:).BoundingBox;

    if imcloseThres>0
        TEMPsmallDil = bwlabeln(imdilate(TEMPsmall,strel('cube',imcloseThres)),26);% Dilated and labelled version of image
    elseif imcloseThres==0
        TEMPsmallDil = bwlabeln(TEMPsmall,26);% only labelled (for when you want to clean up with no dilation)
    end
    r = regionprops3(TEMPsmallDil,'Volume');
    [~,index] = sort(r.Volume(:));

    if isempty(index)==0
        indexkeep = index(end);% find the largest object and keep only that
        MASK=TEMPsmallDil==indexkeep;% make mask with largest region only
        TEMPsmall(MASK==0)=0;% remove everything but the largest object

        TEMPbig(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1) = TEMPsmall.*ii;% add it in to TEMPbig

        Lclean = Lclean + TEMPbig;% add TEMPbig on to Lclean
    end
end

end