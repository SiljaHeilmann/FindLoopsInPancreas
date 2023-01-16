function [BWout] = returnLargestConnComp(BW,conn)
% returns an image with only the largest connected component in it

LABELS = bwlabeln(BW,conn);

r = regionprops3(LABELS,'Volume');
[sorted,index] = sort(r.Volume);
biggestIndex = index(end);


BWout = LABELS==biggestIndex;


end