function [ neg ] = invertbw(im)
% Iverts logical BW image

im = logical(im);

neg = logical((im-1).*-1);

end

