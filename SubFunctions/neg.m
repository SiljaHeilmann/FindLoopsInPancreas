function [ IM_out ] = neg(IM_in )
% returns negative of image

IM_in = -1.*double(IM_in);

IM_out = imrescale(IM_in);

end

