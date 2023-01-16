function [ IM_out ] = imlog(IM_in )
% RESCALES BY LOG TRANSFORM uint8, uin16 or double IMAGE SO THAT FULL RANGE IS USED 
% returns uint8 in case of double...

IM_in = log(double(IM_in));

IM_in(isfinite(IM_in)==0)=0; % remove -Inf, set to zero

IM_out = imrescale(IM_in);

end

