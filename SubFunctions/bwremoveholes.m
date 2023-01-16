function [bw_im_out] = bwremoveholes(bw_im_in,hole_size_thresh,conn)
% Takes BW image and fills holes that are smaller than a certain threshold

% bw_im_in - binary input image, 2D or 3D
% hole_size_thresh - threshold for removing holes
% conn - connectivity of hole region (optional - default is 8 for 2D og 26 for 3D)

dim = size(size(bw_im_in),2); % what is the dimension of bw_im_in?

if ~exist('conn','var') %if conn parameter does not exist, default it to something
    if dim == 2
        conn = 8;
    elseif dim == 3
        conn = 26;
    end
end


neg = bw_im_in.*-1 +1; % invert black and white

neg = bwareaopen(neg,hole_size_thresh,conn); % remove small white objects

bw_im_out = neg.*-1 +1; % invert black and white

end

