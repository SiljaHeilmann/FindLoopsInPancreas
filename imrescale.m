function [ IM_out ] = imrescale(IM_in )
% RESCALES INTENSITY OF uint8, uin16 or double IMAGE SO THAT FULL RANGE IS USED 
% returns uint8 in case of double...

type = class(IM_in)

IM_in = double(IM_in);


switch type
    case 'uint8'
IM_out = uint8(((IM_in-min(IM_in(:)))/(max(IM_in(:))-min(IM_in(:))))*(2^8-1));% cast as uint8
    case 'uint16'
IM_out = uint16(((IM_in-min(IM_in(:)))/(max(IM_in(:))-min(IM_in(:))))*(2^16-1));% cast as uint8
    case 'double'
IM_out = uint8(((IM_in-min(IM_in(:)))/(max(IM_in(:))-min(IM_in(:))))*(2^8-1));% cast as uint8



% switch type
%     case 'uint8'
% IM_out = uint8(((IM_in-min(IM_in(:)))/(prctile(IM_in(:),99)-min(IM_in(:))))*(2^8-1));% cast as uint8
%     case 'uint16'
% IM_out = uint16(((IM_in-min(IM_in(:)))/(prctile(IM_in(:),99)-min(IM_in(:))))*(2^16-1));% cast as uint8
%     case 'double'
% IM_out = uint8(((IM_in-min(IM_in(:)))/(prctile(IM_in(:),99)-min(IM_in(:))))*(2^8-1));% cast as uint8



end

end

