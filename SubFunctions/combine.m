function [col_im] = combine(im,bw)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
col_im = uint8(zeros(size(im,1),size(im,2),3));


neg = (bw-1).*1;
im_forground = im;
im_forground(bw==0)=0;

im_background = im;
im_background(bw==1)=0;

col_im(:,:,1)=im+uint8(uint8(bw).*max(im(:)./2).*0.5)-im_background.*1.3;
col_im(:,:,2)=im+uint8(uint8(bw).*max(im(:)./2).*0.5)-im_forground;
col_im(:,:,3)=im;

%col_im(:,:,1)=im+uint8(uint8(bw).*max(im(:)./2).*0.5)-im_forground; %RED
%col_im(:,:,2)=im+uint8(uint8(bw).*max(im(:)./2).*0.5)-im_background; %GREEN
%col_im(:,:,3)=im;

%col_im(:,:,1)=im_background.*50; %RED
%col_im(:,:,2)=im_forground.*20; %GREEN
%col_im(:,:,3)=zeros(size(im));

end

