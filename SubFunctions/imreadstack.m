function [ tiff_stack ] = imreadstack(filename)

tiff_info = imfinfo(filename); % return tiff structure, one element per image
temp_tiff = imread(filename, 1) ; % read in first image

tiff_stack = uint8(zeros(size(temp_tiff,1),size(temp_tiff,2),size(tiff_info, 1)));

%add tiff to tiff_stack
for ii = 1:size(tiff_info, 1)
    temp_tiff = imread(filename, ii);
    tiff_stack(:,:,ii) = temp_tiff;%cat(3 , tiff_stack, temp_tiff);
end

end

