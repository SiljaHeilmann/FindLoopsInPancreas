function [rgb_stack_out] = labeloverlay3(Im3D,label3D,alpha)
% labels 3D image and returns rgb stack (n x m x 3 x k)

if isempty(alpha)
    alpha = 0.5;
end

label3D(1,1,:) = max(label3D(:)).*ones(1,size(label3D,3));% make sure a patch with the higesth label number is always in current slice (for color choice)

rgb_stack_out = uint8(zeros(size(label3D,1),size(label3D,2),3,size(label3D,3))); % initialize RGB stack

for zz=1:size(rgb_stack_out,4)
    rgb_stack_out(:,:,:,zz) = labeloverlay(Im3D(:,:,zz),label3D(:,:,zz),'Transparency',alpha);% use 2D labeloverlay
end

end

