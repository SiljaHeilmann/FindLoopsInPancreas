function [sub_im] = imcrop3D(im,boundingBox)
%imcrop for 3d stack
%boundingBox is upper left corner coor x,y,z [x y z x_width y_width z_width]
x = floor(boundingBox(1));
y = floor(boundingBox(2));
z = floor(boundingBox(3));
x_width = ceil(boundingBox(4))+1;
y_width = ceil(boundingBox(5))+1;
z_width= ceil(boundingBox(6))+1;

temp = im;

temp(:,:,z+z_width:end)=[]; % remove bottom
temp(:,:,1:z-1)=[]; % remove top

temp(:,x+x_width:end,:)=[]; % remove x side
temp(y+y_width:end,:,:)=[]; % remove y side

temp(:,1:x-1,:)=[]; 
temp(1:y-1,:,:)=[]; 

sub_im = temp;

end

