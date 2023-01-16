function [BW] = holefilt3D(BW,min_hole_2Darea, max_hole_2Darea)
% Only keep 2D black holes in a certain size range in a 3D BW image


% first use normal 3D hole closing (no matter the size)
BW= imfill(BW,'holes');

for zz = 1:size(BW,3)
    temp = BW(:,:,zz);
    temp_neg = imcomplement(temp);
    temp_neg = bwareafilt(temp_neg,[min_hole_2Darea*2 max_hole_2Darea]);
    BW(:,:,zz) = imcomplement(temp_neg);    
end


end

%BW= imfill(BW,'holes');

% for xx = 1:size(BW,2)
%     temp = BW(:,xx,:);
%     temp_neg = imcomplement(temp);
%     temp_neg = bwareafilt(permute(temp_neg,[1 3 2]),[min_hole_2Darea max_hole_2Darea]);
%     BW(:,xx,:) = imcomplement(permute(temp_neg,[1 3 2]));
% end
% 
% BW = imfill(BW,'holes');
% 
% for yy = 1:size(BW,1)
%     temp = BW(yy,:,:);
%     temp_neg = imcomplement(temp);
%     temp_neg = bwareafilt(permute(temp_neg,[2 3 1]),[min_hole_2Darea max_hole_2Darea]);
%     BW(yy,:,:) = imcomplement(permute(temp_neg,[3 1 2]));
% end
% 
% BW = imfill(BW,'holes');
