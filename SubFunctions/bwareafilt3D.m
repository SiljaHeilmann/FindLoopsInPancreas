function BWout = bwareafilt3D(BWin,range)
% remove connected regions where number of voxels is above or below the
% range

BWin = logical(BWin); % make sure BWin has type logical

smallest_remain =  logical(BWin - bwareaopen(BWin,range(1))); 

largest_remain = bwareaopen(BWin,range(2));

BWout = logical(BWin - smallest_remain - largest_remain);

end

