function [zxy] = flipyz(xyz)

zxy = permute(xyz,[3 2 1]);

end

