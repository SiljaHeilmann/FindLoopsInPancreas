function [zyx] = flipxz(xyz)

zyx = permute(xyz,[3 1 2]);

end
