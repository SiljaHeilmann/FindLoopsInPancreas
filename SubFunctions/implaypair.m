function [] = implaypair(mov1,mov2)

if size(mov1,1)==size(mov2,1) && size(mov1,2)==size(mov2,2) && size(mov1,3)==size(mov2,3) && size(mov1,4)==size(mov2,4)

if size(mov1,3)==1 % if third dim is 1 time must be in the fourth dim
   mov1 = permute(mov1,[1 2 4 3]); % permute so that time becomes third dim
end

if size(mov2,3)==1
   mov2 = permute(mov2,[1 2 4 3]); 
end

% initialize empty color (rgb) stack (uint8)
col = uint8(zeros(size(mov1,1),size(mov1,2),size(mov1,3),3)); % time is in third dim

% normalize signal in red channel
col_red = double(mov1);
col_red = uint8(((col_red-min(col_red(:)))./(max(col_red(:))-min(col_red(:)))).*255);

% normalize signal in green channel 
col_green = double(mov2);
col_green = uint8(((col_green-min(col_green(:)))./(max(col_green(:))-min(col_green(:)))).*255);

col(:,:,:,1) = col_red;
col(:,:,:,2) = col_green;


implayS(permute(col,[1 2 4 3]),1)

else
    disp('Error - movies are not same size!')

end
end
