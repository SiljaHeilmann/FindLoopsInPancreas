function [] = implaylabeloverlay(mov1,mov2,alpha)

if size(mov1,1)==size(mov2,1) && size(mov1,2)==size(mov2,2) && size(mov1,3)==size(mov2,3) && size(mov1,4)==size(mov2,4)

if size(mov1,3)==1 % if third dim is 1 time must be in the fourth dim
   mov1 = permute(mov1,[1 2 4 3]); % permute so that time becomes third dim
end

if size(mov2,3)==1
   mov2 = permute(mov2,[1 2 4 3]); 
end

% initialize empty color (rgb) stack (uint8)
col = uint8(zeros(size(mov1,1),size(mov1,2),size(mov1,3),3)); % time is in third dim

mov2(1,1,:) = 1; % trick so that labels have same colors over time...
mov2(1,2,:) = max(mov2(:)); % trick so that labels have same colors over time...

for tt=1:size(mov1,3)  
    col(:,:,tt,:) = labeloverlay(mov1(:,:,tt),mov2(:,:,tt),'Transparency', alpha);% label is the second input
end

col(1,1:2,:,:) = 0; % trick so that labels have same colors over time...

implayS(permute(col,[1 2 4 3]),1)% permute to put time back in 4th dim

else
    disp('Error - movies are not same size!')

end
end
