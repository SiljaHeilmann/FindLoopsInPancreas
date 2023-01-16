function M = draw_line(a)

imshow(a); title('draw outline of epithelial tissue')
%imagesc(a);

set(gca,'nextplot','replacechildren')

h = imfreehand(gca);
%accepted_pos = wait(h);

k = waitforbuttonpress; % k=1 key press, k=0 mouse press

while k==0 % mouse press when unhappy with region picked keeps you in loop
    close all
%    imagesc(a);
    imshow(a); title('draw outline of epithelial tissue')

    hold on;
    set(gca,'nextplot','replacechildren')
    h=imfreehand(gca);
    k = waitforbuttonpress; % k=1 key press, k=0 mouse press
end

p=round(getPosition(h));
close all

M=poly2mask(p(:,1),p(:,2),size(a,1),size(a,2));
