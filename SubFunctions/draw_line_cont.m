function M = draw_line_cont(a,BW)

B = bwboundaries(BW);

imshow(a); title('Draw outline of epithelial tissue - click mouse for do over, klick key when happy')

hold on

if isempty(B)==0
    for kk=1:length(B)
        b=B{kk};
        plot(b(:,2),b(:,1),'-r')
    end
end

set(gca,'nextplot','replacechildren')

h = imfreehand(gca);

k = waitforbuttonpress; % k=1 key press, k=0 mouse press

while k==0 % mouse press when unhappy with region picked keeps you in loop
    close all
%    imagesc(a);
    imshow(a); title('Draw outline of epithelial tissue - click mouse for do over, klick key when happy')
    hold on;
    if isempty(B)==0
        for kk=1:length(B)
            b=B{kk};
            plot(b(:,2),b(:,1),'-r')
        end
    end

    set(gca,'nextplot','replacechildren')
    h = imfreehand(gca);
    k = waitforbuttonpress; % k=1 key press, k=0 mouse press
end

p = round(getPosition(h));
close all

M = poly2mask(p(:,1),p(:,2),size(a,1),size(a,2));
M(BW==1) = 1;
