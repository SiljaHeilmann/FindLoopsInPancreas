function voroimage = voronoizone2(x,y,img)

%[voroimage,sub_image]=voronoizone(x,y,img)


% Voronoi diagram based image zoning

% [voroimage,sub_image]=voronoizone(x,y,img)
% For a defined number of points on an image plane, voronoizone computes
% the voronoi diagram on the image space and divide the image into sub
% images according to the zones.

% Coded by Kalyan S Dash, IIT Bhubaneswar

szImg = size(img);
ONES = ones(szImg(1),szImg(2));
[COLS,ROWS] = meshgrid(1:szImg(2),1:szImg(1));

current_min  = ones(szImg(1),szImg(2)).*10000; % first layer is new current stack, second layer is previous stack
current_index  = zeros(szImg(1),szImg(2)); % first layer is new current stack, second layer is previous stack

for kk=1:length(x)

    DIST = (ROWS-ONES.*y(kk)).^2 +(COLS-ONES.*x(kk)).^2;

    current_min_bw = DIST<current_min; % put one everywhere where a dist is found to be smaller than previous values  
    
    current_index(current_min_bw)  = kk;
    
    current_min(current_min_bw) = DIST(current_min_bw); % put in dist where they where smaller
    
    
end

voroimage = current_index;

end