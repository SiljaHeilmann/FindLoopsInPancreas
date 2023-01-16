function D = tform2D(tform,moving)
%Takes a tform object and returns the corresponding displacement field


    % make a displacement field D that corresponds to tform
    [X,Y] = meshgrid(1:size(moving,2),1:size(moving,1)); % make grids with x and y coordinates of pixels
    [XT,YT] = transformPointsForward(tform,X,Y); % Transform pixel coordinates using tform
    D = zeros(size(X,1),size(X,2),2);
    D(:,:,1) = X-XT; % displacement field is the difference
    D(:,:,2) = Y-YT;


end

