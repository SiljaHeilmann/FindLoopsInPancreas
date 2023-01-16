function [ ] = implayS(movie,side)

switch side
    case 1
        
        handle = implay(movie);
        handle.Visual.ColorMap.UserRangeMin = min(movie(:));
        handle.Visual.ColorMap.UserRangeMax = max(movie(:));
        handle.Visual.ColorMap.UserRange = 1;
        handle.Parent.Position = [0 0 800 1500];
        handle.Visual.Axes.Position = [0 0 500 500];

        
    case 2
        
        
        handle = implay(movie);
        handle.Visual.ColorMap.UserRangeMin = min(movie(:));
        handle.Visual.ColorMap.UserRangeMax = max(movie(:));
        handle.Visual.ColorMap.UserRange = 1;
        handle.Parent.Position = [1500 0 800 1500];
        handle.Visual.Axes.Position = [0 0 500 500];
        
        
end

end

