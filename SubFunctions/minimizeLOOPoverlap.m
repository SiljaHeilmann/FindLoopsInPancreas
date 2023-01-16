function [cyclesNEW,numChanges] = minimizeLOOPoverlap(cycles)
% MATLAB cycle basis function tends to return cycles with overlap - this
% function rearranges cycles pairwise to minimize overlap
% cycles - cell array/list of cycles (list of node index)

numChanges = 0;
ii=0;
while ii<length(cycles)-1
    ii = ii+1; % we iterate through the list of cycles
    jj = ii;   % in order not to go through pairs twice
    while jj <length(cycles) 
        jj = jj+1; % we iterate through the list of cycles

        IoU  =  length(intersect(cycles{ii},cycles{jj}))./length(union(cycles{ii},cycles{jj}));% IoU of current cycle pair

        % if IoU>0 (there is overlap) and intersect is greater than either
        % of the setdiff s there is away to rearrange into two cycles with
        % so go and do that
        if  ii~=jj && IoU>0 && (length(intersect(cycles{ii},cycles{jj}))>length(setdiff(cycles{ii},cycles{jj})) || length(intersect(cycles{ii},cycles{jj})) > length(setdiff(cycles{jj},cycles{ii})))

            if length(setdiff(cycles{ii},cycles{jj}))<length(setdiff(cycles{jj},cycles{ii})) % setdiff(cycles{ii},cycles{jj}) will be the new shared part
                cyclesNew1 = union(setdiff(cycles{ii},cycles{jj}), intersect(cycles{ii},cycles{jj}));
                cyclesNew2 = union(setdiff(cycles{ii},cycles{jj}), setdiff(cycles{jj},cycles{ii}));
            else % setdiff(cycles{jj},cycles{ii}) will be the new shared part (shortest)
                cyclesNew1 = union(setdiff(cycles{jj},cycles{ii}), intersect(cycles{ii},cycles{jj}));
                cyclesNew2 = union(setdiff(cycles{ii},cycles{jj}), setdiff(cycles{jj},cycles{ii}));
            end

            cycles{ii} = cyclesNew1; % update cycles so it contains the new rearranges cycles
            cycles{jj} = cyclesNew2;
            ii = 0; % something was rearranged - start over! 
            numChanges = numChanges + 1; % update counter of number of changes 
            break % break out of jj while loop
        end
    end
end

cyclesNEW = cycles;% set output of function

end