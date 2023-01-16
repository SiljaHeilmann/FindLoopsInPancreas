function [Gnew] = removeSmallCycles(G,thresh)

[cycles, edgecycles ] = cyclebasis(G);

count=0;
indexToRemove =[];
for ii = 1:length(cycles)
    if length(cycles{ii})<=thresh
        count = count +1;
        indexToRemove(count)=ii;
    end
end

smallEdgeCycles = edgecycles(indexToRemove);

bigEdgeCycles = edgecycles;
bigEdgeCycles(indexToRemove)=[];

BEC = [bigEdgeCycles{:}];
SEC = [smallEdgeCycles{:}];

BECuni = unique(BEC);
SECuni = unique(SEC);

inboth = intersect(SECuni,BECuni);
edgesThatCanBeRemoved = setxor(SECuni,inboth);

G.Edges.Remove = ones(size(G.Edges.Dist)); % make new property and set all to 1
G.Edges.Remove(edgesThatCanBeRemoved)=0; % set edges in small cycles that do not overlap with big cycles to 0

EdgesNew = G.Edges(G.Edges.Remove==1,:); % make new edge table with out the edges from small loops

Gnew = graph(EdgesNew,G.Nodes,"omitselfloops"); % make new graph with new edges and nodes from old graph

end