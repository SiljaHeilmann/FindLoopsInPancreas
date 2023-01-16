function [Gsimple] = skel2graphSimple(IMskel)

G = skel2graph(IMskel);

C = centrality(G,"degree");

G.Nodes.Degree = C;
G.Nodes.EndNode = G.Nodes.Degree==1;% will have a one if its an end node
G.Nodes.BranchNode = G.Nodes.Degree>2;% will have a one if its an branch nodes

subNodeIDs = G.Nodes.Index(G.Nodes.BranchNode==1 | G.Nodes.EndNode==1,:);% the nodes we want to use in simplified network bn's and en's

vec1 = single(G.Nodes.Degree(G.Edges.EndNodes(:,1),:)==2);% vec with 1 everywhere where first node has degree 2
vec2 = single(G.Nodes.Degree(G.Edges.EndNodes(:,2),:)==2);% vec with 1 everywhere where second node has degree 2

edgeWeights = single(vec1+vec2 ~= 2);% 1 if above 2 or below 2, 0 if 2... ( since we want edge weigths to be 0 between two nodes that both have degree 2)

G.Edges.Weight = edgeWeights;

d1 = distances(G,subNodeIDs,subNodeIDs,'Method','positive');% distances along paths of nodes with degree 2 will have zero weigth

% make ajacency matrix
A = d1<3;

Aupper = triu(A);
Aupper = Aupper - diag(ones(1,size(A,1)));

Edges = table(zeros(sum(Aupper(:)),2),'VariableNames',{'EndNodes'});

count = 1;
for cc = 1:size(Aupper,1)
    rr = find(Aupper(:,cc));
    if isempty(rr)==0
        Edges.EndNodes(count:count+size(rr,1)-1,:) = [ones(size(rr,1),1).*cc,rr];
        count = count+size(rr,1);
    end
end

Nodes = G.Nodes(subNodeIDs,:);

Gsimple = graph(Edges,Nodes,'omitselfloops');

% figure
% plot(Gsimple,'XData',Gsimple.Nodes.PixelCoor(:,1),'YData',Gsimple.Nodes.PixelCoor(:,2),'ZData',Gsimple.Nodes.PixelCoor(:,3)); hold on



end