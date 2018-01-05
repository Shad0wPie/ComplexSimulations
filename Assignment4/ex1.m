A = load('communityExample.txt');
A = sparse(A);

k = full(sum(A,2));
m = sum(k)/2;
[Q, clusters] = CalculateGraphSplit(A,k,m);
clusters = (clusters + abs(min(clusters)))/2 +1;

fprintf('%.6f\n',Q)

H = plot(graph(A), 'MarkerSize',10, 'NodeLabel', clusters);
title(sprintf('Grapph clustering, Q=%.6f',Q))

colortable = jet(max(clusters));
colors = zeros(length(clusters),3);
for i=1:length(clusters)
    colors(i,:) = colortable(clusters(i),:);
end
H.NodeColor = colors;

%%
FindLargestCluster(A)