graphSize = 400;
c = 4;
p = 0.01;
iterations = 20;

rs = 0:0.1:5;
clustersizes = zeros(iterations,length(rs));

for iteration=1:iterations
    A = GenerateSmallWorld(graphSize, c, p);
    for i=1:length(rs)
        clustersizes(iteration,i) = percolationSIR(A,rs(i));
    end
end
clustersizes = sum(clustersizes,1)/iterations;
plot(rs,clustersizes/graphSize)
title(sprintf("S-R graph (small world with %d nodes and c=%d)", graphSize, c))
xlabel("R")
ylabel("S")

%%

k = full(sum(A,2));
m = sum(k)/2;
[Q, clusters] = CalculateGraphSplit(A,k,m);
clusters = (clusters + abs(min(clusters)))/2 +1;

H = plot(graph(A), 'MarkerSize',10, 'NodeLabel', clusters);
colortable = jet(max(clusters));
colors = zeros(length(clusters),3);
for i=1:length(clusters)
    colors(i,:) = colortable(clusters(i),:);
end
H.NodeColor = colors;
