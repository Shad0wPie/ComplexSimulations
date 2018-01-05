clear all
graphSize = 400;
c = 4;
p = 0.01;
R = 0.1;
randomDensity = 0.01;
preferentialM = 3;
iterations = 10;

fs = 0:0.01:0.6;
clustersizesRandom = zeros(iterations,length(fs));
clustersizesTargeted = zeros(iterations,length(fs));
clustersizesRandomTargeted = zeros(iterations,length(fs));

for iteration=1:iterations
%     A = GenerateSmallWorld(graphSize, c, p);
    A = GenerateRandomGraph(graphSize, randomDensity);
%     A = GeneratePreferential(graphSize, preferentialM);
    for i=1:length(fs)
        clustersizesRandom(iteration,i) = RandomVaccination(A,fs(i));
        clustersizesTargeted(iteration,i) = TargetedVaccination(A,fs(i));
        clustersizesRandomTargeted(iteration,i) = PsudoTargetedVaccination(A,fs(i));
    end
end
clustersizesRandom = sum(clustersizesRandom,1)/iterations;
clustersizesTargeted = sum(clustersizesTargeted,1)/iterations;
clustersizesRandomTargeted = sum(clustersizesRandomTargeted,1)/iterations;

plot(fs,clustersizesRandom/graphSize, 'r')
hold on
plot(fs,clustersizesTargeted/graphSize, 'b')
title(sprintf('Random graph, nodes=%d, density=%.2f', graphSize, randomDensity))
plot(fs,clustersizesRandomTargeted/graphSize, 'g')
hold off
legend(["Random","Targeted","Randomly Targeted"])

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
