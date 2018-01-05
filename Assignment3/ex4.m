A = load('smallWorldExample.txt');

clusteringCoefficient = CalculateClustering(A);
sprintf('Clustering coefficient is %.6f', clusteringCoefficient)

%%
graphSize = size(A,1);

figure(2)
coords = [cos(2*pi*(1:graphSize)/graphSize); sin(2*pi*(1:graphSize)/graphSize)]';
gplot(A,coords)
title(sprintf('Example graph (coefficient=%.6f)',clusteringCoefficient))
axis equal