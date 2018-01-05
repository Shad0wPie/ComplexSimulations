%% Network1
A = load('Network1.txt');
A = sparse(A(:,1),A(:,2),1);

%% Network2
A = load('Network2.txt');
A = sparse(A(:,1),A(:,2),1);

%% Network3
A = load('Network3.txt');
A = sparse(A(:,1),A(:,2),1);

%% Calculate lengths
pathLengths = FindPathLengths(A);

%% Network 1 precalculated
A = load('Network1.txt');
A = sparse(A(:,1),A(:,2),1);
pathLengths = load('network1_path_lengths.mat');
pathLengths = pathLengths.pathLengths;

%%
graphSize = size(A,2);
averageLength = sum(sum(pathLengths))/(graphSize*(graphSize-1));
diameter = max(max(pathLengths));
fprintf("Average path length: %.6f\n", full(averageLength));
fprintf("Diameter: %d\n", full(diameter));

%%
clusteringCoefficient = CalculateClustering(A);
fprintf("Clustering coefficient: %.6f\n", full(clusteringCoefficient));

%%
edges = full(sum(A,2));
ks = min(edges):max(edges);
measured_fs = arrayfun(@(k) sum(edges>=k)/graphSize,ks);

loglog(ks,measured_fs,'b')
title(sprintf('Network 3. clustering=%.6f, avg-path=%.6f, diameter=%d',clusteringCoefficient, averageLength, diameter))
xlabel('k')
ylabel('F(k)')

%%
figure(2)
coords = [cos(2*pi*(1:graphSize)/graphSize); sin(2*pi*(1:graphSize)/graphSize)]';
gplot(A,coords)
title(sprintf('Network 1, nodes: %d, max degree: %d', graphSize,max(edges)))
axis equal
