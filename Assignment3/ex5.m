%% small world example
A = load('smallWorldExample.txt');
A = sparse(A);

%%
graphSize = size(A,2);

pathLengths = full(FindPathLengths(A));

averageLength = sum(sum(pathLengths))/(graphSize*(graphSize-1));
diameter = max(max(pathLengths));

sprintf("Average path length: %.6f", averageLength)
sprintf("Diameter: %d", diameter)