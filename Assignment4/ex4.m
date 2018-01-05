clear all
graphSize = 400;
randomDensity = 0.01;
preferentialM = 3;

iterations = 10;

% A = GeneratePreferential(graphSize, preferentialM);
A = GenerateRandomGraph(graphSize, randomDensity);

k = full(sum(A,2));
averageK = mean(k);

degDist = histc(k, 1:max(k))/graphSize;
excessDegDist = zeros(size(degDist));
for i=1:(length(degDist)-1)
    excessDegDist(i) = (i + 1).*degDist(i+1)/averageK;    
end

h=bar([degDist,excessDegDist]);

l = cell(1,2);
l{1}='degree distribution'; l{2}='excess degree distribution';
legend(h,l);
title(sprintf('Degree distributions for random (nodes=%d, density=%.2f)',graphSize, randomDensity))
% title(sprintf('Degree distributions for preferential (nodes=%d, m=%d)',graphSize, preferentialM))
