clear all
graphSize = 3000;
p = 0.1;
c = 4;

A = sparse(graphSize, graphSize);
for i=1:graphSize
    for neighbour=1:c/2
        upperIndex = i+neighbour;
        if upperIndex > graphSize
            upperIndex = upperIndex - graphSize;
        end
        lowerIndex = i-neighbour;
        if lowerIndex < 1
            lowerIndex = lowerIndex + graphSize;
        end
        A(i,upperIndex) = 1;
        A(i,lowerIndex) = 1;
    end
end

for i=1:graphSize*c
    if rand < p
        ind1 = fix(rand*graphSize)+1;
        ind2 = fix(rand*graphSize)+1;
        A(ind1,ind2) = 1;
    end
end

A(logical(eye(size(A)))) = 0;

%%
figure(2)
coords = [cos(2*pi*(1:graphSize)/graphSize); sin(2*pi*(1:graphSize)/graphSize)]';
gplot(A,coords)
title(sprintf('Connectivity graph (n=%d, c=%d, p=%f)',graphSize,c,p))
axis equal