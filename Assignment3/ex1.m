graphSize = 10;
density = 0.01;

A = sprand(graphSize, graphSize, density);
A = triu(A);
A = A + A';
A(A~=0) = 1;
A(logical(eye(size(A)))) = 0;

edges = full(sum(A,2));
ks = min(edges):max(edges);
measured_ps = arrayfun(@(k) sum(edges==k)/graphSize,ks);
theoretical_ps = arrayfun(@(k) dist(k,graphSize,density),ks);

figure(1)
hold on
plot(ks,measured_ps,'b')
plot(ks,theoretical_ps, 'r')
title(sprintf('Degree distribution plot (n=%d, p=%f)',graphSize,density))
xlabel('k')
ylabel('P(k)')
legend('Measured', 'Theoretical')
hold off
%%
figure(2)
coords = [cos(2*pi*(1:graphSize)/graphSize); sin(2*pi*(1:graphSize)/graphSize)]';
gplot(A,coords)
title(sprintf('Connectivity graph (n=%d, p=%f)',graphSize,density))
axis equal