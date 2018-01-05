%% Parameters
n0 = 5;
d0 = 0.1;
m = 2;
timesteps = 50;
%% Generate random graph
A = sparse(ones(n0,n0));
A(logical(eye(size(A)))) = 0;

%% Preferential growth
graphSize = n0;

for t=1:timesteps
    edges = full(sum(A,2));
    connections = randsample(1:graphSize, m, true, edges);
    for connection=connections
        A(graphSize+1, connection) = 1;
        A(connection, graphSize+1) = 1;
    end
    graphSize = graphSize + 1;
end

%%
figure(2)
coords = [cos(2*pi*(1:graphSize)/graphSize); sin(2*pi*(1:graphSize)/graphSize)]';
gplot(A,coords)
title(sprintf('Connectivity graph (n=%d, m=%d, n0=%d)',graphSize,m,n0))
axis equal

%%
theoretical_func = @(k)(2*m^2*k.^(-2));

edges = full(sum(A,2));
kmin = min(edges);
kmax = max(edges);

ks = kmin:kmax;
ks_smooth = linspace(kmin,kmax);
measured_fs = arrayfun(@(k) sum(edges>=k)/graphSize,ks);

figure(1)
loglog(ks,theoretical_func(ks),'b')
hold on 
loglog(kmin:kmax,measured_fs, 'ro')
hold off
title(sprintf('Inverse cumulative degree distribution (n=%d, m=%f)',graphSize,m))
xlabel('k')
ylabel('F(k)')
% legend('Measured', 'Theoretical')
