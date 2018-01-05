function A = GeneratePreferential(graphSize, m)

    A = sparse(ones(m,m));
    A(logical(eye(size(A)))) = 0;

    timesteps = graphSize-m;
    % Preferential growth
    graphSize = m;

    for t=1:timesteps
        edges = full(sum(A,2));
        connections = randsample(1:graphSize, m, true, edges);
        for connection=connections
            A(graphSize+1, connection) = 1;
            A(connection, graphSize+1) = 1;
        end
        graphSize = graphSize + 1;
    end

end

