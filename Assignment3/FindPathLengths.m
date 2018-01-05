function pathLengths = FindPathLengths(adjacencyMatrix)
    % Floyd–Warshall
    
    graphSize = size(adjacencyMatrix,1);
    pathLengths = adjacencyMatrix+Inf;
    pathLengths(adjacencyMatrix~=0) = 1;
    pathLengths(logical(eye(size(pathLengths)))) = 0;

    for k = 1:graphSize
        k
        newLengths = pathLengths(:,k) + pathLengths(k,:);
        pathLengths = min(newLengths, pathLengths);
    end
    
end