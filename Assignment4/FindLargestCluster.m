function largestClusterSize = FindLargestCluster(A)

    graphSize = size(A,1);
    
    clusterMap = zeros(graphSize,1);
    
    currentCluster = 1;
    while any(clusterMap==0)
        startNode = find(clusterMap==0,1);
        exploreNodeCluster(startNode)
        currentCluster = currentCluster + 1;
    end
    
    function exploreNodeCluster(currentNode)
        neighbours = find(A(currentNode,:));
        clusterMap(currentNode) = currentCluster;
        for neighbour=neighbours
            if clusterMap(neighbour)
                continue
            end
            exploreNodeCluster(neighbour)
        end
    end

    largestClusterSize = max(histc(clusterMap, unique(clusterMap)));
end

