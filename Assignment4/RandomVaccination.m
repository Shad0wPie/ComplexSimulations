function largestClusterSize = RandomVaccination(A, f)    
    graphSize = size(A,1);
    
    numNodesToRemove = fix(graphSize*f);
    indices = randperm(graphSize);
    nodesToRemove = indices(1:numNodesToRemove);
    A(nodesToRemove,:) = [];
    A(:,nodesToRemove) = [];

    largestClusterSize = FindLargestCluster(A);
end

