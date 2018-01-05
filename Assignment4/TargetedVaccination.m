function largestClusterSize = TargetedVaccination(A, f)    
    graphSize = size(A,1);
    
    numNodesToRemove = fix(graphSize*f);
    ranks = full(sum(A,2));
    [sortedValues,sortIndex] = sort(ranks,'descend');
    nodesToRemove = sortIndex(1:numNodesToRemove);
    A(nodesToRemove,:) = [];
    A(:,nodesToRemove) = [];
    
    largestClusterSize = FindLargestCluster(A);

end

