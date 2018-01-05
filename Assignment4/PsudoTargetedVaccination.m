function largestClusterSize = PsudoTargetedVaccination(A, f)    
    graphSize = size(A,1);
    numNodesToRemove = fix(graphSize*f);
    
    nodesToRemove = zeros(1,numNodesToRemove);
    
    allNodes = 1:graphSize;
    
    i=1;
    while i <= numNodesToRemove
        notRemovedNodes = allNodes(~ismember(allNodes,nodesToRemove));
        index = fix(rand*length(notRemovedNodes)) + 1;
        neighbourIndices = find(A(index,:));
        if isempty(neighbourIndices)
            continue
        end
        notRemovedNeighbours = neighbourIndices(~ismember(neighbourIndices,nodesToRemove));
        if isempty(notRemovedNeighbours)
            continue
        end
        randIndex = fix(rand*length(notRemovedNeighbours)) + 1;
        randomNeighbour = notRemovedNeighbours(randIndex);
        nodesToRemove(i) = randomNeighbour;
        i=i+1;
    end
    
    A(nodesToRemove,:) = [];
    A(:,nodesToRemove) = [];
    
    largestClusterSize = FindLargestCluster(A);

end

