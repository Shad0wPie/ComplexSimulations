function [Q,clusters] = CalculateGraphSplit(A, trueK, trueM)
    k = full(sum(A,2));
    m = sum(k)/2;
    B = A - trueK*trueK'/(2*trueM);
    B = B - diag(sum(B,1));
%     B(logical(eye(size(B)))) = B(logical(eye(size(B)))) - (k - trueK*(2*m/(2*trueM)));
    [V,D] = eig(B);
    eigenValues = diag(D);
    [sortedValues,sortIndex] = sort(eigenValues(:),'descend');
    maxIndex = sortIndex(1);
    largestEigenVector = V(:,maxIndex);
    s = ones(size(largestEigenVector));
    s(largestEigenVector<=0) = -1;
    Q = 1/(4*trueM)*s'*B*s;
    if Q > 0 && length(s) > abs(sum(s))
        % perform the split and also split the subgraphs
        A1 = A(s>0,s>0);
        trueK1 = trueK(s>0);
        [Q1, cluster1] = CalculateGraphSplit(A1,trueK1, trueM);
        A2 = A(s<0,s<0);
        trueK2 = trueK(s<0);
        [Q2, cluster2] = CalculateGraphSplit(A2,trueK2, trueM);
        
        clusters = s(:);
        if Q1 > 0
            clusters(s>0) = clusters(s>0) + cluster1 + abs(min(cluster1));
        end
        if Q2 > 0
            clusters(s<0) = clusters(s<0) + cluster2 - max(cluster2);
        end
%         fprintf('%.6f\n',Q);
        Q = Q + Q1 + Q2;
%         fprintf('%.6f\n',Q);
    else
        % don't do the split
        clusters = ones(size(largestEigenVector));
        Q=0;
    end

end

