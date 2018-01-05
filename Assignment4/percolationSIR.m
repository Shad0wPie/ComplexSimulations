function largestClusterSize = percolationSIR(A, R)
    graphSize = size(A,1);
    for i=1:graphSize
        for j=i:graphSize
            psi = 1/(1+1/R);
            if rand < 1- psi
                A(i,j)=0;
                A(j,i)=0;
            end
        end
    end
    largestClusterSize = FindLargestCluster(A);
end