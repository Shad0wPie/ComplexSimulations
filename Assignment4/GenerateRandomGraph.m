function A = GenerateRandomGraph( graphSize, density )

A = sprand(graphSize, graphSize, density);
A = triu(A);
A = A + A';
A(A~=0) = 1;
A(logical(eye(size(A)))) = 0;



end

