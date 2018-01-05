function clustering = CalculateClustering( A )
%CALCULATECLUSTERING Summary of this function goes here
%   Detailed explanation goes here

    distinctTriangles = trace(A^3)/6;
    
    ranks = sum(full(A),2);
    connectedTriplets = sum(ranks.*(ranks-1)/2);
    
    clustering = 3 * distinctTriangles/connectedTriplets;
end

