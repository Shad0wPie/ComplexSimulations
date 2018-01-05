function A = GenerateSmallWorld(graphSize, c, p)

    A = sparse(graphSize, graphSize);
    for i=1:graphSize
        for neighbour=1:c/2
            upperIndex = i+neighbour;
            if upperIndex > graphSize
                upperIndex = upperIndex - graphSize;
            end
            lowerIndex = i-neighbour;
            if lowerIndex < 1
                lowerIndex = lowerIndex + graphSize;
            end
            A(i,upperIndex) = 1;
            A(i,lowerIndex) = 1;
        end
    end

    for i=1:graphSize*c
        done = false;
        while ~done
            if rand < p
                ind1 = fix(rand*graphSize)+1;
                ind2 = fix(rand*graphSize)+1;
                if ind1 ~= ind2
                    A(ind1,ind2) = 1;
                    A(ind2,ind1) = 1;
                    done = true;
                end
            else
                done = true;
            end
        end 
    end

    A(logical(eye(size(A)))) = 0;

end
