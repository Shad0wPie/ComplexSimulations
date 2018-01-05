function y = dist(k,n,p)
    y = nchoosek(n-1,k)*p^k*(1-p)^(n-1-k);
end

