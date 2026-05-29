function total_parts = partitionNk(N,k) 
    if N>1
       total_parts = partitions(N,ones(1,k));
    elseif N==1
        total_parts = eye(k);
    elseif N==0
        total_parts = zeros(1,k);
    end
end