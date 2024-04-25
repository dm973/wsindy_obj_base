function t3 = prodtags(t1,t2)
    if isempty(t2)
        if ~isempty(t1)
            if isequal(class(t1),'double')
                t1 = mat2cell(t1,ones(size(t1,1),1),size(t1,2));
            end
            t3 = t1;
        else
            t3 = [];
        end
    elseif isempty(t1)
        if isequal(class(t2),'double')
            t2 = mat2cell(t2,ones(size(t2,1),1),size(t2,2));
        end
        t3 = t2;
    else
        if isequal(class(t1),'double')
            t1 = mat2cell(t1,ones(size(t1,1),1),size(t1,2));
        end
        if isequal(class(t2),'double')
            t2 = mat2cell(t2,ones(size(t2,1),1),size(t2,2));
        end
        t3 = cell(length(t1),length(t2));
        for i=1:length(t1)
            for j=1:length(t2)
                t3{i,j} = [t1{i};t2{j}];
            end
        end
        t3=reshape(t3',[],1);
    end
end
