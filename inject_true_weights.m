function w_true = inject_true_weights(WS,true_nz_weights)
    w_true = WS.reshape_w('w',WS.weights*0);
    for i=1:WS.numeq
        tags = true_nz_weights{i};
        for j = 1:length(WS.lib(i).terms)
            tt = WS.lib(i).terms{j};
            try
                dt = tt.linOp.difftags;
            catch
                dt = zeros(1,WS.ndims);
            end
            ind = ismember(tags(:,1:end-1),[tt.ftag,dt],'rows');
            if any(ind)
                w_true{i}(j) = tags(ind,end); 
            end
        end
    end
    w_true = cell2mat(w_true);
end