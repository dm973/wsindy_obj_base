function [lib_true,true_S] = true_lib(nstates,true_nz_weights)
    true_S = cellfun(@(tnz)[(1:size(tnz,1))' tnz(:,end)],true_nz_weights,'Un',0);
    lib_true = repelem(library('nstates',nstates),length(true_nz_weights),1);
    for i=1:length(true_nz_weights)
        for j=1:size(true_nz_weights{i},1)
            lib_true(i).add_terms(term('ftag',true_nz_weights{i}(j,1:nstates),'linOp',true_nz_weights{i}(j,nstates+1:end-1)));
        end
    end
end