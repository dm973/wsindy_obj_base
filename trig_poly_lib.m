function [lib,true_S] = trig_poly_lib(polys,trigs,x_diffs,nstates,ndims,numeq,true_nz_weights)
    tags_poly = get_tags(polys,[],nstates);
    % tags_1 = tags_1(tags_1(:,3)>0,:);
    tags_trig = get_tags([],trigs,nstates);
    tags = prodtags(tags_poly,tags_trig);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    diff_tags = diff_tags(diff_tags(:,end)==0,:);
    % diff_tags = diff_tags(sum(diff_tags>0,2)<=1,:);
    true_S = cell(numeq,1);
    for i=1:size(diff_tags,1)
        for j=1:length(tags)
            if ~and(sum(diff_tags(i,:))>0,isequal(tags{j},zeros(1,nstates)))
                lib.add_terms(term('ftag',tags{j},'linOp',diff_tags(i,:)));
                if exist('true_nz_weights','var')
                    if ~isempty(true_nz_weights)
                        S = cell2mat(cellfun(@(tnz) ismember(tnz(:,1:end-1),[tags{j} diff_tags(i,:)],'rows'), true_nz_weights, 'Un',0));
                        for jj=1:size(S,2)
                            if any(S(:,jj))
                                true_S{jj} = [true_S{jj};[length(lib.terms) true_nz_weights{jj}(S(:,jj),end)]];
                            end
                        end
                    end
                end
            end
        end
    end
end