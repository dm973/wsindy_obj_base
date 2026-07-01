%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this function creates a library suitable for PDE discovery with all 
%%%% combinations of poly and trig terms from the lists 'polys' and 'trigs',
%%%% tensor producted with the all spatial differential operators contained 
%%%% in the list array 'x_diffs'.
%%%% 'custom_add' is a cell array of custom terms to be added that fall 
%%%% outside this poly-trig library. custom_remove_f is a cell array of
%%%% functions of the form @(tag) bool(tag) such that bool maps tag to a
%%%% boolean, and if bool returns true, the term is removed. This is
%%%% strictly to remove terms from the poly-trig portion of the library.
%%%% custom_remove_t is a matrix of poly-trig tags to be explicitly removed

function lib = get_lib_pde(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t)    
    nstates = Uobj.nstates;
    ndims = Uobj.ndims;
    
    tags = get_tags(polys,trigs,nstates);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    diff_tags = diff_tags(diff_tags(:,end)==0,:);
    for j=1:size(tags,1)
        for i=1:size(diff_tags,1)
            if all([~and(sum(diff_tags(i,:))>0,...
                    isequal(tags(j,:),zeros(1,nstates))),...
                    ~cellfun(@(b)b([tags(j,:) diff_tags(i,:)]),custom_remove_f),...
                    ~ismember_rows([tags(j,:) diff_tags(i,:)],custom_remove_t)])
                lib.add_terms(term('ftag',tags(j,:),'linOp',diff_tags(i,:)));
            end
        end
    end
    lib.add_terms(custom_add);
    
end