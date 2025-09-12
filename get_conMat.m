function C = get_conMat(lib,connected_terms)

%%% connected_terms is a cell array. Each cell is a list of terms whose
%%% weights W=[w1;...;wn] are connected by coefficients C=[c1;...;cn]
%%% according to W=Cv, where v is the scalar weight corresponding to the
%%% grouped terms in a given cell. Terms in each cell are listed in the form
%%% {eq, term_id, coef} to specify the equation the term appears in,
%%% the term itself as a tag or a string equal to output of term.get_str
%%% coef is the value of c.
%%% The output is a connectivity matrix C such that W=CV where V are the
%%% new (reduced weights).

%%% examples of connected terms for NLS:
% connected_terms ={ {{1,'D^[2,0]x2.^1',1}, {2,'D^[2,0]x1.^1',-1}}, ...
%                   {{1,'D^[0,0]x1.^2.*x2.^1',1}, {2,'D^[0,0]x1.^1.*x2.^2',-1}} };

% connected_terms ={ {{1,[0 1 2 0],1}, {2,[1 0 2 0],-1}}, ...
%                    {{1,[2 1 0 0],1}, {2,[1 2 0 0],-1}}, ...
%                    {{1,[0 3 0 0],1}, {2,[3 0 0 0],-1}}, ...
%                 };

    lib_sizes = arrayfun(@(L)length(L.terms),lib);
    nt = sum(lib_sizes);
    C = eye(nt);
    remove_cols = false(nt,1);
    strs=arrayfun(@(L)cellfun(@(t)t.get_str,L.terms,'un',0),lib,'un',0);
    tags=arrayfun(@(L)L.tags,lib,'un',0);
    for i=1:length(connected_terms)
        tt = connected_terms{i};
        foo = [];
        for j=1:length(tt)
            eq=tt{j}{1};
            term_id=tt{j}{2};
            coef=tt{j}{3};
            if ischar(term_id)
                ind=find(cellfun(@(s)isequal(term_id,s),strs{eq}));
            elseif isnumeric(term_id)
                ind=find(cellfun(@(t)isequal(term_id,t),tags{eq}));
            end
            if ~isempty(ind)
                foo = [foo,sum(lib_sizes(1:eq-1))+ind];
                C(foo(end),foo(1)) = coef;
                if length(foo)>1
                    remove_cols(foo(end)) = true;
                end
            end
        end
    end
    C = C(:,~remove_cols);
end