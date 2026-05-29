% addpath(genpath('C:/Users/385570/Desktop/repos/wsindy_obj_base/'))
% 
% n = 2; 
% d = 2;
% tags = [1 2 1 0; 2 0 0 0];
% tags = get_tags(0:4,[],2)*[eye(2) zeros(2)];
% tags = tags(tags(:,2)<3,:);
% 
% [lib,true_S] = trig_poly_lib(0:3,[],0:2,n,d,1,{@(t)t(2)>3,@(t)t(1)<2},[],[]);
% tags = cell2mat(lib.tags(:));
% 
% s = [0.1,0.2];
% [Ainv,Alarge,supp,lib] = getAinv(tags,s);

% [tags_cell,sub_inds,original_tags] = get_sub_tags(true_nz_weights{eqs(eq)}(:,1:end-1),3);
% supp=get_original_support(tags_cell,sub_inds,original_tags);

function [A,Ainv,Alarge,supp,lib,tags] = getAinv(tags,s)
    n = length(s);
    [tags_cell,sub_inds,original_tags] = get_sub_tags(tags,n);
    Alarge = 1;
    for d=1:length(tags_cell)
        [Afun,Alarge_d,~] = get_momentA(tags_cell{d},n);
        Afun_mat = Afun(s(:));
        if d==1
            A = Afun_mat;
            Ainv = inv(Afun_mat);
        else
            A = blkdiag(A,Afun_mat);
            Ainv = blkdiag(Ainv,inv(Afun_mat));
        end
        Alarge = blkdiag(Alarge,Alarge_d);
    end
    supp=get_original_support(tags_cell,sub_inds,original_tags);
    [lib,tags] = getA_lib(tags_cell,original_tags);
end

function [lib,tags] = getA_lib(tags_cell,original_tags)
    n = size(tags_cell{1},2);
    tags = [];
    for d=1:length(tags_cell)
        tags = [tags;[tags_cell{d},repmat(original_tags{d}{1},size(tags_cell{d},1),1)]];
    end
    lib = library('nstates',n);
    lib.add_terms(tags);
end

function supp=get_original_support(tags_cell,sub_inds,original_tags)
    supp = []; L=0;
    for i=1:length(tags_cell)
        for j=1:size(original_tags{i}{2},1)
            supp = [supp;L+find(sub_inds{i}==j)];
        end
        L = L+ length(sub_inds{i});
    end
end

function [tags_cell,sub_cell,original_tags] = get_sub_tags(tags,n)
     diff_tags = tags(:,n+1:end);
     [~,b,c] = unique(diff_tags,'rows');
     num_ops = max(c);
     tags_cell = cell(num_ops,1);
     sub_cell = cell(num_ops,1);
     original_tags = cell(num_ops,1);
     for i=1:max(c)
         tags_temp = tags(c==i,1:n);
         tags_all = [];
         for j=1:size(tags_temp,1)
             tags_all = [tags_all;get_sub_tag(tags_temp(j,:))];
             tags_all = unique(tags_all,'rows');
         end         
         % get_tags(0:max(sum(tags_temp,2)),[],n);
         % for nn=1:n
         %     tags_all = tags_all(tags_all(:,nn)<=max(tags_temp(:,nn)),:);
         % end
         [~,sub_cell{i}] = ismember(tags_all,tags_temp,'rows');
         tags_cell{i} = tags_all;
         original_tags{i} = {diff_tags(b(i),:),tags_temp};
     end
end

function tags = get_sub_tag(tag)
    tags = cell(length(tag),1);
    s = cell(length(tag),1); 
    for i=1:length(tag)
        tags{i} = tag(i):-2:0;
    end
    [s{:}] = ndgrid(tags{:});
    tags = cell2mat(cellfun(@(s)s(:),s(:)','un',0));
end

function [Afun,Alarge,As] = get_momentA(tags,n)
    tags = tags(:,1:n);
    X = str2sym(strcat('x',num2str((1:n)')));
    S = str2sym(strcat('s',num2str((1:n)')));
    syms x s Alarge
    m_fun = @(p) prod(p-1:-2:1)*s^p;
    As = cell(n,1);
    for nn=1:n
        syms A
        p = max(tags(:,nn));
        for i=1:p+1
            for j=i:2:p+1
                A(i,j) = nchoosek(j-1,i-1)*m_fun(j-i)*x^(i-1);
            end
        end
        As{nn} = subs(A,[x,s],[X(nn),S(nn)]);
    end    
    
    for tt= 1:size(tags,1)
        tag = tags(tt,:);
        row_ind = find(ismember(tags,tag,'rows'));
        prod_tags = 1;
        for nn=1:n
            prod_tags=expand(prod_tags*sum(As{nn}(:,tag(nn)+1)));
        end
        [a,b]=coeffs(prod_tags,X);
        tag_reps = zeros(length(b),n);
        for k=1:length(b)
            for nn=1:n
                tag_reps(k,nn) = subs(diff(b(k),1,X(nn)),X,ones(n,1));
            end
            col_ind = find(ismember(tags,tag_reps(k,:),'rows'));
            Alarge(col_ind,row_ind) = a(k);
        end
    end
    Afun = @(sigs) eval(subs(Alarge,S,sigs));
end

function [lib,true_S] = trig_poly_lib(polys,trigs,x_diffs,nstates,ndims,numeq,custom_remove_f,custom_remove_t,true_nz_weights)
    tags_poly = get_tags(polys,[],nstates);
    tags_trig = get_tags([],trigs,nstates);
    tags = prodtags(tags_poly,tags_trig);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    if ~isempty(diff_tags)
        diff_tags = diff_tags(diff_tags(:,end)==0,:);
        true_S = cell(numeq,1);
        for i=1:size(diff_tags,1)
            for j=1:length(tags)
                if all([~and(sum(diff_tags(i,:))>0,...
                        isequal(tags{j},zeros(1,nstates))),...
                        ~cellfun(@(b)b([tags{j} diff_tags(i,:)]),custom_remove_f),...
                        ~ismember_rows([tags{j} diff_tags(i,:)],custom_remove_t)])
                    lib.add_terms(term('ftag',tags{j},'linOp',diff_tags(i,:)));


                    %%% get true weights if specified
                    if exist('true_nz_weights','var')
                        if ~isempty(true_nz_weights)
                            S = cellfun(@(tnz) ismember(tnz(:,1:end-1),[tags{j} diff_tags(i,:)],'rows'), true_nz_weights(:), 'Un',0);
                            for jj=1:length(S)
                                if any(S{jj})
                                    true_S{jj} = [true_S{jj};[length(lib.terms) true_nz_weights{jj}(S{jj},end)]];
                                end
                            end
                        end
                    end

                end
            end
        end
    else
        true_S = [];
    end
end