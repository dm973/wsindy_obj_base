%%% defines a library object simply used to gather a set of terms. Operations
%%% consist mainly of mapping term methods over a list of terms, with some exceptions
%%% (such as trimlib) which performs a library trimming based on taylor expansion of the terms

classdef library < handle

    properties
        terms
        tags
        nstates
    end

    methods
        
        function obj = library(varargin)
            p = inputParser;
            addParameter(p,'tags',{});
            addParameter(p,'nstates',1);
            addParameter(p,'terms',{});
            addParameter(p,'polys',[]);
            addParameter(p,'gradon',1);
            addParameter(p,'neg',0);
            addParameter(p,'boolT',@(tags)true(size(tags,1),1));
            parse(p,varargin{:})
            obj.tags = {};
            obj.terms = {};
            if ~isempty(p.Results.terms)
                obj = add_terms(obj,p.Results.terms,p.Results.gradon);
            end

            if ~isempty(p.Results.tags)
                obj = add_terms(obj,p.Results.tags,p.Results.gradon);
            end

            if ~isempty(obj.terms)
                obj.nstates = obj.terms{1}.nstates;
            else
                obj.nstates = p.Results.nstates;
            end

            if ~isempty(p.Results.polys)
                tags_poly = get_tags(p.Results.polys,[],obj.nstates,'neg',p.Results.neg,'boolT',p.Results.boolT);
                obj = add_terms(obj,tags_poly,p.Results.gradon);
            end

        end

    end

    methods
        function fs = get_fHandles(obj,inds)
            fs = cellfun(@(t)t.fHandle,obj.terms(inds),'uni',0);
        end

        function obj = add_terms(obj,tags,gradon,dup)
            if ~exist('gradon','var')
                gradon = 1;
            end
            if ~exist('dup','var')
                dup = 0;
            end
            if isequal(class(tags),'double') % if collection of pw terms, then add
                tags = mat2cell(tags,ones(size(tags,1),1),size(tags,2));
            elseif isequal(class(tags),'function_handle')
                tags = term('fHandle',tags);
            end
            if isequal(class(tags),'cell')
                for j=1:length(tags)
                    if or(dup,~any(cellfun(@(t)isequal(t,tags{j}),obj.tags)))
                        if isequal(class(tags{j}),'double')
                            obj.terms = [obj.terms,{term('ftag',tags{j},'gradon',gradon)}];
                            obj.tags = [obj.tags,tags(j)];
                    elseif isequal(class(tags{j}),'function_handle')
                            obj.terms = [obj.terms,{term('fHandle',tags{j},'gradon',gradon)}];
                            obj.tags = [obj.tags,tags(j)];
                        elseif any(cellfun(@(x) isequal(x,'absterm'),superclasses(tags{j})))
                            obj.terms = [obj.terms,{tags{j}}];
                            obj.tags = [obj.tags,tags(j)];
                        end
                    end
                end
            elseif any(cellfun(@(x) isequal(x,'absterm'),superclasses(tags)))
                for j=1:length(tags)
                    % try
                    %     check = functions(tags(j).fHandle).function;
                    % catch
                    %     check = NaN;
                    % end
                    % if or(dup,~any(cellfun(@(t)isequal(functions(t.fHandle).function,check),obj.terms)))
                        obj.terms = [obj.terms,{tags(j)}];
                        obj.tags = [obj.tags,{tags(j).ftag}];
                    % end
                end
            end
            obj.nstates = obj.terms{1}.nstates;
        end

        function Theta = evalterms(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta = cell2mat(cellfun(@(tm) tm.evalterm(dat),obj.terms(S),'uni',0));
        end

        function Theta_grad = evalGradterms(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta_grad = cellfun(@(tm) tm.evalgrads(dat),obj.terms(S),'uni',0);
        end

        function Theta_grad = grad2mat(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta_grad = obj.evalGradterms(dat,S);
            Theta_grad = cellfun(@(tm) reshape([tm{:}],[],1), Theta_grad,'uni',0);
            Theta_grad = cell2mat(Theta_grad);
        end

        function Theta_H = evalHvec(obj,dat,J)
            Theta_grad = grad2mat(obj,dat);
            Theta_H = J*Theta_grad;
        end

        function obj = complib(obj,t2s,t1s)
            for j=1:length(t2s)
                for i=1:length(t1s)
                    obj.add_terms(compterm(t2s{j},t1s{i}));
                end
            end
        end

        function obj = trimlib(obj,pt,ord,tol,toggleH,togglepar)
            if ~exist('toggleH','var')
                toggleH=0;
            end
%             q = sym('q');
            args = num2cell(pt);
            tays = cell(1,length(obj.tags));
            ns = obj.nstates;
            termz = obj.terms;
            parfor(j=1:length(obj.tags),togglepar)
                tays{j} = termz{j}.get_taylor(pt,ord,toggleH);
            end
            tays = cell2mat(tays);
            sp = tays~=0;
            Xsparse = unique(sp','rows')';
            ii = cell(size(Xsparse,2),1);
            As = {};
            if toggleH==1
                ord = ord -1;
            end
            for j=1:size(Xsparse,2)
                ii = find(~any(sp-Xsparse(:,j),1));
                if length(ii)>1
                    G = tays(:,ii);
                    for l=1:obj.nstates
                        k = find(Xsparse((l-1)*ord+(1:ord),j),1);
                        if ~isempty(k)
                            G((l-1)*ord+(1:ord),:) = G((l-1)*ord+(1:ord),:)./G((l-1)*ord+k,:);
                        end
                    end
                    %%% compute ||G_i-G_j||, note the first non-zero entries already agree
                    A = real(sqrt((vecnorm(G)').^2 - 2*(G')*G + (vecnorm(G).^2)));
                    A = A-triu(A);
                    S = ones(size(A));
                    S = logical(S - triu(S));
                    A = and(A < tol,S);
                    if any(A(:))
                        As = [As,{{Xsparse(:,j),ii,G,A}}];
                    end
                end
            end
            
            rm_inds = [];
            for i=1:length(As)
                max_inds = [];
                for j=1:obj.nstates
                    foo_ind = 1+ord-find(flipud(As{i}{1}((j-1)*ord+1:j*ord)),1);
                    max_inds = [max_inds foo_ind];
                end
                while any(As{i}{end}(:))
                    r={[],[]};
                    [r{:}] = find(As{i}{end});
                    r = [r{:}];
                    r = unique(r(:));
                    [~,a] = max(sum(abs(As{i}{3}(max_inds,r))));
                    rm = r(a);
                    rm_inds = [rm_inds As{i}{2}(rm)];
                    As{i}{2} = As{i}{2}([1:rm-1 rm+1:end]);
                    As{i}{3} = As{i}{3}(:,[1:rm-1 rm+1:end]);
                    As{i}{4} = As{i}{4}(:,[1:rm-1 rm+1:end]);
                    As{i}{4} = As{i}{4}([1:rm-1 rm+1:end],:);
                end
            end
            
            J = length(obj.terms);
            obj.terms = obj.terms(setdiff(1:J,rm_inds)); 
            obj.tags = obj.tags(setdiff(1:J,rm_inds)); 
        end

    end

end