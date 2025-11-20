classdef wendy_model_qc < wendy_model
    properties
        Cfacs
    end

    methods
        function obj = wendy_model_qc(dat,lib,tf,statcorrect,varargin)
            obj = obj@wendy_model(dat,lib,tf,statcorrect,varargin{:});
            obj.Cfacs = [];
        end
    end

    methods

        function obj = get_Cfacs(obj)
            R0 = obj.get_R0;
            obj.Cfacs = cell(obj.ntraj,obj.numeq,obj.numeq);
            for nt=1:obj.ntraj
                for n=1:obj.numeq
                    n1=size(obj.Gs{nt}{n},1);
                    n2=size(obj.Gs{nt}{n},2)+1;
                    for m=1:n
                        m1=size(obj.Gs{nt}{m},1);
                        m2=size(obj.Gs{nt}{m},2)+1;
                        obj.Cfacs{nt,n,m} = repmat( {sparse(n1,m1)}, n2,m2);
                    end
                end
            end

            for nt=1:obj.ntraj
                for n=1:obj.numeq
                    for j=0:size(obj.Gs{nt}{n},2)
                        if j==0
                            term_temp1 = obj.lhsterms{n};
                        else
                            term_temp1 = obj.lib(n).terms{j};
                        end
                        Ltilde = obj.get_Ltilde(term_temp1,nt,n);
                        Ltilde = Ltilde*R0;
                        for m=1:n
                            if m==n
                                jj_max = j;
                            else
                                jj_max = size(obj.Gs{nt}{m},2);
                            end
                            for jj=0:jj_max
                                if jj==0
                                    term_temp2 = obj.lhsterms{m};
                                else
                                    term_temp2 = obj.lib(m).terms{jj};
                                end
                                if any( abs(term_temp1.ftag) & abs(term_temp2.ftag) ) % assuming noise is independent between compartments (R0 is block diagonal)
                                    Ltilde2 = obj.get_Ltilde(term_temp2,nt,m);
                                    obj.Cfacs{nt,n,m}{j+1,jj+1} = Ltilde*Ltilde2';
                                end
                            end
                        end
                    end
                end
            end
        end

        function L = get_Ltilde(obj,tt,nt,m)
            A = tt.diffmat(obj.dat(nt));
            A = cell2mat(A(:)');
            V = obj.tf{nt}{m}.get_testmat(tt.linOp);
            L = V*A;
            % A = tt.evalgrads(obj.dat(nt));
            % V = obj.tf{nt}{m}.get_testmat(tt.linOp);
            % A = cellfun(@(g)V.*g(:)',A,'un',0);
            % L = cell2mat(A(:)');
        end

        function obj = get_cov(obj,w)
            if exist('w','var')
                obj.weights = w;
                if isempty(obj.Cfacs)
                    obj.get_Cfacs;
                end
                toggle_compute = 1;
            elseif and(~isempty(obj.weights),isempty(obj.cov))
                if isempty(obj.Cfacs)
                    obj.get_Cfacs;
                end
                toggle_compute = 1;
            else
                toggle_compute = 0;
            end
            if toggle_compute
                if any(obj.weights)
                    w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib)); 
                    w = cellfun(@(w)[-1;w],w,'un',0);
                    C = [];
                    for nt=1:obj.ntraj
                        Cnt = [];
                        N = sum(cellfun(@(G)size(G,1),obj.Gs{nt}));
                        for n=1:obj.numeq
                            p = 0;             % current write column (0-based)
                            for m=1:n
                                Cm = 0;
                                for j=1:size(obj.Gs{nt}{n},2)+1
                                    if m==n
                                        jj_max = j;
                                    else
                                        jj_max = size(obj.Gs{nt}{m},2)+1;
                                    end
                                    for jj=1:jj_max
                                        Cm_new = (w{n}(j)*w{m}(jj))*obj.Cfacs{nt,n,m}{j,jj};
                                        Cm = sparse(Cm + Cm_new);
                                        if and(m==n,j~=jj)
                                            Cm = sparse(Cm + Cm_new');
                                        end
                                    end
                                end           
                                if m == 1
                                    m_row = size(Cm,1);
                                    Cn = sparse(m_row, N);  % pre-size to final width (zeros are free)
                                end
                                k = size(Cm,2);         % width of this block
                                Cn(:, p + (1:k)) = Cm;   % write into the slice
                                p = p + k;
                            end
                            Cnt = cat(1,Cnt,Cn);
                        end
                        Cnt = tril(Cnt);
                        Cnt = Cnt + Cnt.'-spdiags(diag(Cnt),0,size(Cnt,1),size(Cnt,2));
                        C = blkdiag(C,Cnt);
                    end
                    obj.cov = C;
                else
                    obj.cov = speye(sum(cellfun(@(LS) sum(cellfun(@(L)size(L,1),LS)), obj.L0)));
                end
            end
        end

        function gC = get_gradC(obj)
            w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib));
            w = cellfun(@(w)[-1;w],w,'un',0);
            gC = cellfun(@(w)repmat({sparse(size(obj.cov,1),size(obj.cov,1))},length(w)-1,1),w,'un',0);
            for nt = 1:obj.ntraj
                p=0; % equation row block 
                for n=1:obj.numeq
                    for j=1:length(w{n})-1 % gradients with respect to RHS coefficients, exclude LHS
                        q=0; % equation column block
                        for m=1:obj.numeq % get blocks for each matrix
                            if m<=n % get block as normal
                                Cm = 0;
                                for jj=1:size(obj.Gs{nt}{m},2)+1
                                    if m<n
                                        Cm_new = w{m}(jj)*obj.Cfacs{nt,n,m}{j+1,jj};
                                    else
                                        if jj<=j+1
                                            Cm_new = w{m}(jj)*obj.Cfacs{nt,n,m}{j+1,jj};
                                        else
                                            Cm_new = w{m}(jj)*obj.Cfacs{nt,n,n}{jj,j+1}.';
                                        end
                                    end
                                    Cm = sparse(Cm + Cm_new);
                                end
                            else % get transpose of block below diag
                                Cm = 0;
                                for jj=1:size(obj.Gs{nt}{m},2)+1
                                    Cm_new = w{m}(jj)*obj.Cfacs{nt,m,n}{jj,j+1}.';
                                    Cm = sparse(Cm + Cm_new);
                                end
                            end
                            if m == 1
                                m_row = size(Cm,1);
                            end
                            k = size(Cm,2);         % width of this block
                            gC{n}{j}(p+(1:m_row), q + (1:k)) = Cm;   % write into the slice
                            q = q + k;
                        end
                    end
                    p = p + m_row;
                end
            end
            gC = cellfun(@(ggC)cellfun(@(g) g+g.',ggC,'un',0),gC,'un',0);
            gC = vertcat(gC{:})';
        end

    end
end


