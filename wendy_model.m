classdef wendy_model < wsindy_model
    properties
        Hfac
        H
        biasfac
        bias
        statcorrect
    end

    methods
        function obj = wendy_model(dat,lib,tf,statcorrect,varargin)
            obj = obj@wsindy_model(dat,lib,tf,varargin{:});
            obj.Hfac = {};
            obj.H = [];
            obj.biasfac = {};
            obj.bias = [];
            obj.statcorrect = statcorrect;
        end
    end

    methods

        function obj = get_Hfac(obj)
            % disp(['getting covariance factors...'])
            obj.Hfac = repmat({arrayfun(@(L)cell(length(L.terms),1),obj.lib,'uni',0)},obj.ntraj,1);
            S = obj.get_supp;
            for j=1:obj.ntraj
                for i=1:obj.numeq
                    for k=1:length(obj.lib(i).terms)
                        if and(ismember(k,S{i}),isempty(obj.Hfac{j}{i}{k}))
                            hess = obj.lib(i).terms{k}.get_hess;
                            if any(arrayfun(@(tm)tm.coeff~=0,hess))
                                hess = arrayfun(@(tm)...
                                    spdiags(reshape(tm.evalterm(obj.dat(j)),[],1),0,prod(obj.dat(j).dims),prod(obj.dat(j).dims)),...
                                    hess,'uni',0);
                                V = obj.tf{j}{i}.get_testmat(obj.lib(i).terms{k}.linOp);
                                A = cell2mat(hess(:)');
                                obj.Hfac{j}{i}{k} = V*A;
                            else
                                obj.Hfac{j}{i}{k} = sparse(length(obj.bs{1}{i}),prod(obj.dat(j).dims)*obj.nstates^2);
                            end
                        end
                    end
                end 
            end
            % disp(['completed.'])
        end

        function [R,Rsig] = get_H_R(obj)
            N = obj.nstates;
            R = sparse(0,0);
            for n=1:obj.ntraj
                sigmas = obj.dat(n).estimate_sigma;
                Rsig = zeros(N^2,N^2);
                for i=1:N % first component column index
                    for j=1:N % second component column index
                        for k=1:N % first component row index
                            for l=1:N % second component row index
                                if all([i==k,j==l,i~=j])
                                    Rsig((i-1)*N+j,(k-1)*N+l) = sigmas{i}^2*sigmas{j}^2;
                                elseif all([i==j,j==k,k==l])
                                    Rsig((i-1)*N+j,(k-1)*N+l) = 2*sigmas{i}^4;
                                elseif all([i==l,j==k,i~=j])
                                    Rsig((i-1)*N+j,(k-1)*N+l) = sigmas{i}^2*sigmas{j}^2;
                                end
                            end
                        end
                    end
                end
                R = blkdiag(R,kron(Rsig,speye(prod(obj.dat(n).dims))));
            end
        end

        function obj = get_H(obj)
            obj.get_Hfac;
            w = obj.reshape_w; 
            S = obj.get_supp;
            obj.H = arrayfun(@(s)cell(obj.numeq,1),(1:obj.ntraj)','uni',0);
            for i=1:obj.ntraj
                for j=1:obj.numeq
                    obj.H{i}{j} = w{j}(S{j}(1))*obj.Hfac{i}{j}{S{j}(1)};
                    for k=2:length(S{j})
                        obj.H{i}{j} = obj.H{i}{j} + w{j}(S{j}(k))*obj.Hfac{i}{j}{S{j}(k)};
                    end
                end
            end
            obj.H = cellfun(@(L)cell2mat(L),obj.H,'un',0);
            obj.H = blkdiag(obj.H{:});
            obj.H = 1/4*(obj.H*(obj.get_H_R*obj.H'));
        end

        function obj = get_biasfac(obj)
            if isempty(obj.biasfac)
                obj.biasfac = repmat({arrayfun(@(L)cell(length(L.terms),1),obj.lib,'uni',0)},obj.ntraj,1);
            end

            S = obj.get_supp;
            for j=1:obj.ntraj
                obj.dat(j).get_R0;
                for i=1:obj.numeq
                    for k=1:length(obj.lib(i).terms)
                        if and(ismember(k,S{i}),isempty(obj.biasfac{j}{i}{k}))
                            lap = obj.lib(i).terms{k}.get_lap;
                            Y = arrayfun(@(d2) d2.evalterm(obj.dat(j)),lap(:),'uni',0);
                            Y = cellfun(@(d2) d2(:),Y,'uni',0);
                            Y = cell2mat(Y);
                            Y = obj.dat(j).R0*Y;
                            x = prod(obj.dat(j).dims);
                            msk = spdiags(ones(x,obj.nstates),0:x:x*(obj.nstates-1),x,length(Y));
                            V = obj.tf{j}{i}.get_testmat(obj.lib(i).terms{k}.linOp);
                            obj.biasfac{j}{i}{k} = V*(msk*Y);
                        end
                    end
                end 
            end
            % disp(['completed.'])
        end

        function obj = get_bias(obj)
            if isempty(obj.biasfac)
                obj.get_biasfac;
            end
            w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib)); 
            obj.bias = cellfun(@(bc)cellfun(@(b)b*0,bc,'uni',0),obj.bs,'uni',0);
            S = obj.get_supp;
            for i=1:obj.ntraj
                for j=1:obj.numeq
                    for k=1:length(S{j})
                        obj.bias{i}{j} = obj.bias{i}{j} + w{j}(S{j}(k))*obj.biasfac{i}{j}{S{j}(k)};
                    end
                end
            end
            obj.bias = cellfun(@(b)cell2mat(b),obj.bias,'un',0);
            obj.bias = -1/2*cell2mat(obj.bias(:));
        end

        function obj = add_weights(obj,w,varargin)
            IP = inputParser;
            addRequired(IP,'w');
            addParameter(IP,'toggle_cov',0);
            parse(IP,w,varargin{:});
            w = IP.Results.w;
            toggle_cov = IP.Results.toggle_cov;
            if toggle_cov==1
                obj.get_cov(w);
                if obj.statcorrect(2)>0
                    obj.get_bias;
                else
                    obj.bias = sparse(length(obj.b{1}),1);
                end
            else
                obj.weights = w;
                obj.get_features;
            end
        end

        function [G,b,RT] = apply_cov(obj,G,b,diag_reg)
            if obj.statcorrect(2)==1
                if isempty(obj.bias)
                    obj.get_bias;
                end
                b = b - obj.bias;
            end

            if obj.statcorrect(1)>0
                if isempty(obj.cov)
                    obj.get_cov;
                end
                RT = obj.cov;
                if obj.statcorrect(1)>1
                    if obj.dat(1).estimate_sigma{1}>obj.statcorrect(1)-1
                        obj.get_H;
                        RT = RT + obj.H;
                    end
                end
                check = 0;
                while check == 0
                    try
                        RT = RT + (diag_reg/(1-diag_reg)*mean(diag(RT)))*speye(size(obj.cov,1));
                        RT = chol( RT )';
                        RT = sqrt(1-diag_reg)*RT;
                        check = 1;
                    catch
                        if diag_reg==0
                            diag_reg=10^-16;
                        else
                            diag_reg = diag_reg*10;
                        end
                        fprintf('\nincreasing Cov regularization to %0.2e\n',diag_reg)
                    end
                end
                G = RT \ G;
                b = RT \ b;
            else
                RT = speye(length(b),length(b))*norm(obj.res)/sqrt(length(b)-1);
            end
        end
   
        function A = apply_gradC(obj,v)
            R0 = obj.get_R0;
            w = R0*(obj.L'*v);
            v_cell = obj.vec2cell(v,cellfun(@(b)length(b),obj.bs{1}));
            S = obj.get_supp;
            A = obj.G{1}*0;
            indc = 1;
            indr = 1;
            for j=1:obj.numeq
                for k=1:length(S{j})
                    A(:,indc) = obj.L*(R0*(obj.L1{1}{j}{S{j}(k)}'*v_cell{j}));
                    foo = obj.L1{1}{j}{S{j}(k)}*w;
                    A(indr:indr+length(foo)-1, indc) = A(indr:indr+length(foo)-1, indc) + foo;
                    indc = indc+1;
                end
                indr = indr + length(v_cell{j});
            end
        end

        function w_cell = vec2cell(obj,w,length_vec)
            ind=1;
            numeq = length(length_vec);
            w_cell = cell(numeq,1);
            for nn=1:numeq
                w_cell{nn} = w(ind:ind+length_vec(nn)-1);
                ind = ind+length_vec(nn);
            end
        end


    end
end