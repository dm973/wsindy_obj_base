classdef wendy_model < wsindy_model
    properties
        Hfac
        H
        biasfac
        bias
    end

    methods
        function obj = wendy_model(dat,lib,tf,varargin)
            % written as t2(t1(x))
            obj = obj@wsindy_model(dat,lib,tf,varargin{:});
            obj.Hfac = {};
            obj.H = [];
            obj.biasfac = {};
            obj.bias = [];
        end
    end

    methods


       % build covariance 
        function obj = get_Hfac(obj)
            % disp(['getting covariance factors...'])
            if ~exist('coarsen','var')
                coarsen = 1;
            end

            if isempty(obj.L0)
                obj.L0 = repmat({cell(obj.numeq,1)},obj.ntraj,1);
                if obj.toggleH
                    obj.L1 = repmat({repmat(arrayfun(@(L)cell(length(L.terms),1),obj.lib,'uni',0),obj.nstates,1)},obj.ntraj,1);
                else
                    obj.L1 = repmat({arrayfun(@(L)cell(length(L.terms),1),obj.lib,'uni',0)},obj.ntraj,1);
                end
            end
            S = obj.get_supp;
            if isempty(obj.features)
                obj.get_features;
            end
            for j=1:obj.ntraj
                for i=1:obj.numeq
                    if isempty(obj.L0{j}{i})
                        grads = obj.lhsterms{i}.diffmat(obj.dat(j));
                        grads = cell2mat(grads(:)');
                        V = obj.tf{j}{i}.get_testmat(obj.lhsterms{i}.linOp);
                        obj.L0{j}{i} = V*grads(:,1:obj.coarsen_L:end);
                    end
                    if obj.toggleH
                        for k=1:length(obj.lib.terms)
                            if and(ismember(k,S{1}),isempty(obj.L1{j}{i}{k}))
                                grads = obj.features{i}{ismember(S{1},k)}.diffmat(obj.dat(j));
                                V = obj.tf{j}{i}.get_testmat(obj.features{i}{ismember(S{1},k)}.linOp);
                                A = cell2mat(grads(:)');
                                obj.L1{j}{i}{k} = V*A(:,1:obj.coarsen_L:end);
                            end
                        end
                    else
                        for k=1:length(obj.lib(i).terms)
                            if and(ismember(k,S{i}),isempty(obj.L1{j}{i}{k}))
                                grads = obj.lib(i).terms{k}.diffmat(obj.dat(j));
                                V = obj.tf{j}{i}.get_testmat(obj.lib(i).terms{k}.linOp);
                                A = cell2mat(grads(:)');
                                obj.L1{j}{i}{k} = V*A(:,1:obj.coarsen_L:end);
                            end
                        end
                    end
                end 
            end
            % disp(['completed.'])
        end

        function obj = get_H(obj)
            w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib)); 
            S = obj.get_supp;
            if obj.toggleH
                S = repmat(S,1,obj.nstates);
                w = repmat(w,1,obj.nstates);
            end
            obj.L = obj.L0;
            for i=1:obj.ntraj
                for j=1:obj.numeq
                    for k=1:length(S{j})
                        obj.L{i}{j} = obj.L{i}{j} - w{j}(S{j}(k))*obj.L1{i}{j}{S{j}(k)};
                    end
                end
            end
            obj.L = cellfun(@(L)cell2mat(L),obj.L,'un',0);
            if obj.multitest == 1
                obj.L = cell2mat(obj.L);
            else
                obj.L = blkdiag(obj.L{:});
            end
        end

        function obj = get_biasfac(obj)
            if isempty(obj.biasfac)
                obj.biasfac = repmat({arrayfun(@(L)cell(length(L.terms),1),obj.lib,'uni',0)},obj.ntraj,1);
            end

            S = obj.get_supp;
            for j=1:obj.ntraj
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
            addParameter(IP,'toggle_cov',1);
            parse(IP,w,varargin{:});
            w = IP.Results.w;
            toggle_cov = IP.Results.toggle_cov;
            if toggle_cov==1
                obj.get_cov(w);
                obj.get_bias;
            else
                obj.weights = w;
                obj.get_features;
            end
        end

        function obj = get_cov(obj,w)
            if obj.multitest == 1
                obj.dat(1).get_R0;
                R0 = obj.dat(1).R0(1:obj.coarsen_L:end,1:obj.coarsen_L:end);
            else
                R0 = [];
                for i=1:obj.ntraj
                    obj.dat(i).get_R0;
                    R0 = blkdiag(R0,obj.dat(i).R0(1:obj.coarsen_L:end,1:obj.coarsen_L:end));
                end
            end
            if exist('w','var')
                s_old = obj.get_supp;
                obj.weights = w;
                s_new = obj.get_supp;
                if isempty(obj.L1)
                    obj.get_Lfac;
                else
                    if ~isempty(s_old)
                        if ~all(cellfun(@(so,sn)all(ismember(sn,so)),s_old,s_new))                     
                            obj.get_Lfac;
                        end
                    else
                        obj.get_Lfac;    
                    end
                end
                toggle_compute = 1;
            elseif and(~isempty(obj.weights),isempty(obj.cov))
                obj.get_Lfac;
                toggle_compute = 1;
            else
                toggle_compute = 0;
            end
            if toggle_compute
                if any(obj.weights)
                    obj.get_L;
                    obj.cov = (obj.L*R0)*(obj.L');
                else
                    obj.cov = speye(sum(cellfun(@(LS) sum(cellfun(@(L)size(L,1),LS)), obj.L0)));
                end
            end
        end

        function [G,b,RT] = apply_cov(obj,G,b,diag_reg)
            if isempty(obj.cov)
                obj.get_cov;
            end
            if isempty(obj.bias)
                obj.get_bias;
            end
            check = 0;
            while check == 0
                try
                    RT = obj.cov + (diag_reg/(1-diag_reg)*mean(diag(obj.cov)))*speye(size(obj.cov,1));
                    RT = chol( RT )';
                    RT = sqrt(1-diag_reg)*RT;
                    check = 1;
                catch
                    diag_reg = diag_reg*10;
                    fprintf('\nincreasing Cov regularization to %0.2e\n',diag_reg)
                end
            end
            G = RT \ G;
            b = RT \ (b - obj.bias);
        end
    end
end