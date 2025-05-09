classdef wsindy_model < handle

    properties
        dat
        ndims
        ntraj
        nstates
        numeq
        lhsterms
        lib
        Ji
        toggleH
        tf
        Gs
        bs
        G
        b
        weights
        cov
        L
        L0
        L1
        coarsen_L
        subinds
        features
        diag_reg
        catm
        nzs
        cv
        multitest
    end

    methods

        function obj = wsindy_model(dat,lib,tf,varargin)
            obj.dat = dat;
            obj.lib = lib(:);
            obj.ndims = obj.dat(1).ndims;
            obj.ntraj = length(obj.dat);
            obj.nstates = obj.dat(1).nstates;
            
            E = eye(obj.nstates);
            default_toggleH = 0;
            default_Ji = [];
            default_diag_reg = 10^-10;
            default_toggle_Gb = 1;
            default_catm = 'component';
            default_cv = [];
            default_coarsen_L = 1;

            p = inputParser;
            addRequired(p,'dat');
            addRequired(p,'lib');
            addRequired(p,'tf');
            addParameter(p,'lhsterms',[]);
            addParameter(p,'toggleH',default_toggleH);
            addParameter(p,'Ji',default_Ji);
            addParameter(p,'diag_reg',default_diag_reg);
            addParameter(p,'opt',default_diag_reg);
            addParameter(p,'toggle_Gb',default_toggle_Gb);
            addParameter(p,'catm',default_catm);
            addParameter(p,'cv',default_cv);
            addParameter(p,'coarsen_L',default_coarsen_L);
            addParameter(p,'multitest',0);

            parse(p,dat,lib,tf,varargin{:})

            obj.lhsterms = p.Results.lhsterms;
            if isempty(obj.lhsterms)
                obj.lhsterms = arrayfun(@(i)term('ftag',E(i,:),'linOp',[zeros(1,obj.ndims-1) 1]),(1:obj.nstates)','uni',0);
            elseif isequal(class(obj.lhsterms),'double')
                obj.lhsterms = arrayfunvec(obj.lhsterms,@(L)term('ftag',L(1:obj.nstates),'linOp',L(obj.nstates+1:end)),2,0);
            elseif isequal(class(obj.lhsterms),'diffOp')
                obj.lhsterms = num2cell(obj.lhsterms);
            elseif any(cellfun(@(x) isequal(x,'absterm'),superclasses(obj.lhsterms)))
                obj.lhsterms = num2cell(obj.lhsterms);
            end

            obj.numeq = length(obj.lhsterms);
            obj.toggleH = p.Results.toggleH;
            obj.Ji = p.Results.Ji;
            obj.diag_reg = p.Results.diag_reg;
            obj.catm = p.Results.catm;
            obj.cv = p.Results.cv;
            obj.coarsen_L = p.Results.coarsen_L;
            obj.multitest = p.Results.multitest;
            
            if and(obj.toggleH,isempty(obj.Ji))
                obj.Ji = kron(eye(obj.nstates/2),[[0 1];[-1 0]]);
            end

            obj.Gs = {};
            obj.bs = {};
            obj.cov = [];
            obj.L0 = {};
            obj.L1 = {};

            obj.features = {};
            obj.weights = [];

            if and(length(obj.lib)~=obj.numeq,~obj.toggleH)
                obj.lib = repelem(obj.lib,1,obj.numeq);
            end

            if or(isequal(class(tf),'testfcn'),ismember('testfcn',superclasses(tf)))
                obj.tf = repelem({repelem({tf},1,obj.numeq)},obj.ntraj,1);
            elseif isequal(class(tf),'cell')
                if or(isequal(class(tf{1}),'testfcn'),ismember('testfcn',superclasses(tf{1})))
                    if length(tf)==obj.numeq
                        obj.tf = repelem({tf},obj.ntraj,1);
                    elseif length(tf)==obj.ntraj
                        obj.tf = cellfun(@(tf)repelem({tf},1,obj.numeq),tf,'uni',0);
                    end
                else
                    obj.tf = tf;
                end
            else
                obj.tf = [];
            end
            % obj.tf = cellfun(@(t) t(:)', obj.tf(:), 'uni',0);

            if ~isempty(obj.cv)
                for j=1:obj.numeq
                    lib_foo = library();
                    lib_foo.complib(obj.lib(j).terms,obj.cv);
                    obj.lib(j) = lib_foo;
                    obj.lhsterms{j} = compterm(obj.lhsterms{j},obj.cv);
                end
            end

            if p.Results.toggle_Gb
                obj.get_Gb;
            end
            
        end

    end

    methods

        % build linear system 
        function Theta = get_theta(obj)
            if obj.toggleH
                Theta = cell(obj.ntraj,1);
                for j=1:obj.ntraj
                    Theta_temp = obj.lib.grad2mat(obj.dat(j).Uobs);
                    Ji_kron = kron(obj.Ji,speye(size(Theta_temp,1)/obj.nstates));
                    Theta{j} = Ji_kron*Theta_temp;
                end
            end
        end

        function obj = get_Gb(obj)
            if and(~isempty(obj.lib),~isempty(obj.tf))
                obj.Gs = cell(obj.ntraj,1);
                obj.bs = cell(obj.ntraj,1);
                for j=1:obj.ntraj
                    obj.bs{j} = cellfun(@(tf,gt) tf.test(obj.dat(j),gt), obj.tf{j}(:), obj.lhsterms(:),'uni',0);
                    if obj.toggleH
                        Theta = obj.lib.grad2mat(obj.dat(j).Uobs);
                        Ji_kron = kron(obj.Ji,speye(size(Theta,1)/obj.nstates));
                        Theta = Ji_kron*Theta;
                        J = length(obj.lib.terms);
                        obj.Gs{j} = cellfun(@(Th,v) conv2(v,1,Th,'valid'), ...
                            mat2cell(Theta,ones(obj.nstates,1)*size(Theta,1)/obj.nstates,J),...
                            cellfun(@(tf) tf.Cfs{1}(1,:)', obj.tf{j}(:),'uni',0),'uni',0);
                        obj.Gs{j} = arrayfun(@(i)obj.Gs{j}{i}(unique(min(obj.tf{j}{i}.subinds{1},obj.tf{j}{i}.dims(1)-2*obj.tf{j}{i}.rads(1))),:), 1:length(obj.Gs{j}),'uni',0)';
                    else
                        obj.Gs{j} = arrayfun(@(tf,lib) cell2mat(cellfun(@(tm) tf.test(obj.dat(j),tm),lib.terms,'uni',0)),[obj.tf{j}{:}]',obj.lib(:),'uni',0);
                    end
                end
            end
        end

        function obj = get_Gb_ij(obj,j,i)
            if isempty(obj.Gs)
                obj.get_Gb;
            else
                obj.bs{j}{i} = obj.tf{j}{i}.test(obj.dat(j),obj.lhsterms{i});
                if obj.toggleH
                    Theta = obj.lib.grad2mat(obj.dat(j).Uobs);
                    Ji_kron = kron(obj.Ji,speye(size(Theta,1)/obj.nstates));
                    Theta = Ji_kron*Theta;
                    J = length(obj.lib.terms);
                    Theta  = mat2cell(Theta,ones(obj.nstates,1)*size(Theta,1)/obj.nstates,J);
                    obj.Gs{j}{i} = conv2(obj.tf{j}{i}.Cfs{1}(1,:),1,Theta{i},'valid');
                    obj.Gs{j}{i} = obj.Gs{j}{i}(unique(min(obj.tf{j}{i}.subinds{1},obj.tf{j}{i}.dims(1)-2*obj.tf{j}{i}.rads(1))),:);
                else
                    obj.Gs{j}{i} = cell2mat(cellfun(@(tm) obj.tf{j}{i}.test(obj.dat(j),tm),obj.lib(i).terms,'uni',0));
                end
            end
        end

        function obj = cat_Gb(obj,varargin)

            default_cator = obj.catm;
            p = inputParser;
            addParameter(p,'cat',default_cator);
            parse(p,varargin{:})
            cator = p.Results.cat;

            if isempty(obj.Gs)
                obj.get_Gb;
            end

            if or(obj.toggleH,cator==-1)
                obj.G = cellfun(@(Gs) cell2mat(Gs), obj.Gs,'uni',0);
                obj.b = cellfun(@(bs) cell2mat(bs), obj.bs,'uni',0);
                obj.G = {cell2mat(obj.G)};
                obj.b = {cell2mat(obj.b)};
            elseif isequal(cator,'component')
                obj.catm = 'component';
                obj.b = repmat({[]},obj.numeq,1);
                obj.G = repmat({[]},obj.numeq,1);
                for n=1:obj.numeq
                    for m=1:obj.ntraj
                        obj.b{n} = [obj.b{n};obj.bs{m}{n}];
                        obj.G{n} = [obj.G{n};obj.Gs{m}{n}];
                    end
                end
            elseif isequal(cator,'blkdiag')
                obj.catm = 'blkdiag';
                obj.G = cellfun(@(Gs) blkdiag(Gs{:}), obj.Gs,'uni',0);
                obj.b = cellfun(@(bs) cell2mat(bs), obj.bs,'uni',0);
                obj.G = {cell2mat(obj.G)};
                obj.b = {cell2mat(obj.b)};
            end

        end

        % update / transform weights 
        function obj = add_weights(obj,w,varargin)
            IP = inputParser;
            addRequired(IP,'w');
            addParameter(IP,'toggle_cov',0);
            parse(IP,w,varargin{:});
            w = IP.Results.w;
            toggle_cov = IP.Results.toggle_cov;
            if toggle_cov==1
                obj.get_cov(w);
            else
                obj.weights = w;
                obj.get_features;
            end
        end

        function params = reshape_w(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            if ~isempty(IP.Results.w)
                params = cellfun(@(w)w,reshape_cell(IP.Results.w,arrayfun(@(L)length(L.terms),obj.lib)), 'uni',0);
            else
                params = {};
            end
        end

        function params = get_params(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            params = obj.reshape_w('w',IP.Results.w);
            if ~isempty(params)
                params = cellfun(@(w)w(w~=0),params,'uni',0);
            end
        end

        function supp = get_supp(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            params = obj.reshape_w('w',IP.Results.w);
            if ~isempty(params)
                supp = cellfun(@(w)find(w~=0),params, 'uni',0);
            else
                supp = {};
            end
        end

        % get / display model features 
        function obj = get_features(obj)
            obj.features = cell(obj.numeq,1);
            % disp(['getting model features...'])
            if obj.toggleH
                perm = obj.Ji*(1:obj.nstates)';
                inds = find(obj.weights);
                for n=1:obj.nstates
                    ftemp = cell(1,length(inds));
                    for j=1:length(inds)
                        ftemp{j} = obj.lib.terms{inds(j)}.gradterms(abs(perm(n)));
                        ftemp{j}.coeff = sign(perm(n));
                        ftemp{j} = ftemp{j}.get_grads;
                    end
                    obj.features{n} = ftemp;
                end
            else
                if isvector(obj.weights)
                    w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib));
                end         
                for n=1:obj.numeq
                    obj.features{n} = obj.lib(n).terms(w{n}~=0);
                end
            end
            % disp(['completed.'])
        end

        function [rhs,params,params_diff,feats,feats_diff] = get_rhs(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            weights_cell = cellfun(@(w)w,reshape_cell(IP.Results.w,arrayfun(@(L)length(L.terms),obj.lib)), 'uni',0);
            params = cellfun(@(p)p(p~=0),weights_cell,'uni',0);
            supp = obj.get_supp('w',IP.Results.w);

            feats = cell(obj.numeq,1);
            for n=1:obj.numeq
                for j=1:length(supp{n})
                    if isempty(obj.cv)
                        feats{n}{j} = obj.lib(n).terms{supp{n}(j)}.get_rhs;
                    else
                        tt = term('ftag',obj.lib(n).terms{supp{n}(j)}.t2.ftag,'fHandle',obj.lib(n).terms{supp{n}(j)}.t2.fHandle,'linOp',obj.lib(n).terms{supp{n}(j)}.linOp);
                        feats{n}{j} = tt.get_rhs;
                    end
                end
            end
            rhs_f = @(x)rhs_fun(feats,params,x);
            t_diff_ords = cellfun(@(lhs)lhs.linOp.difftags(end),obj.lhsterms);
            if range(t_diff_ords)==0
                diff_ord = t_diff_ords(1);
                if obj.nstates == obj.numeq
                    rhs = @(x)[reshape(x(obj.nstates+1:diff_ord*obj.nstates),[],1);rhs_f(x)];
                elseif obj.nstates == 2*obj.numeq
                    rhs = @(x)[reshape(x(obj.nstates/2+1:obj.nstates),[],1);rhs_f(x)];
                else                    
                    rhs = rhs_f;
                end
            else
                disp(['cannot form rhs, diff orders not matching'])
            end
            % if all(cellfun(@(lhs)isequal(lhs.linOp.difftags,1),obj.lhsterms))
            %     rhs = @(x)rhs_fun(feats,params,x);
            % elseif all(cellfun(@(lhs)isequal(lhs.linOp.difftags,2),obj.lhsterms))
            %     rhs_f = @(x)rhs_fun(feats,params,x);
            %     if obj.nstates == obj.numeq
            %         rhs = @(x)[reshape(x(obj.nstates+1:2*obj.nstates),[],1);rhs_f(x)];
            %     elseif obj.nstates == 2*obj.numeq
            %         rhs = @(x)[reshape(x(obj.nstates/2+1:obj.nstates),[],1);rhs_f(x)];
            %     end
            % elseif all(cellfun(@(lhs)isequal(lhs.linOp.difftags,3),obj.lhsterms))
            %     rhs_f = @(x)rhs_fun(feats,params,x);
            %     if obj.nstates == obj.numeq
            %         rhs = @(x)[reshape(x(obj.nstates+1:3*obj.nstates),[],1);rhs_f(x)];
            %     end
            % end

            params_diff={}; feats_diff={};

        end

        function str_c = disp_mod(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            w = IP.Results.w;
            obj.get_features;
            w = obj.get_params('w',w);
            % w = obj.get_params;
            % w{:}
            str_c = cell(obj.numeq,1);
            for i=1:obj.numeq
                str={};
                for j=1:length(obj.features{i})
                    s = obj.features{i}{j}.get_str;
                    str=[str;{['(',num2str(w{i}(j),'%.2e'),')*',s]}];
                end
                str_c{i}=str;
            end

        end

        % build covariance 
        function obj = get_Lfac(obj,coarsen)
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
        end

        function obj = get_L(obj)
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

        function R0 = get_R0(obj)
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
        end

        function obj = get_cov(obj,w)
            R0 = obj.get_R0;
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
                    % if obj.toggleH
                    %     obj.L = cellfun(@(L)cell2mat(L),obj.L0,'uni',0);
                    %     obj.L = cellfun(@(x,y) x - y, obj.L,add_L1(obj),'uni',0);
                    %     if obj.multitest == 1
                    %         obj.L = cell2mat(obj.L);
                    %     else
                    %         obj.L = blkdiag(obj.L{:});
                    %     end
                    % else
                        obj.get_L;
                        obj.cov = obj.L*(R0*obj.L');
                    % end
                else
                    obj.cov = speye(sum(cellfun(@(LS) sum(cellfun(@(L)size(L,1),LS)), obj.L0)));
                end
            end
        end

        function [G,b,RT] = apply_cov(obj,G,b,diag_reg)
            if isempty(obj.cov)
                obj = obj.get_cov;
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
            b = RT \ b;
        end

        % row trimming

        function obj = trim_rows(obj,varargin)
            p = inputParser;
            addParameter(p,'meth','std');
            addParameter(p,'trim_factor',10);
            parse(p,varargin{:})
            for i=1:obj.ntraj
                for j=1:obj.numeq
                    reduced_row_num =  floor(length(obj.bs{i}{j})/p.Results.trim_factor);        
                    if isequal(p.Results.meth,'std')
                        [~,inds_keep] = sort(std([obj.bs{i}{j} obj.Gs{i}{j}],[],2),'descend');
                    elseif isequal(p.Results.meth,'std_corner')
                                % min_rows = floor(length(b_0)/1)
                        % tol_cp = 0.5;
                        % [F,x] = histcounts(GG,floor(sqrt(size(GG,1))),'Normalization','pdf');
                        % dat = fliplr(cumsum(fliplr(F)));
                        % dat = dat/norm(dat);
                        % [mu,res_cp] = findchangepts(dat,'statistic','linear');
                        % inds_keep = GG>x(mu);
                        % if or(length(find(inds_keep))<min_rows,sqrt(res_cp)>tol_cp)
                        %     inds_keep = logical(ones(size(b_0)));
                        % else
                        %     figure(10)
                        %     findchangepts(dat,'statistic','linear')
                        % end
                    end
                    
                    inds_keep = ismember((1:length(obj.bs{i}{j}))',inds_keep(1:min(end,reduced_row_num)));

                    obj.Gs{i}{j} = obj.Gs{i}{j}(inds_keep,:);
                    obj.bs{i}{j} = obj.bs{i}{j}(inds_keep,:);
                    obj.L0{i}{j} = obj.L0{i}{j}(inds_keep,:);
                    S = obj.get_supp;
                    for k=1:length(S{j})
                        obj.L1{i}{j}{S{j}(k)} = obj.L1{i}{j}{S{j}(k)}(inds_keep,:);
                    end
                end
            end
            obj.cat_Gb;
        end

        % performance metrics
        function res = res(obj,meth)
            if ~exist('meth','var')
                if ~isempty(obj.weights)
                    % res = cat(1,obj.b{:}) - cat(1,obj.G{:})*obj.weights;
                    res = (cell2mat(obj.b) - blkdiag(obj.G{:})*obj.weights);
                else
                    res = [];
                end
            elseif isequal(meth,'sepcomp')
                res = cell(obj.numeq,1);
                if ~isempty(obj.weights)
                    s = obj.get_supp;
                    p = obj.get_params;
                    for j=1:obj.numeq
                        res{j} = (obj.b{j}-obj.G{j}(:,s{j})*p{j})/norm(obj.b{j});
                        disp(norm(res{j}))
                    end
                end
            end                
        end

        function res = res_ji(obj,j,i)
            if ~isempty(obj.weights)
                res = obj.bs{j}{i} - obj.Gs{j}{i}*obj.weights;
            else
                res = [];
            end
        end

        function T = termmags(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            w = IP.Results.w;
            if isequal(obj.catm ,'component')
                T = cellfun(@(w,s,G,b)vecnorm(G(:,s).*w')/norm(b),obj.get_params('w',w),obj.get_supp,obj.G,obj.b,'uni',0);
            elseif isequal(obj.catm ,'blkdiag')
                T = vecnorm(obj.G{1}.*w')/norm(obj.b{1});
                T = T(obj.weights~=0);
            end
        end

        function T = termproj(obj,varargin)
            IP = inputParser; 
            addParameter(IP, 'w', obj.weights);
            parse(IP, varargin{:});
            w = IP.Results.w;
            if isequal(obj.catm ,'component')
                T = cellfun(@(w,s,G,b)abs(b'*G(:,s))/norm(b)./vecnorm(G(:,s)),obj.get_params('w',w),obj.get_supp,obj.G,obj.b,'uni',0);
            elseif isequal(obj.catm ,'blkdiag')
                T = abs(obj.b{1}'*obj.G{1})/norm(obj.b{1})./vecnorm(obj.G{1});
                T = T(w~=0);
            end
        end

        % development features
        function x = get_adjoint_rhs(obj,ref_data)
            x = [];
            for i=1:obj.ntraj
                for j=1:obj.nstates
                    x = [x;obj.tf{i}{1}.get_testmat(0)*(ref_data(i).Uobs{j}-obj.dat(i).Uobs{j})];
                end
            end

        end

        function L = get_adjoint_L(obj)
            L0w = cell(obj.ntraj,1);
            L1w = cell(obj.ntraj,1);
            w = reshape_cell(obj.weights,arrayfun(@(L)length(L.terms),obj.lib)); 
            S = obj.get_supp;                
            for i=1:obj.ntraj
                L0_ = cell(size(obj.L1{i}));
                L1_ = cell(size(obj.L1{i}));
                for j=1:length(obj.L1{i})
                    dims = num2cell(size(obj.L0{i}{j}));
                    L1_{j} = sparse(dims{:});
                    for k=1:length(S{j})
                        L1_{j} = L1_{j} + w{j}(S{j}(k))*obj.L1{i}{j}{S{j}(k)};
                    end
                    B = zeros(size(L1_{j},1)*obj.nstates,size(L1_{j},2)/obj.nstates);
                    ind=0;
                    for k=1:obj.nstates
                        B(ind+(1:size(L1_{j},1)),:) = L1_{j}(:,(k-1)*size(L1_{j},2)/obj.nstates + (1:size(L1_{j},2)/obj.nstates) );
                        ind = ind+size(L1_{j},1);
                    end
                    L1_{j}=B;
                    B = zeros(size(L0_{j},1)*obj.nstates,size(L0_{j},2)/obj.nstates);
                    ind=0;
                    for k=1:obj.nstates
                        B(ind+(1:size(obj.L0{i}{j},1)),:) = obj.L0{i}{j}(:,(k-1)*size(obj.L0{i}{j},2)/obj.nstates + (1:size(obj.L0{i}{j},2)/obj.nstates) );
                        ind = ind+size(obj.L0{i}{j},1);
                    end
                    L0_{j}=B;
                end
                L0w{i} = cell2mat(L0_(:)');
                L1w{i} = cell2mat(L1_(:)');
            end
            L = cellfun(@(x,y) x - y, L0w,L1w,'uni',0);
            L = blkdiag(L{:});
        end

        function res0 = get_theoryres(obj,i,j)
            if isempty(obj.Gs)
                obj.get_Gb;
            end
            sig = obj.dat(i).estimate_sigma;
            Vn = norm(obj.tf{i}{j}.get_testmat(obj.lhsterms{j}.difftags),'fro');
            res0 = Vn*sig{j};
            b0 = norm(obj.bs{i}{j});
            res0 = res0/b0;
        end

        function [Fs,Js,Etags] = get_functional_forms(obj,varargin)

            p = inputParser;
            addParameter(p,'d',[]);
            addParameter(p,'w',[]);
            parse(p,varargin{:})

            d = p.Results.d;
            if ~isempty(d)
                d1 = digits;
                digits(d)
            end
            
            w = p.Results.w;
            if ~isempty(w)
                w = obj.reshape_w('w',w);
            else
                w = obj.reshape_w;
            end

            Xc = sym2cell(str2sym(strcat('x',num2str((1:obj.nstates)'))));
            Fs = cell(obj.numeq,1);
            Js = cell(obj.numeq,1);
            Etags = cell(obj.numeq,1);
            for i=1:obj.numeq
                X = obj.lib(i).terms;
                difftags = [];
                for j = 1:length(X)
                    try
                        difftags = [difftags;X{j}.linOp.difftags];
                    catch
                        difftags = [difftags;[0 0]];
                    end
                end
                
                uniq_tag = unique(difftags,'rows');
                fs = cell(size(uniq_tag,1),1);
                js = cell(size(uniq_tag,1),obj.nstates);
                for j=1:size(uniq_tag,1)
                    f = 0; J = num2cell(zeros(1,obj.nstates));
                    for k = find(ismember(difftags,uniq_tag(j,:),'rows'))'
                        f = f + sym(w{i}(k),'d')*X{k}.fHandle(Xc{:});
                        for m = 1:obj.nstates
                            J{m} = J{m} + sym(w{i}(k),'d')*X{k}.gradterms(m).fHandle(Xc{:});
                        end
                    end
                    fs{j} = simplify(f);
                    js(j,:) = J;
                    for m=1:obj.nstates
                        js{j,m} = simplify(js{j,m});
                    end
                end
                
                for j=1:size(uniq_tag,1)
                    fs{j} = matlabFunction(fs{j},'vars',Xc);
                    for m=1:obj.nstates
                        js{j,m} = matlabFunction(js{j,m},'vars',Xc);
                    end
                end
                Etags{i} = uniq_tag;
                Fs{i} = fs;
                Js{i} = js;


                if ~isempty(d)
                    digits(d1);
                end
                
            end

        end

        function [f,J,s,Fs] = get_functional_forms_vec(obj,varargin)

            p = inputParser;
            addParameter(p,'w',[]);
            parse(p,varargin{:})
            
            w = p.Results.w;
            if ~isempty(w)
                w = obj.reshape_w('w',w);
            else
                w = obj.reshape_w;
            end
        
            Fs = cell(obj.numeq,1);
            for i=1:obj.numeq
                X = obj.lib(i).terms;
                difftags = [];
                for j = 1:length(X)
                    try
                        difftags = [difftags;X{j}.linOp.difftags];
                    catch
                        difftags = [difftags;[0 0]];
                    end
                end
                
                uniq_tag = unique(difftags,'rows');
                fs = cell(size(uniq_tag,1),1);
                for j=1:size(uniq_tag,1)
                    inds = find(ismember(difftags,uniq_tag(j,:),'rows'))';
                    tags = zeros(length(inds),obj.nstates);
                    tags_J = zeros(length(inds),obj.nstates,obj.nstates);
                    W = zeros(length(inds),1);
                    W_J = zeros(length(inds),obj.nstates);
                    for k = 1:length(inds)
                        W(k) = w{i}(inds(k));
                        tags(k,:) = X{inds(k)}.ftag;
                        for ll=1:obj.nstates
                            tags_J(k,ll,:) = X{inds(k)}.gradterms(ll).ftag;
                        end
                        W_J(k,:) = W(k)*arrayfun(@(gt)gt.coeff,X{inds(k)}.gradterms);
                    end
                    fs{j} = {uniq_tag(j,:), tags, W, tags_J, W_J};
                end
                Fs{i} = fs;
            end 
            f = @(q) obj.eval_ff_vec(Fs,[1,0],q);
            s = @(q) obj.eval_ff_vec(Fs,[0,0],q);
            J = @(q) obj.eval_Jf_vec(Fs,[1,0],q);        
        end

        function X = eval_ff_vec(obj,Fs,difftag,q)
            X = q*0;
            for i=1:obj.numeq
                for j=1:length(Fs{i})
                    if isequal(Fs{i}{j}{1},difftag)
                        X(:,i) = obj.eval_ff(Fs{i}{j}{2},Fs{i}{j}{3},q);
                    end
                end
            end
        end

        function JX = eval_Jf_vec(obj,Fs,difftag,q)
            if isvector(q)
                q = q(:)';
            end
            JX = zeros(size(q,1),obj.nstates,obj.nstates);
            for i=1:obj.numeq
                for j=1:length(Fs{i})
                    num_terms = size(Fs{i}{j}{2},1);
                    if isequal(Fs{i}{j}{1},difftag)
                        for k=1:obj.nstates
                            JX(:,i,k) = obj.eval_ff(reshape(Fs{i}{j}{4}(:,k,:),num_terms,obj.nstates),Fs{i}{j}{5}(:,k),q);
                        end
                    end
                end
            end
        end

        function out = eval_ff(obj,tags,W,q)
            if size(q,2)~=size(tags,2)
                disp(['error: mismatch btw data dim and tems'])
            end
            out = zeros(size(q,1),1);
            for j=1:length(W)
                if W(j)~=0
                    out = out + W(j)*prod(q.^tags(j,:),2);
                end
            end
        end

        function diff_tags = get_diffs(obj)
            diff_tags = arrayfun(@(L)L.get_diffs(obj.ndims),obj.lib(:),'un',0);
            diff_tags = cell2mat(diff_tags);
            diff_tags = unique(diff_tags,'rows');
        end

        function A = get_stability_matrix(obj,varargin)

            p = inputParser;
            addParameter(p,'w',[]);
            parse(p,varargin{:})            
            w = p.Results.w;
            if isempty(w)
                w = obj.weights;
            end
            [~,~,~,Fs] = obj.get_functional_forms_vec('w',w);
            A = @(q,k)0;
            diff_tags = obj.get_diffs;
            
            %%% note this assumes left-hand side is first time derivative lhs_ord = 1;
            for k=1:size(diff_tags,1)
                diff_tag = diff_tags(k,:);
                A = @(q,k) A(q,k) + prod((1i*k).^diff_tag(1:end-1))*squeeze(obj.eval_Jf_vec(Fs,diff_tag,q));
            end

        end

        function [Acell,k2,ks,kmesh,kmesh_uns] = get_homog_state_matrices(obj,u0,varargin)
            p = inputParser;
            addParameter(p,'w',[]);
            parse(p,varargin{:})            
            w = p.Results.w;
            if isempty(w)
                w = obj.weights;
            end
            A = obj.get_stability_matrix('w',w);

            k_max = obj.dat.dims(1:end-1)/2;
            ks = arrayfun(@(k) (0:k),k_max,'un',0);
            kmesh = cell(obj.dat.ndims-1,1);
            [kmesh{:}] = ndgrid(ks{:});
            kmesh_uns = cell2mat(cellfun(@(g)g(:),kmesh(:)','un',0));
            k2 = vecnorm(kmesh_uns,2,2);
            kmesh = kmesh_uns*cellfun(@(g)2*pi/range(g),obj.dat.grid(1:end-1)') ;

            Acell = arrayfun(@(k)A(u0,kmesh(k,:)),1:size(kmesh,1),'un',0);
        end

        function [unstable_modes,eigs_all,k2,ks,kmesh] = get_homog_state_stability(obj,u0,varargin)
            p = inputParser;
            addParameter(p,'w',[]);
            parse(p,varargin{:})            
            w = p.Results.w;
            if isempty(w)
                w = obj.weights;
            end

            [Acell,k2,ks,kmesh,kmesh_uns] = obj.get_homog_state_matrices(u0,'w',w);
            
            eigs_all = zeros(size(kmesh,1),obj.nstates);
            unstable_modes = {};
            for k = 1:size(kmesh,1)
                [V,E] = eigs(Acell{k});
                E = diag(E);
                eigs_all(k,:) = E;
                if any(real(E)>0)
                    ind = find(real(E)>0);
                    unstable_modes = [unstable_modes;{V(:,ind),kmesh_uns(k,:),E(ind)}];
                end
            end
        end

        function plot_homog_state_stability(obj,u0,fignum,varargin)

            p = inputParser;
            addParameter(p,'w',[]);
            parse(p,varargin{:})            
            w = p.Results.w;
            if isempty(w)
                w = obj.weights;
            end
            [~,eigs_all,k2,ks] = obj.get_homog_state_stability(u0,'w',w);
            maxk2_unstable = max(k2(max(real(eigs_all),[],2)>0));
            most_unstable = max(real(eigs_all),[],2);
            [~,most_unstable] = max(most_unstable);
            most_unstable = k2(most_unstable);
            
            figure(fignum);clf
            subplot(3,1,1)
                histogram(real(eigs_all(:)))
                title('real eigs')
            subplot(3,1,2)
                histogram(imag(eigs_all(:)))
                title('imag eigs')
            subplot(3,1,3)
                if obj.dat.ndims==2
                    plot(ks{1},max(real(eigs_all),[],2),'o-',ks{1},ks{1}*0,'--')
                    if ~isempty(maxk2_unstable)
                        xlim([0 maxk2_unstable+2])
                    end
                elseif obj.dat.ndims==3
                    surf(ks{:},reshape(max(real(eigs_all),[],2),length(ks{1}),[]),'edgecolor','none')
                    view([0 90])
                    colorbar
                end
            
            if isempty(maxk2_unstable)    
                title('no unstable modes')
            else
                title(['max unstable mode:',num2str(maxk2_unstable), '; most unstable mode:',num2str(most_unstable)])
            end
        end

    end

end