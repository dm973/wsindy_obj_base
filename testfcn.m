%%% defines a test function object that evaluates terms weakly at a 
%%% given dataset. That is, given a term object tm and wsindy_data 
%%% object, Uobj testfcn.test(Uobj,tm) discretizes the functional
%%% T : B(U)->R^k given by T(tm(Uobj)). When T is linear and k=1, 
%%% this is a duality pairing <T,tm(Uobj)>. Details on the quadrature rule
%%% are defined in Uobj, integration by parts occurs depending on the structure
%%% of tm.

classdef testfcn < handle

    properties
        stateind
        nstates
        ndims
        dims
        phifuns
        rads
        Cfs
        Cfsfft
        subinds
        dv
        meth
        param
        avg
        mtmin
        mtmax
        pre_conv
        Kmax
    end

    methods
        
        function obj = testfcn(dat,varargin)

            if isequal(class(dat),'double')
                dat = wsindy_data(dat,[]);
            end

            obj.nstates = dat.nstates;
            obj.ndims = dat.ndims;
            obj.dims = dat.dims;

            default_phifuns = {@(t) exp(9.*((t.^2-1).^(-1)+1))};
            default_rads = [];
            default_subinds = [];
            default_meth = 'FFT';
            default_param = 2;
            default_avg = 0;
            default_stateind = 1;
            default_mtmin = repelem(3,obj.ndims,1);
            default_mtmax = cellfun(@(g) floor((length(g)-1)/2), dat.grid);
            default_pre_conv = [];
            default_Kmax = inf;

            p = inputParser;
            addRequired(p,'dat');
            addParameter(p,'phifuns',default_phifuns);
            addParameter(p,'rads',default_rads);
            addParameter(p,'meth',default_meth);
            addParameter(p,'param',default_param);
            addParameter(p,'avg',default_avg);
            addParameter(p,'subinds',default_subinds);
            addParameter(p,'stateind',default_stateind);
            addParameter(p,'mtmin',default_mtmin);
            addParameter(p,'mtmax',default_mtmax);
            addParameter(p,'pre_conv',default_pre_conv);
            addParameter(p,'Kmax',default_Kmax);
            parse(p,dat,varargin{:})

            obj.dv = cellfun(@(g)mean(diff(g)),dat.grid);
            obj.phifuns = p.Results.phifuns;
            obj.rads = p.Results.rads;
            obj.Cfs = cell(obj.ndims,1);
            obj.Cfsfft = cell(obj.ndims,1);
            obj.subinds = p.Results.subinds;
            obj.meth = p.Results.meth;
            obj.param = p.Results.param;
            if length(obj.param)~=obj.ndims
                obj.param = repelem(obj.param(1),1,obj.ndims);
            end
            obj.avg = p.Results.avg;
            obj.stateind = p.Results.stateind;
            obj.mtmin = p.Results.mtmin;
            obj.mtmax = p.Results.mtmax;
            obj.pre_conv = p.Results.pre_conv;
            obj.Kmax= p.Results.Kmax;

            if isequal(class(obj.phifuns),'function_handle')
                obj.phifuns = {obj.phifuns};
            elseif isequal(obj.phifuns,'delta')
                obj.phifuns = {@(x)x==0};
            elseif isequal(obj.phifuns,'pp')
                obj.phifuns = {'pp'};
            elseif isequal(obj.phifuns,'Gauss')
                obj.phifuns = {'Gauss'};
            elseif isempty(obj.phifuns)
                obj.phifuns = default_phifuns;
            end

            if length(obj.phifuns)~=obj.ndims
                obj.phifuns = repelem(obj.phifuns,obj.ndims,1);
            end

            for j=1:obj.ndims
                if isequal(obj.phifuns{j},'delta')
                    obj.phifuns{j} = @(x)x==0;
                end
            end

            obj.get_rads(dat);
            obj.get_subinds(dat);

        end

    end

    methods

        function vec = test(obj,dat,term)
            if isequal(class(term),'double')
                Cfsinds = term;
                obj = obj.add_difftags(Cfsinds);
                vec = obj.op(dat,Cfsinds);
            elseif or(isequal(class(term),'diffOp'),any(cellfun(@(x) isequal(x,'diffOp'),superclasses(term))))
                Cfsinds = term.difftags;
                obj = obj.add_difftags(Cfsinds);
                vec = obj.op(dat.Uobs{term.stateind},Cfsinds);
            elseif isequal(class(term.linOp),'diffOp')
                L = term.linOp;
                Cfsinds = L.difftags;
                obj = obj.add_difftags(Cfsinds);
                vec = obj.op(term.evalterm(dat),Cfsinds);
            elseif isequal(class(term),'addterm')
                vec = obj.test(dat,term.t1)+obj.test(dat,term.t2);
            elseif isempty(term.linOp)
                Cfsinds = zeros(obj.ndims,1);
                obj = obj.add_difftags(Cfsinds);
                vec = obj.op(term.evalterm(dat),Cfsinds);
            end
            vec = vec(:);
        end

        function V = get_testmat(obj,difftag)
            if isequal(class(difftag),'diffOp')
                difftag = difftag.difftags;
            elseif length(difftag)~=obj.ndims
                difftag = zeros(obj.ndims,1);
            end
            obj.add_difftags(difftag);
            V = 1;
            for i=1:obj.ndims
                Vtemp = spdiags(repmat(fliplr(obj.Cfs{i}(difftag(i)+1,:)),obj.dims(i)-size(obj.Cfs{i},2)+1,1),...
                        0:size(obj.Cfs{i},2)-1, obj.dims(i)-size(obj.Cfs{i},2)+1, obj.dims(i));
                Vtemp = Vtemp(unique(min(obj.subinds{i},end)),:);
                V = kron(Vtemp,V);
            end
        end

        function obj = get_rads(obj,dat)
            obj.rads = zeros(obj.ndims,1);
            for i=1:obj.ndims
                if any([isempty(obj.phifuns{i}) isempty(obj.meth) isempty(obj.param)])
                    mt = obj.mtmin(i);
                else
                    if isequal(obj.meth,'direct')
                        mt = obj.param(i);
                    elseif isequal(obj.meth,'FFT')
                        dat.getcorners;
                        if isequal(class(obj.phifuns{i}),'function_handle')
                            mt = get_tf_support(obj.phifuns{i},dat.dims(i),obj.param(i),dat.ks(obj.stateind,i));
                        else
                            [obj.phifuns{i},mt] = obj.get_phi_handle(dat.ks(obj.stateind,i),dat.dims(i),obj.phifuns{i},obj.param{i});
                        end
                    elseif isequal(obj.meth,'timefrac')
                        mt = floor(length(dat.grid{i})*obj.param(i));
                    elseif isequal(obj.meth,'mode1frac')
                        mt = dat.get_mode(obj.stateind,obj.param{i}(1),i);
                        mt = ceil(mt*obj.param{i}(2));
                    elseif isequal(obj.meth,'mtmin')
                        mt = obj.param(i)*obj.mtmin;
                    end
                    mt = min(max(obj.mtmin(i),mt),obj.mtmax(i));
                end
                obj.rads(i) = mt;
            end
            if obj.avg
                obj.rads = obj.rads*0 + floor(mean(obj.rads));
            end
        end

        function obj = get_subinds(obj,dat)
            if isequal(class(obj.subinds),'double')
                if isempty(obj.subinds)
                    obj.subinds = arrayfun(@(i)1:max(1,floor(obj.rads(i)/4)):dat.dims(i)-2*obj.rads(i),1:dat.ndims,'uni',0);
                elseif isscalar(obj.subinds)
                    if obj.subinds>=0
                        obj.subinds = arrayfun(@(i)1:obj.subinds:dat.dims(i)-2*obj.rads(i),1:dat.ndims,'uni',0);
                    else
                        obj.subinds = arrayfun(@(i)1:max(1,floor(obj.rads(i)/-obj.subinds)):dat.dims(i)-2*obj.rads(i),1:dat.ndims,'uni',0);
                    end
                elseif length(obj.subinds)==dat.ndims
                    obj.subinds = arrayfun(@(i)1:obj.subinds(i):dat.dims(i)-2*obj.rads(i),1:dat.ndims,'uni',0);
                elseif isvector(obj.subinds)
                    obj.subinds = repmat({obj.subinds},1,dat.ndims);
                end
            end
            obj.subinds = cellfun(@(s)s(1:min(end,obj.Kmax)),obj.subinds,'uni',0);
        end
    
        function obj = get_Cfsfft(obj)
            for k=1:obj.ndims
                [m,n] = size(obj.Cfs{k});
                if size(obj.Cfsfft{k},1)~=m
                    obj.Cfsfft{k} = fft([obj.Cfs{k}(:,end) zeros(m,obj.dims(k)-n) obj.Cfs{k}(:,1:end-1)],[],2);
                end
            end
        end

        function obj = add_difftags(obj,Cfsinds)
            for k=1:obj.ndims
                if size(obj.Cfs{k},1)-1<Cfsinds(k)
                    Cfs_new = obj.phi_weights(k,size(obj.Cfs{k},1):Cfsinds(k));
                    if ~isempty(obj.pre_conv)
                        obj.Cfs{k} = [obj.Cfs{k};arrayfunvec(Cfs_new,@(y)conv(obj.pre_conv(:)',y,'valid'),2)];
                    else
                        obj.Cfs{k} = [obj.Cfs{k};Cfs_new];
                    end
                end
%                 obj.Cfs{k} = obj.Cfs{k}/norm(obj.Cfs{k}(1,:),2);
            end
            obj = obj.get_Cfsfft;
        end

        function X = op(obj,X,Cfsinds)
            if isvector(X)
                X = conv(X(:),obj.Cfs{1}(Cfsinds+1,:),'valid');
                X = X(unique(min(obj.subinds{1},obj.dims(1)-2*obj.rads(1))));
            else
                cols = arrayfun(@(k) obj.Cfsfft{k}(Cfsinds(k)+1,:), 1:obj.ndims, 'uni',0);
                for k=1:obj.ndims
                    col_fft = cols{k}(:);                    
                    if obj.ndims == 1
                        shift = [1 2];
                        shift_back = shift;
                    else
                        shift = circshift(1:obj.ndims,1-k);
                        shift_back=circshift(1:obj.ndims,k-1);
                    end                
                    X = ifft(col_fft.*fft(permute(X,shift)));
                    inds = cell(obj.ndims,1);
                    inds{1} = unique(min(obj.subinds{k},obj.dims(k)-2*obj.rads(k))); 
                    inds(2:obj.ndims) = repmat({':'},obj.ndims-1,1);
                    X = X(inds{:});
                    X = permute(X,shift_back);
                end
                X = real(X);
            end
        end

        function Cfs = phi_weights(obj,k,diffs)
            xf = linspace(-1,1,2*obj.rads(k)+1);
            if ~isequal(functions(obj.phifuns{k}).function,'@(x)x==0')
                x = xf(2:end-1);
                Cfs = zeros(length(diffs),2*obj.rads(k)+1);
                syms y;
                for j=1:length(diffs)
                    Df = matlabFunction(diff(obj.phifuns{k}(y),diffs(j)));
                    Cfs(j,2:end-1) = fillmissing(Df(x),'constant',Df(eps));
                    inds = find(isinf(abs(Cfs(j,:))));
                    for i=1:length(inds)
                        Cfs(j,inds(i)) = Df(xf(inds(i))-sign(xf(inds(i)))*eps);
                    end
                end
            else
                if obj.rads(k)>=1
                    fdcoeffs = fdcoeffF(max(diffs),0,xf);
                    Cfs = fdcoeffs(:,diffs+1)'.*(-1).^diffs';
                elseif obj.rads(k)>0
                    Cfs = fliplr(eye(max(diffs)+1));
                    Cfs = Cfs(diffs+1,:);
                else
                    Cfs = ones(max(diffs)+1,1);
                    Cfs = Cfs(diffs+1);
                end
            end
            if obj.rads(k)>=1
                Cfs = Cfs.*(obj.rads(k)*obj.dv(k)).^(-diffs(:))*obj.dv(k);
            end
            % Cfs = Cfs/norm(Cfs(1,:),1);
        end

        function [phifun,m,p] = get_phi_handle(obj,k,N,phi_class,param)
            if length(param)~=3
                tau = 10^-10;
                tauhat = 2;
                maxd = 1;
            else
                tau = param(1);
                tauhat = param(2);
                maxd = param(3);
            end
            if isequal(phi_class,'pp')
                l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
                mstar1 = sqrt(3)/pi*N/2/k*tauhat;
                mstar2 = 1/pi*tauhat*(N/2)/k*sqrt(log(exp(1)^3/tau^8));
                m = floor(min(fzero(@(m)l(m,k,N), [mstar1 mstar2]),(N-1)/2));
                p = max(maxd+1,ceil(log(tau)/log(1-(1-1/m)^2)));
                phifun = @(x) (1-x.^2).^p;
            elseif isequal(phi_class,'Gauss')
                m = floor(min(1+N*tauhat/2/pi/k*sqrt(-2*log(tau)),(N-1)/2));
                p = 2*pi*k/tauhat/N;
                phifun = @(x) exp(-(m*p*x).^2/2); % p = dx/sig
            end 
        end
    end

end