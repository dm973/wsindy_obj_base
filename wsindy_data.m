%%% defines an object that contains the data used to learn a model
%%% all operations independent of model library, optimizer, etc are
%%% included here, and properties of the data (regular vs. irregular grid)
%%% are encoded here, so that downstream tasks are agnostic.

classdef wsindy_data < handle
    properties
        Uobs
        noise
        grid
        nstates
        ndims
        dims
        isUniform
        R0
        ks
        dv
        sigmas
        scales
    end

    methods
        function obj = wsindy_data(Uobs,grid,varargin)

            obj.Uobs = Uobs;
            if isequal(class(obj.Uobs),'double')
                obj.Uobs = mat2cell(obj.Uobs,size(obj.Uobs,1),ones(1,size(obj.Uobs,2)));
            end

            obj.grid = grid;
            if isequal(class(obj.grid),'double')
                if ~isempty(obj.grid)
                    x = size(obj.grid);
                    if x(1)<x(2)
                        obj.grid = obj.grid';
                    end
                    obj.grid = mat2cell(obj.grid,size(obj.grid,1),ones(1,size(obj.grid,2)));
                else
                    obj.grid = {(1:length(obj.Uobs{1}))'};
                end
            end

            default_isUniform = [];
            p = inputParser;
            addRequired(p,'Uobs');
            addRequired(p,'grid');
            addParameter(p,'isUniform',default_isUniform);
            addParameter(p,'sigmas',[]);
            parse(p,Uobs,grid,varargin{:})

            obj.nstates = length(obj.Uobs);
            obj = obj.get_dims;
            obj.ndims = length(obj.grid);
            obj.noise = [];
            if isempty(p.Results.isUniform)
                obj = obj.setUniformity;
            else
                obj.isUniform = p.Results.isUniform;
            end
            obj.dv = cellfun(@(g) mean(diff(g)), obj.grid);
            obj.sigmas = p.Results.sigmas;
            if isempty(obj.sigmas)
                obj.sigmas = cell(obj.nstates,1);
            end
            obj.set_scales(1);

        end

    end

    methods

        function obj = coarsen(obj,s,d)

            if ~exist('d','var')
                if isscalar(s)
                    d=1:obj.ndims;
                    s = s + d*0;
                else
                    if length(s)<obj.ndims
                        s = [repmat(s(1),1,obj.ndims-1),s(end)];
                    end
                    d=1:min(length(s),obj.ndims);
                end
            else
                d = d(1:min(length(s),obj.ndims));
            end

            for i=1:length(d)
                if s(i)>0
                    obj = coarsen_dim(obj,s(i),d(i));
                else
                    s_temp = max(floor(obj.dims(i)/(-s(i))),1);
                    obj = coarsen_dim(obj,s_temp,d(i));
                end
            end
            obj = obj.get_dims;
        end

        function obj = trimend(obj,s,d)
            if length(s)==1
                if ~exist('d','var')
                    d=1;
                end
                obj = trimend_dim(obj,s,d);
            else
                if ~exist('d','var')
                    d=1:length(s);
                end
                for i=1:length(s)
                    if length(d)~=length(s)
                        obj = trimend_dim(obj,s(i),i);
                    else
                        obj = trimend_dim(obj,s(i),d(i));
                    end
                end
            end
            obj = obj.get_dims;
        end

        function obj = trimstart(obj,s,d)
            if length(s)==1
                if ~exist('d','var')
                    d=1;
                end
                obj = trimstart_dim(obj,s,d);
            else
                if ~exist('d','var')
                    d=1:length(s);
                end
                for i=1:length(s)
                    if length(d)~=length(s)
                        obj = trimstart_dim(obj,s(i),i);
                    else
                        obj = trimstart_dim(obj,s(i),d(i));
                    end
                end
            end
            obj = obj.get_dims;
        end

        function obj = addnoise(obj,nr,varargin)

            if nr~=0
                default_seed = 'shuffle';
                default_dist = 0;
                default_alg = 0;
                default_uniform = 0;
    
                p = inputParser;
                addRequired(p,'obj');
                addRequired(p,'nr');
                addParameter(p,'seed',default_seed);
                addParameter(p,'dist',default_dist);
                addParameter(p,'alg',default_alg);
                addParameter(p,'uniform',default_uniform);
                parse(p,obj,nr,varargin{:})
    
                rng(p.Results.seed);
                obj.noise = cell(obj.nstates,2);

                if p.Results.uniform == 1
                    U_cat = cell2mat(obj.Uobs(:)');
                    [U_cat,N,~,sigma] = gen_noise(U_cat,p.Results.nr,p.Results.dist,p.Results.alg);
                    for n=1:obj.nstates
                        obj.noise(n,:) = {N(:,n),sigma};
                        obj.Uobs{n} = U_cat(:,n);
                    end
                else
                    for n=1:obj.nstates
                        [U,N,~,sigma] = gen_noise(obj.Uobs{n},p.Results.nr,p.Results.dist,p.Results.alg);
                        obj.noise(n,:) = {N,sigma};
                        obj.Uobs{n} = U;
                    end
                end
            end

        end

        function m = get_mode(obj,stateind,modenum,dim)
            x = permute(obj.Uobs{stateind},[dim 1:dim-1 dim+1:max(2,obj.ndims)]);
            x = x - mean(x);
            x = cat(1,x,flipud(x));
            xf = reshape(abs(fft(x)),2*obj.dims(dim),[]);
            xf = mean(xf(1:floor(end/2),:),2);

            Ufft = cumsum(flipud(xf));
            NN = length(Ufft);
            xx = 0:(NN-1);
            Ufft = Ufft/max(abs(Ufft))*NN;
            m = diff(Ufft([1 end]))/diff(xx([1 end]));
            Rtheta = [[1 m];[-m 1]]/sqrt(1+m^2);
            U_tilt = Rtheta(2,:)*[xx(:)';Ufft(:)'];

            [a,b] = findpeaks(-fliplr(U_tilt),'MinPeakDistance',8);
            if  modenum==0
                [~,ii]=max(a);
                loc = b(ii);
            else
                loc = b(min(modenum,end));
            end
            m = obj.dims(dim)/((loc-1)/2);

%             figure(50)
%             findpeaks(-fliplr(U_tilt),'MinPeakDistance',8)
%             hold on
%             plot([loc,loc],[min(-U_tilt) max(-U_tilt)])
%             hold off
%             drawnow

%             [~,I] = sort(xf,'descend');
%             m = obj.dims(dim)/((I(modenum)-1)/2);


        end

        function obj = get_R0(obj)
            if isempty(obj.R0)
                for j=1:obj.nstates
                    if isempty(obj.sigmas{j})
                        obj.sigmas{j} = obj.estimate_sigma{j};
                    end
                end
                M = prod(obj.dims);
                obj.R0 = spdiags(kron(reshape([obj.sigmas{:}].^2,[],1),ones(M,1)),0,M*obj.nstates,M*obj.nstates);
            end
        end

        function plotDyn(obj,varargin)

            default_ax = gca;
            p = inputParser;
            addParameter(p,'ax',default_ax);
            parse(p,varargin{:})
            
            if obj.ndims==1
                for n=1:obj.nstates
                    plot(p.Results.ax,obj.grid{1},obj.Uobs{n},'markersize',5,'linewidth',2)
                    hold on
                    if n==1
                        title(num2str(obj.dims))
                    end
                end
                hold off
                legend
            end

        end

        function plotPhase(obj,varargin)

            default_ax = gca;
            default_noisy = 1;
            default_coord = reshape(1:obj.nstates,2,[]);
            p = inputParser;
            addParameter(p,'ax',default_ax);
            addParameter(p,'coord',default_coord);
            addParameter(p,'noisy',default_noisy);
            parse(p,varargin{:})
        
            for j=1:size(p.Results.coord,2)
                if or(p.Results.noisy,isempty(obj.noise))
                    plot(p.Results.ax,obj.Uobs{p.Results.coord(1,j)},obj.Uobs{p.Results.coord(2,j)},'markersize',5,'linewidth',2)
                else
                    plot(p.Results.ax,obj.Uobs{p.Results.coord(1,j)}-obj.noise{p.Results.coord(1,j),1},obj.Uobs{p.Results.coord(2,j)}-obj.noise{p.Results.coord(2,j),1},'markersize',5,'linewidth',2)
                end
                hold on
            end
            hold off
        end

        function sig = estimate_sigma(obj,varargin)

            default_ord = min(floor((min(obj.dims)-1)/2),6);
            p = inputParser;
            addParameter(p,'ord',default_ord);
            addParameter(p,'alg','AWGN');
            addParameter(p,'set',false);
            parse(p,varargin{:})

            ord = p.Results.ord;
            alg = p.Results.alg;
 
            C = fdcoeffF(ord,0,-ord:ord);
            filter = C(:,end);
            filter = filter/norm(filter,2);
            if obj.ndims > 1
                [~,ind] = max(size(obj.Uobs{1}));
                sig = cellfun(@(U)rms(reshape(convn(permute(U,[ind 1:ind-1 ind+1:obj.ndims]),filter(:),'valid'),[],1)),obj.Uobs,'uni',0);
            else
                if isequal(alg,'AWGN')
                    sig = cellfun(@(U)rms(conv(U,filter(:),'valid')),obj.Uobs,'uni',0);
                elseif isequal(alg,'logn')
                    sig = cellfun(@(U)rms(conv(log(U),filter(:),'valid')),obj.Uobs);
                    sig = sqrt((exp(sig.^2)-1).*exp(2*cellfun(@(U)mean(log(U)),obj.Uobs)+(sig.^2)));
                    sig = num2cell(sig);
                end

            end

            if p.Results.set
                obj.sigmas = sig;
            end

        end

        function x0 = get_x0(obj,rhs,varargin)

            default_ind = 1;
            default_Tf = obj.dims(1)/4;
            default_wd_test = 10:60;
            default_polydeg = 2;
            default_tol = 10^-6;
            default_Hfun = [];
            default_meth = 'ind';

            p = inputParser;
            addRequired(p,'rhs');
            addParameter(p,'ind',default_ind);
            addParameter(p,'Hfun',default_Hfun);
            addParameter(p,'Tf',default_Tf);
            addParameter(p,'wd_test',default_wd_test);
            addParameter(p,'polydeg',default_polydeg);
            addParameter(p,'tol',default_tol);
            addParameter(p,'meth',default_meth);
            parse(p,rhs,varargin{:})

            Tf = p.Results.Tf;
            ind = p.Results.ind;
            wd_test = p.Results.wd_test;
            polydeg = p.Results.polydeg;
            tol = p.Results.tol;
            Hfun = p.Results.Hfun;
            meth = p.Results.meth;
            if and(isempty(Hfun),isequal(meth,'H=const'))
                disp('No Hamiltonian found, switching to default x0 getter')
                meth = default_meth;
            end

            if isequal(meth,'grid_search')
                err = zeros(length(wd_test),1);
                x0_reduced = zeros(length(wd_test),obj.nstates);
                options_ode_sim = odeset('RelTol',tol,'AbsTol',tol*ones(1,obj.nstates));
                for ll=1:length(wd_test)
                    c = arrayfun(@(i)polyfit(0:wd_test(ll)-1,obj.Uobs{i}(1:wd_test(ll)),polydeg),1:obj.nstates,'uni',0);
                    x0_reduced(ll,:) = cellfun(@(x)x(end),c);
                    [~,xH0_true]=ode15s(@(t,x)rhs(x),obj.grid{1}(1:Tf),x0_reduced(ll,:),options_ode_sim);
                    err(ll) = mean(arrayfun(@(i) norm(obj.Uobs{i}(1:size(xH0_true,1)) - xH0_true(:,i),inf),1:obj.nstates));
%                     figure(ll)
%                     for i=1:obj.nstates
%                         plot(1:Tf,obj.Uobs{i}(1:Tf),'k-o',1:Tf,xH0_true(:,i),'b-s',1:wd_test(ll),polyval(c{i},0:wd_test(ll)-1))
%                         hold on
%                     end
%                     hold off
%                     drawnow
                end
                [~,a]=min(err);
                x0 = x0_reduced(a,:);
            elseif isequal(meth,'trap')
                x0 = mean(cell2mat(obj.Uobs)-cumtrapz(obj.grid{1}, arrayfunvec(cell2mat(obj.Uobs)',rhs,1)'));
            elseif isequal(meth,'H=const')
                dat = Hfun_true(obj.Uobs{:});
                ii = find(abs(dat-mean(dat))<0.1*range(dat),1);
                x0 = arrayfun(@(i)obj.Uobs{i}(ii),1:obj.nstates);
            elseif isequal(meth,'ind')
                x0 = arrayfun(@(i)obj.Uobs{i}(ind),1:obj.nstates);
                if ~isempty(obj.noise)
                    for i=1:length(x0)
                        x0(i) = x0(i) - obj.noise{i,1}(ind);
                    end
                end
            end

        end

        function [Ufft,xx] = get_fft(obj,stateind,dim)
            xobs = permute(obj.Uobs{stateind},[dim 1:dim-1 dim+1:max(2,obj.ndims)]);
%                     xobs = xobs - mean(xobs);
%                     xobs = cat(1,xobs,flipud(xobs));
            t = obj.grid{dim}(:);
%                     t = [t;t(end)+mean(diff(t))+t-t(1)];
        
            T = length(t);
            xx = (0:floor(length(t)/2)-1)'*(2*pi)/range(t);
            NN = length(xx);

            Ufft = mean(reshape(abs(fft(xobs)),T,[]),2) /sqrt(NN);
            Ufft = Ufft(1:NN,:);
        end

        function obj = getcorners(obj,toggle_plot)
            if ~exist('toggle_plot','var')
                toggle_plot=0;
            end
        
            obj.ks = zeros(obj.nstates,obj.ndims);
            for j=1:obj.nstates
                for i=1:obj.ndims 
                    [Ufft,xx] = obj.get_fft(j,i);
                    if toggle_plot
                        tplot = toggle_plot+(j-1)*obj.ndims+i-1;
                    else
                        tplot = 0;
                    end
                    Ufft(1)=0;
                    obj.ks(j,i) = getcorner(Ufft,xx,tplot);
                end
            end        
        end

        function obj = plotFFTtf(obj,tf,startfig)
            for j=1:obj.nstates
                for i=1:obj.ndims
                    figure(startfig+1+(j-1)*obj.ndims+i-1)
                    [Ufft,xx] = obj.get_fft(j,i);
                    NN = length(xx);
                    k = obj.ks(j,i);
                    Ufft = Ufft/max(Ufft);
                    tf_dat = abs(tf{min(end,j)}.Cfsfft{i}(1,1:NN));
                    tf_dat = tf_dat/max(tf_dat);
                    semilogy(0:NN-1,Ufft,'k-',0:NN-1,tf_dat,'r-',k,Ufft(k),'gd','markersize',10)
                    legend({'Ufft','tf','k'})
                end
            end     
        end

        function obj = undo_scales(obj)
            if ~isempty(obj.scales)
                obj.Uobs = cellfun(@(U,s)U*s,obj.Uobs,num2cell(obj.scales(1:obj.nstates)),'un',0);
                obj.grid = cellfun(@(x,s)x*s,obj.grid,num2cell(obj.scales(obj.nstates+1:end)),'un',0);
                obj.dv = obj.dv(:)'.*obj.scales(obj.nstates+1:end);
                obj.scales = ones(1,obj.nstates+obj.ndims);
            end
        end

        function obj = set_scales(obj,scl,varargin)
            p = inputParser;
            addParameter(p,'nrm',[]);
            addParameter(p,'val',2);
            addParameter(p,'lib',[]);
            addParameter(p,'tf',[]);
            parse(p,varargin{:})
            
            nrm = p.Results.nrm;
            val = p.Results.val;
            if isscalar(val)
                val = repmat({val},1,obj.nstates+obj.ndims);
            end
            lib = p.Results.lib;
            tf = p.Results.tf;

            if isequal(length(scl),obj.nstates+obj.ndims)
                scales = scl;
            else
                scales = [ones(1,obj.nstates) ones(1,obj.ndims)];
                try            
                    pd = max(cell2mat(cellfun(@(ttf) cellfun(@(p) p(end),ttf.param),tf(:),'un',0)),[],1);
                    md = max(cell2mat(cellfun(@(ttf) ttf.rads(:)',tf(:),'un',0)),[],1);
                    dx = obj.dv(:)';
                    betad = max(cell2mat(arrayfun(@(L) max(cell2mat(L.tags(:)),[],1), lib(:), 'un',0)),[],1);
                    betad = betad(1:obj.nstates);
                    ad = max(cell2mat(arrayfun(@(L) max(cell2mat(cellfun(@(tt)tt.linOp.difftags,L.terms(:),'un',0)),[],1), lib(:), 'un',0)),[],1);
                    scales_u = arrayfun(@(b,i) (norm(obj.Uobs{i}(:).^b)/norm(obj.Uobs{i}(:)))^(1/b),betad,1:obj.nstates);
                    scales_x = arrayfun(@(p,m,d,a) 1./((nchoosek(p,a/2)*factorial(a))^(1/a)/(m*d)),pd,md,dx,ad);
                    scales = [scales_u,scales_x];
                catch
                    if ~isempty(nrm)
                        scales(1:obj.nstates) = cellfun(@(U,v)norm(U(:),nrm)/v,obj.Uobs,val(1:obj.nstates));
                        scales(obj.nstates+1:end) = cellfun(@(x,v)norm(x(:),nrm)/v,obj.grid,val(obj.nstates+1:end));
                    else        
                        if or(isequal(scl,'Ubar'),~isequal(length(scl),obj.nstates+obj.ndims))
                            if isempty(scl)
                                scales(1:obj.nstates) = cellfun(@(U)mean(abs(U(:))),obj.Uobs);
                                scales(obj.nstates+1:end) = cellfun(@(x)mean(abs(x(:))),obj.grid);
                            elseif isequal(scl,'Ubar')
                                scales(1:obj.nstates) = cellfun(@(U)mean(abs(U(:))),obj.Uobs);
                            end
                        end
                    end
                end
            end
            obj.Uobs = arrayfun(@(i)obj.Uobs{i}/scales(i),1:obj.nstates,'un',0);
            obj.grid = arrayfun(@(i)obj.grid{i}/scales(obj.nstates+i),1:obj.ndims,'un',0);
            obj.dv = obj.dv(:)'./scales(obj.nstates+1:end);
            obj.scales = scales;
        end

    end

    methods (Hidden = true)

        function obj = get_dims(obj)
            obj.dims = cellfun(@(x)length(x),obj.grid(:));
        end

        function obj = coarsen_dim(obj,s,d)
            if ~exist('d','var')
                d=1;
            end
            s = max(s,1);
            inds = cell(1,obj.ndims);
            inds{d} = 1:s:size(obj.Uobs{1},d); 
            inds([1:d-1 d+1:obj.ndims]) = repmat({':'},obj.ndims-1,1);
            obj.Uobs = cellfun(@(U) U(inds{:}), obj.Uobs, 'uni',0);
            for i=1:obj.ndims
                if isvector(obj.grid{i})
                    if i==d
                        obj.grid{d} = obj.grid{d}(1:s:end);
                    end                    
                else 
                    obj.grid{i} = obj.grid{i}(inds{:});
                end
            end
        end

        function obj = trimend_dim(obj,s,d)
            if ~exist('d','var')
                d=1;
            end
            inds = cell(1,obj.ndims);
            inds{d} = 1:min(s,size(obj.Uobs{1},d)); 
            inds([1:d-1 d+1:obj.ndims]) = repmat({':'},obj.ndims-1,1);
            obj.Uobs = cellfun(@(U) U(inds{:}), obj.Uobs, 'uni',0);
            for i=1:obj.ndims
                if isvector(obj.grid{i})
                    if i==d
                        obj.grid{d} = obj.grid{d}(1:min(s,end));
                    end                    
                else 
                    obj.grid{i} = obj.grid{i}(inds{:});
                end
            end
        end

        function obj = trimstart_dim(obj,s,d)
            if ~exist('d','var')
                d=1;
            end
            inds = cell(1,obj.ndims);
            inds{d} = min(s,size(obj.Uobs{1},d)):size(obj.Uobs{1},d); 
            inds([1:d-1 d+1:obj.ndims]) = repmat({':'},obj.ndims-1,1);
            obj.Uobs = cellfun(@(U) U(inds{:}), obj.Uobs, 'uni',0);
            for i=1:obj.ndims
                if isvector(obj.grid{i})
                    if i==d
                        obj.grid{d} = obj.grid{d}(inds{d});
                    end                    
                else 
                    obj.grid{i} = obj.grid{i}(inds{:});
                end
            end
        end
 
        function obj = setUniformity(obj)
            obj.isUniform = cellfun(@(t) max(abs(diff(t(:))-mean(diff(t(:)))))<10^-12,obj.grid);
        end


    end

end
