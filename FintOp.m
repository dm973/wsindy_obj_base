classdef FintOp < linearOp % Fredholm integral operator
    properties
        f_in % inner function of state
        K % function of grid 
        Kmat
        ndims
        dims
        svdtol
    end

    methods
        function obj = FintOp(varargin)
            p = inputParser;
            addParameter(p,'K',@(varargin) varargin{1}*0);
            addParameter(p,'f_in',[]);
            addParameter(p,'svdtol',1e-6);
            parse(p,varargin{:})

            obj.f_in = p.Results.f_in;
            if isempty(obj.f_in)
                obj.nstates = 1;
                obj.f_in = term(ftag,zeros(1,obj.nstates));
            else
                obj.nstates = obj.f_in.nstates;
            end
            obj.K = p.Results.K;
            obj.Kmat = [];
            obj.svdtol = p.Results.svdtol;

            obj.dims = [];
            obj.ndims = [];
        end

    end

    methods

        function Y = evaltermLinOp(obj,dat)
            Y = obj.evalterm(dat);
        end

        function Y = evalterm(obj,dat)
            if ~isequal(dat.dims,obj.dims)
                obj.get_Kmat(dat);
            end
            Y = obj.f_in.evalterm(dat);
            Y = obj.convNDfftsketch(Y);
        end

        function obj = get_Kmat(obj,dat)
            obj.dims = dat.dims;
            obj.ndims = length(obj.dims);
            dx = cellfun(@(g) mean(diff(g)), dat.grid(1:end-1));
            xK = dat.grid(1:obj.ndims-1);
            for k=1:obj.ndims-1
                xK{k} = dat.grid{k}(:) - min(dat.grid{k});
                inds = repmat({1},1,obj.ndims); inds{k} = 2*obj.dims(k)-1;
                xK{k} = reshape([-flipud(xK{k}(2:end));xK{k}],inds{:});
            end
            [xK{:}] = ndgrid(xK{:});
            obj.Kmat = dx*obj.K(xK{:});    
            if obj.ndims ==3        
                try 
                    [Uk,Sk,Vk] = svdsketch(obj.Kmat,obj.svdtol);
                catch
                    [Uk,Sk,Vk] = svdsketch(obj.Kmat);
                end
                trunc = find(vecnorm(diag(Sk) - triu(repmat(diag(Sk),1,size(Sk,1))))/norm(diag(Sk))<obj.svdtol,1);
                Sk = Sk(1:trunc,1:trunc);
                Uk = Uk(:,1:trunc);
                Vk = Vk(:,1:trunc);
                obj.Kmat = struct('Uk',Uk,'Sk',Sk,'Vk',Vk);
            end
        end

        function KconvY = convNDfftsketch(obj,Y)
            if obj.ndims==2
                KconvY = conv2(Y,obj.Kmat,'same');
            elseif obj.ndims==3
                KconvY = zeros(size(Y));
                for i=1:trunc
                    KconvY = KconvY + obj.Kmat.Sk(i,i)*permute(convn(...
                                                        permute(convn(...
                            Y,obj.Kmat.Uk(:,i),'same'),[2 1 3]),obj.Kmat.Vk(:,i),'same'),[2 1 3]);
                end
            end
        end

        function s = get_str(obj)
            f_in_tag = obj.f_in.get_str;
            s = [functions(obj.K).function,'*',f_in_tag];
        end

        function m = get_scale(obj,scales)
            m = 1;
            if ~isempty(scales)
                m = obj.f_in.get_scale(scales);
                m = m*prod(scales(obj.nstates+1:end));
            end
        end

        function K = diffmat(obj,dat,varargin)

            p = inputParser;
            addParameter(p,'toggle_noisefree',false);
            parse(p,varargin{:})


            D0 = obj.f_in.diffmat(dat);

            k = obj.Kmat;
            m = dat.dims(1); n = dat.dims(2);
            p = length(k);
            
            % 1D convolution matrix for one column
            c = [k; zeros(m-p,1)];
            r = [k(1); zeros(m-1,1)];
            T = toeplitz(c, r);
            
            % Keep 'same' rows
            center = ceil(p/2);
            T = T(center:center+m-1, :);
            
            % Full 2D operator
            X = kron(speye(n), sparse(T));
            K = cellfun(@(A) X*A,D0,'un',0);
        end


        % function [rhs,t,Kweights] = get_rhs(obj,tau,N)
        %     t = (0:tau:(N-1)*tau)';
        %     Kweights = obj.K(t)*tau;
        %     Kweights([1,end]) = Kweights([1,end])/2;
        %     rhs = @(y) dot(obj.f_in.evalterm(y),Kweights);            
        % end

        % function plot_kernel(obj,dat)
        %     plot(dat.grid{1},obj.Kgrid)
        % end

    end

end