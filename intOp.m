classdef intOp < linearOp
    properties
        f_in % inner function of state
        K % function of time
        ndims
    end

    methods
        function obj = intOp(ndims,varargin)
            p = inputParser;
            addRequired(p,'ndims');
            addParameter(p,'K',[]);
            addParameter(p,'nstates',1);
            addParameter(p,'f_in',[]);
            parse(p,ndims,varargin{:})

            obj.nstates = p.Results.nstates;
            obj.ndims = p.Results.ndims;

            obj.f_in = p.Results.f_in;
            if isempty(obj.f_in)
                obj.f_in = term(ftag,zeros(1,obj.nstates));
            end

            obj.K = p.Results.K;
            if isempty(obj.K)
                obj.K = @(varargin) varargin{1}*0;
            end
        end

    end

    methods

        function Y = evalterm(obj,dat)
            
            Y = obj.f_in.evalterm(dat);
            
            % ndg  = cell(obj.ndims,1);
            % [ndg{:}] = ndgrid(dat.grid{:});
            % Kgrid = obj.K(ndg{:});

            Kgrid = obj.K(dat.grid{1});
            K_tri = zeros(dat.dims, dat.dims);
            for i=1:dat.dims
                K_tri(i,1:i) = Kgrid(end:-1:end-i+1);
            end
            Y = K_tri*Y(:)*dat.dv;
            
            % Y = convn(Y,Kgrid,'same');

        end

        function s = get_str(obj)
            s = ['intOp'];
            % 
            % if obj.stateind==1
            %     s = ['(d/dt)^',num2str(obj.difftags),'(x)'];
            % elseif obj.stateind==2
            %     s = ['(d/dt)^',num2str(obj.difftags),'(z)'];
            % end

        end

    end

end
