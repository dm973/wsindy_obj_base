classdef gridterm < absterm
    properties
        ndims
    end

    methods
        function obj = gridterm(varargin)
            p = inputParser;
            addParameter(p,'ftag',[]);
            addParameter(p,'fHandle',[]);
            addParameter(p,'linOp',[]);
            addParameter(p,'gradon', 1);
            addParameter(p,'nstates', 1);
            addParameter(p,'ndims', 2);
            addParameter(p,'dat', []);
            parse(p,varargin{:})

            dat = p.Results.dat;
            if ~isempty(dat)
                obj.ndims = dat.ndims;
                obj.nstates = dat.nstates;
            else
                obj.ndims = p.Results.ndims;
                obj.nstates = p.Results.nstates;
            end

            obj.ftag = zeros(1,obj.nstates);
            obj.fHandle = p.Results.fHandle;
            obj.linOp = p.Results.linOp;
            obj.gradon = p.Results.gradon;

            obj.gradterms = repelem(term('gradon',0,'nstates',obj.nstates),1,obj.nstates);
            for j=1:obj.nstates
                xstr = reshape(strcat(',x',num2str((1:obj.nstates)'))',[],1)';
                obj.gradterms(j) = term('fHandle',eval(['@(',xstr(2:end),')x1*0']),'gradon',0,'nstates',obj.nstates);
            end

        end
    end

    methods
        function Y = evalterm(obj,dat)
            g = dat.grid;
            g = cellfun(@(g) g(:), g, 'un', 0);
            for i=2:obj.ndims
                g{i} = permute(g{i},fliplr(1:i));
            end
            % g = arrayfun(@(i) permute(g{i},fliplr(1:max(2,i))), [1:obj.ndims], 'un',0);
            Y = obj.fHandle(g{:});
        end

        function Y = evaltermLinOp(obj,dat)
            Y = obj.evalterm(dat);
        end

        function s = get_str(obj)
            s = functions(obj.fHandle).function;
        end

        function m = get_scale(obj,scales)
            m = 1;
            if and(~isempty(scales),~isempty(obj.ftag))
                for i=0:obj.ndims-1
                    m = m*obj.ftag(end-i)^scales(end-i);
                end
            end
        end

        function Y = diffmat(obj,dat)
            Y = repelem({sparse(prod(dat.dims),prod(dat.dims))},1,obj.nstates);
        end

    end

end
