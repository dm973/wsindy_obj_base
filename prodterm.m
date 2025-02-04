classdef prodterm < term
    properties
        t1
        t2
    end

    methods
        function obj = prodterm(t1,t2,varargin)
            obj = obj@term('gradon',0);
            default_gradon = 0;
            default_gradterms = {};
            default_linOp = [];
            p = inputParser;
            addRequired(p,'t1');
            addRequired(p,'t2');
            addParameter(p,'gradon',default_gradon);
            addParameter(p,'gradterms',default_gradterms);
            addParameter(p,'linOp',default_linOp);
            parse(p,t1,t2,varargin{:})

            obj.nstates = t1.nstates;
            obj.gradterms = p.Results.gradterms;
            obj.gradon = p.Results.gradon;
            obj.t1 = t1;
            obj.t2 = t2;
            obj.linOp = p.Results.linOp;
            obj.ftag = {'prod',t1.ftag,t2.ftag};

            if ~isempty(obj.linOp)
                if isequal(class(obj.linOp),'double')
                    obj.linOp = diffOp(obj.linOp,'nstates',obj.nstates);
                end
            end

            if and(obj.gradon,isempty(obj.gradterms))
                obj.get_grads;
            end
        end
    end

    methods
        function Y = evalterm(obj,dat)
            Y = obj.t1.evaltermLinOp(dat).*obj.t2.evaltermLinOp(dat);
        end

        function obj = get_grads(obj)
            obj.gradterms = repmat(addterm.empty,1,obj.nstates);
            if obj.t1.gradon==0
                obj.t1.get_grads;
            end
            if obj.t2.gradon==0
                obj.t2.get_grads;
            end
            for j=1:obj.nstates
                obj.gradterms(j) = addterm(prodterm(obj.t1,obj.t2.gradterms(j),'gradon',0),prodterm(obj.t1.gradterms(j),obj.t2,'gradon',0));
            end
        end

        function s = get_str(obj)
            s = ['(',obj.t1.get_str,')*(',obj.t2.get_str,')'];
            if isequal(class(obj.linOp),'diffOp')
                    s = ['(d/dt)^',num2str(obj.linOp.difftags),s];
            end
        end
        
        function m = get_scale(obj,scales)
            m = obj.t1.get_scale(scales)*obj.t2.get_scale(scales);
        end

        function rhs = get_rhs(obj)
            if isempty(obj.linOp)
                rhs1 = obj.t1.get_rhs;
                rhs2 = obj.t2.get_rhs;
                rhs = @(varargin)rhs1(varargin{:}).*rhs2(varargin{:});
            else
                %%% note, assumes t1,t2 have no linOp
                rhs1 = obj.t1.get_rhs;
                rhs2 = obj.t2.get_rhs;
                t1_new = term('ftag',obj.t1.ftag,'fHandle',obj.t1.fHandle,'linOp',obj.linOp,'nstates',obj.nstates);
                rhs1_new = t1_new.get_rhs;
                t2_new = term('ftag',obj.t2.ftag,'fHandle',obj.t2.fHandle,'linOp',obj.linOp,'nstates',obj.nstates);
                rhs2_new = t2_new.get_rhs;
                rhs = @(varargin)rhs1(varargin{:}).*rhs2_new(varargin{:})+rhs2(varargin{:}).*rhs1_new(varargin{:});
            end
        end

    end

end
