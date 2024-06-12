classdef termx < absterm
    properties
        ndims
    end

    methods
        function obj = termx(varargin)
            p = inputParser;
            addParameter(p,'ftag',[]);
            addParameter(p,'fHandle',@(varargin)varargin{1}*0);
            addParameter(p,'ndims',1);
            parse(p,varargin{:})

            obj.ndims = p.Results.ndims;
            obj.ftag = p.Results.ftag;
            obj.fHandle = p.Results.fHandle;

            if isempty(obj.ftag)
                obj.ftag = zeros(1,obj.ndims);
            end
        end

    end

    methods

        function Y = evalterm(obj,dat) % only fHandle
            if isequal(class(dat),'cell')
                Y = obj.fHandle(dat{:});
            elseif isequal(class(dat),'wsindy_data')
                g = dat.grid;
                g{1} = g{1}(:);
                for i=2:length(g)
                    e = num2cell(ones(1,i-1));
                    g{i} = reshape(g{i},e{:},length(g{i}));
                end
                Y = obj.fHandle(g{:});            
            elseif isequal(class(dat),'double')
                Y = obj.fHandle(dat);
            end
        end

        function s = get_str(obj)
            s = '';
            if ~isempty(obj.fHandle)
                s = [s,obj.get_str_0];
            end
            if isequal(class(obj.linOp),'diffOp')
                s = ['(d/dt)^',num2str(obj.linOp.difftags),s];
            end
        end

        function s = get_str_0(obj)
            s = functions(obj.fHandle).function;
            % s = strrep(s,'x1','x');
            % s = strrep(s,'x2','z');
            % s = strrep(s,'x3','dx');
            % s = strrep(s,'x4','dz');
            for i=1:obj.nstates
                s = strrep(s,['(2.2204e-16+abs(x',num2str(i),'))'],['x',num2str(i)]);
            end
            % s = strrep(s,'(2.2204e-16+abs(x1))','x1');
            % s = strrep(s,'(2.2204e-16+abs(z))','z');
            % s = strrep(s,').^-1','.^-1');
            % s = strrep(s,').^-2','.^-2');
            % s = strrep(s,'((x))','(x)');
            % s = strrep(s,'((z))','(z)');
            % s = strrep(s,'.*','');
            % s = strrep(s,'.^','^');
            s = strrep(s,')1*',')');
            s = s(strfind(s,')')+1:end);
        end

        function rhs = get_rhs(obj)
            if isequal(class(obj.linOp),'diffOp')
                ns = find([obj.gradterms.coeff]);
                if and(length(obj.linOp.difftags)==1,length(ns)==1)
                    % if obj.ftag(ns)==1
                    %     rhs = @(varargin) varargin{ns+obj.nstates};
                    % else
                        e1 = repmat({0},1,ns-1); 
                        e2 = repmat({0},1,obj.nstates-ns);
                        f = W2S(@(x)obj.fHandle(e1{:},x,e2{:}),obj.linOp.difftags);
                        rhs = @(varargin) f([varargin{ns:obj.nstates:end}]);
                    % end
                % elseif and(length(obj.linOp.difftags)==1,length(find(obj.ftag))==2)
                %     inds = find(obj.ftag);
                %     e1 = repmat({0},1,ns-1); e2 = repmat({0},1,obj.nstates-ns);
                %     f = W2S(@(x)obj.fHandle(e1{:},x,e2{:}),obj.linOp.difftags);
                %     rhs = @(varargin) f([varargin{ns:obj.nstates:end}]);                    
                end
            else; rhs = @(varargin) obj.fHandle(varargin{1:obj.nstates});
            end
        end


    end

    methods (Hidden = true)
        function obj = set_fHandle(obj)
            if obj.nstates == 0
                a = strfind(functions(obj.fHandle).function,')');
                args = functions(obj.fHandle).function(1:a(1));
                obj.nstates = 1+length(strfind(args,','));
            end
        end

        function obj = set_tag(obj)
            if and(isequal(class(obj.ftag),'double'),size(obj.ftag,1)==1)
                obj.nstates = length(obj.ftag);
                obj.get_fHandle;
            elseif and(isequal(class(obj.ftag),'double'),size(obj.ftag,1)>1)
                foo = obj.ftag;
                obj = prodterms(term('ftag',obj.ftag(1,:)),term('ftag',obj.ftag(2,:)),obj.gradon);
                for j=3:size(foo,1)
                    obj = prodterms(obj,term('ftag',foo(j,:)),obj.gradon);
                end
            end
        end

        function obj = get_fHandle(obj)
            vars = arrayfun(@(i)['x',num2str(i)],1:obj.nstates,'uni',0);
            expr = [];
            for i=1:obj.nstates
                expr = [expr,vars{i}];
                if i<obj.nstates
                    expr = [expr,','];
                end
            end

            if and(all(real(obj.ftag)==imag(obj.ftag)),any(obj.ftag))
                str = arrayfun(@(i)[num2str(real(obj.ftag(i))),'*',vars{i}],1:obj.nstates,'uni',0);
                for j=1:obj.nstates
                    if ~(str{j}(1)=='-')
                        str{j}=['+',str{j}];
                    end
                end
                str = strcat(str{:});
                obj.fHandle = eval(['@(',expr,')',num2str(obj.coeff),'*','cos(',str,')']);
            elseif and(all(real(obj.ftag)==-imag(obj.ftag)),any(obj.ftag))
                str = arrayfun(@(i)[num2str(real(obj.ftag(i))),'*',vars{i}],1:obj.nstates,'uni',0);
                for j=1:obj.nstates
                    if ~(str{j}(1)=='-')
                        str{j}=['+',str{j}];
                    end
                end
                str = strcat(str{:});
                obj.fHandle = eval(['@(',expr,')',num2str(obj.coeff),'*','sin(',str,')']);
            else
                fUnits = arrayfun(@(k)unitfstr(obj,k),obj.ftag,'uni',0);
                fstr = [];
                for i=1:obj.nstates
                    fstr = [fstr,strrep(fUnits{i},'x',vars{i})];
                    if i<obj.nstates
                        fstr = [fstr,'.*'];
                    end
                end
                obj.fHandle = eval(['@(',expr,')',num2str(obj.coeff),'*',fstr]);
            end
        end

        function f = unitfstr(obj,k)
            if imag(k)==0
                if and(k>=0,floor(k)==k)
                    f = ['x.^',num2str(k)];
                elseif and(k>=0,floor(k)~=k)
                    f = ['abs(x).^',num2str(k)];
                elseif k<0
                    if k>=-0.1
                        f = ['log(',num2str(obj.tol),'+abs(x))'];
                    else
                        f = ['(',num2str(obj.tol),'+abs(x)).^',num2str(k)];
                    end
                end
            else
                if imag(k)>0
                    f = ['cos(',num2str(imag(k)),'*x)'];
                else
                    f = ['sin(',num2str(-imag(k)),'*x)'];
                end
            end
        end

        function kp = unitgradtag(obj,k)
            if imag(k)==0
                if or(k>0,k<-0.1)
                    kp = k-1;
                    c = k;
                elseif and(k>=-0.1,k<0)
                    kp = -1;
                    c = 1;
                else
                    kp = 0;
                    c = 0;
                end
            else
                kp = -k;
                c = -imag(k);
            end
            kp=[kp;c];
        end

        function C = prodcell(obj,C)

            d = size(C{1});
            if any(d==1)
                ind = find(d==1,1);
            else
                ind = length(d)+1;
            end
            C = prod(cat(ind,C{:}),ind);
        end

        function m = get_scale(obj,scales)
            m = 1;
            if isequal(class(obj.ftag),'double')
                if all(isreal(obj.ftag))
                    m=m*prod(scales(1:obj.nstates).^obj.ftag);
                end
            end
            if ~isempty(obj.linOp)
                if isequal(class(obj.linOp),'diffOp')
                    m=m/prod(scales(obj.nstates+1:end).^obj.linOp.difftags);
                end
            end
        end


    end

end
