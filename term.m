classdef term < absterm
    properties
        coeff
        tol
    end

    methods
        function obj = term(varargin)
            default_linOp = [];
            default_ftag = [];
            default_fHandle = [];
            default_gradterms = [];
            default_nstates = 0;
            default_coeff = 1;
            default_gradon = 1;
            default_tol = eps;

            p = inputParser;
            addParameter(p,'ftag',default_ftag);
            addParameter(p,'linOp',default_linOp);
            addParameter(p,'fHandle',default_fHandle);
            addParameter(p,'gradterms',default_gradterms);
            addParameter(p,'nstates',default_nstates);
            addParameter(p,'coeff',default_coeff);
            addParameter(p,'gradon',default_gradon);
            addParameter(p,'tol',default_tol);
            parse(p,varargin{:})

            obj.ftag = p.Results.ftag;
            obj.linOp = p.Results.linOp;
            obj.fHandle = p.Results.fHandle;
            obj.gradterms = p.Results.gradterms;
            obj.nstates = p.Results.nstates;
            obj.coeff = p.Results.coeff;
            obj.gradon = p.Results.gradon;
            obj.tol = p.Results.tol;


            if ~isempty(obj.ftag)
                obj = obj.set_tag;
            elseif ~isempty(obj.fHandle)
                obj = obj.set_fHandle;
            end

            obj.set_linOp(obj.linOp);

            % if ~isempty(obj.linOp)
            %     if isequal(class(obj.linOp),'double')
            %         if ~isequal(obj.linOp,0)
            %             obj.linOp = diffOp(obj.linOp,'nstates',obj.nstates);
            %         else
            %             obj.linOp = [];
            %         end
            %     end
            % end

            % if and(obj.gradon,isempty(obj.gradterms))
            %     obj = obj.get_grads;
            % end

        end

    end

    methods

        function Y = evalterm(obj,dat) % only fHandle
            if isequal(class(dat),'cell')
                Y = obj.fHandle(dat{:});
            elseif isequal(class(dat),'wsindy_data')
                s = dat.scales;
                if and(obj.get_scale(s) == 1,~isempty(s))
                    s = num2cell(s);
                    U = dat.Uobs;
                    U = cellfun(@(u,s)u*s,U,s(1:dat.nstates),'un',0);
                    Y = obj.fHandle(U{:});
                else
                    Y = obj.fHandle(dat.Uobs{:});            
                end
            elseif isequal(class(dat),'double')
                Xcell = obj.dat2cell(dat);
                Y = obj.fHandle(Xcell{:});
            end
        end

        function Y = evaltermLinOp(obj,dat)
            if isempty(obj.linOp)
                Y = obj.evalterm(dat);
            else         
                if isequal(class(obj.linOp),'diffOp')
                    if isempty(obj.linOp.Dmats)
                        obj.linOp.get_Dmats(dat);
                    end
                end
                if dat.ndims==1
                    Y = obj.linOp.evalterm(obj.evalterm(dat));
                else
                    Y = obj.linOp.evalterm({obj.evalterm(dat)});
                end
            end
        end

        function Y = evalgrads(obj,dat)
            if isempty(obj.gradterms)
                obj = obj.get_grads;
                obj.gradon = 1;
            end
            Y = arrayfun(@(g)g.evalterm(dat),obj.gradterms,'uni',0);
        end

        function Y = diffmat(obj,dat)
            Y = cellfun(@(g) spdiags(g(:),0,length(g(:)),length(g(:))),evalgrads(obj,dat),'uni',0);
        end

        function tt = get_taylor(obj,pt,ord,toggleH)
            tt = [];
            args = num2cell(pt);
            if and(isequal(class(obj.ftag),'double'),size(obj.ftag,1) == 1)
                bool = imag(obj.ftag) == 0;
            else 
                bool = zeros(obj.nstates,1);
            end
            for i=1:obj.nstates
                if bool(i)
                    t1 = arrayfun(@(j)nchoosek(obj.ftag(i),j).*pt(i).^j,(0:obj.ftag(i))')*obj.fHandle(args{:})/pt(i)^obj.ftag(i);
                    t1 = t1(max(end-ord+1,1):end);
                else
                    q = sym('q');
                    arg_temp = args;
                    arg_temp{i} = q+1;
                    t1 = sym2poly(taylor(obj.fHandle(arg_temp{:}),q,0,'order',ord))';
                end
                if length(t1)<ord
                    t1 = [zeros(ord-length(t1),1);t1];
                end
                if toggleH==1
                    tt = [tt;flipud(t1(1:end-1))];
                else
                    tt = [tt;flipud(t1)];
                end
            end
        end

        function obj = get_grads(obj)
            obj.gradterms = repelem(term('gradon',0,'nstates',obj.nstates),1,obj.nstates);
            if and(~isempty(obj.ftag),isequal(class(obj.ftag),'double'))
                if and(all(real(obj.ftag)==imag(obj.ftag)),any(obj.ftag))
                    gradtags = real(obj.ftag)*(1-1i);
                    for j=1:obj.nstates
                        obj.gradterms(j) = term('ftag',gradtags,'coeff',-real(obj.ftag(j))*obj.coeff,'gradon',0,'nstates',obj.nstates,'linOp',obj.linOp);
                    end
                elseif and(all(real(obj.ftag)==-imag(obj.ftag)),any(obj.ftag))
                    gradtags = real(obj.ftag)*(1+1i);
                    for j=1:obj.nstates
                        obj.gradterms(j) = term('ftag',gradtags,'coeff',real(obj.ftag(j))*obj.coeff,'gradon',0,'nstates',obj.nstates,'linOp',obj.linOp);
                    end
                else
                    kp = arrayfun(@(k)obj.unitgradtag(k),obj.ftag,'uni',0);
                    gradtags = repmat(obj.ftag,obj.nstates,1);
                    gradtags(logical(eye(obj.nstates))) = cellfun(@(kp)kp(1),kp);
                    for j=1:obj.nstates
                        obj.gradterms(j) = term('ftag',gradtags(j,:),'coeff',kp{j}(2)*obj.coeff,'gradon',0,'nstates',obj.nstates,'linOp',obj.linOp);
                    end
                end
            elseif ~isempty(obj.fHandle)
                X = str2sym(strcat('x',num2str((1:obj.nstates)')));
                Xc = sym2cell(X);
                for j=1:obj.nstates
                    g = diff(obj.fHandle(Xc{:}),X(j),1);
                    % if all(arrayfun(@(i)diff(g,X(i),1)==0,1:obj.nstates)) %g==0
                    if isSymType(g,'number') %g==constant
                        xstr=reshape(strcat(',x',num2str((1:obj.nstates)'))',[],1)';
                        obj.gradterms(j) = term('fHandle',eval(['@(',xstr(2:end),')x1*0+double(g)']),'gradon',0,'nstates',obj.nstates,'linOp',obj.linOp);
                    else
                        obj.gradterms(j) = term('fHandle',matlabFunction(g,'vars',Xc),'gradon',0,'nstates',obj.nstates,'linOp',obj.linOp);
                    end                    
                end
            end

        end

        function s = get_str(obj)
            s = '';
            if ~isempty(obj.fHandle)
                s = [s,obj.get_str_0];
            end
            if isequal(class(obj.linOp),'diffOp')
                s = ['D^[',strrep(num2str(obj.linOp.difftags),'  ',','),']',s];
            end
            % if obj.coeff~=1
            %     s = [num2str(obj.coeff),'*',s];
            % end
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
            x = strfind(s,')')+1;
            s = s(x(1):end);
        end

        function rhs = get_rhs(obj)
            if isequal(class(obj.linOp),'diffOp')
                ns = find([obj.gradterms.coeff]);
                if and(isscalar(obj.linOp.difftags),isscalar(ns))
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

        function obj = set_coeff(obj,coeff)
            obj.coeff = coeff;
            obj.get_fHandle;
            if obj.gradon == 1
                obj = obj.get_grads;
            end
        end

        function obj = set_linOp(obj,linOp)

            if ~isempty(linOp)
                if isequal(class(linOp),'double')
                    if ~isequal(linOp,0)
                        obj.linOp = diffOp(linOp,'nstates',obj.nstates);
                    else
                        obj.linOp = [];
                    end
                end
            end

            if obj.gradon == 1
                obj = obj.get_grads;
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
            
            if isequal(obj.ftag,zeros(1,obj.nstates))
                obj.fHandle = eval(['@(',expr,')',num2str(obj.coeff),'+0*',vars{1}]);
            else
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
                    inds_cos = all([real(obj.ftag)==imag(obj.ftag);obj.ftag~=0],1);
                    inds_sin = all([real(obj.ftag)==-imag(obj.ftag);obj.ftag~=0],1);
                    inds_other = find(all([~inds_cos;~inds_sin;obj.ftag~=0],1));
    
                    fUnits = arrayfun(@(k)unitfstr(obj,k),obj.ftag(inds_other),'uni',0);
                    fstr = [];
                    for i=1:length(fUnits)
                        fstr = [fstr,strrep(fUnits{i},'x',vars{inds_other(i)})];
                        if i<length(fUnits)
                            fstr = [fstr,'.*'];
                        end
                    end
                    if ~isempty(fstr)
                        str_fun= ['@(',expr,')',num2str(obj.coeff),'*',fstr];
                    else
                        str_fun= ['@(',expr,')',num2str(obj.coeff)];
                    end
                    if length(find(inds_cos))>length(find(inds_sin))
                        str = arrayfun(@(i)[num2str(real(obj.ftag(i))),'*',vars{i}],find(inds_cos),'uni',0);
                        for j=1:length(find(inds_cos))
                            if ~(str{j}(1)=='-')
                                str{j}=['+',str{j}];
                            end
                        end
                        str = strcat(str{:});
                        str_fun = [str_fun,'*','cos(',str,')'];
                    elseif length(find(inds_cos))<length(find(inds_sin))
                        str = arrayfun(@(i)[num2str(real(obj.ftag(i))),'*',vars{i}],find(inds_sin),'uni',0);
                        for j=1:length(find(inds_sin))
                            if ~(str{j}(1)=='-')
                                str{j}=['+',str{j}];
                            end
                        end
                        str = strcat(str{:});
                        str_fun = [str_fun,'*','sin(',str,')'];
                    % else
                    %     str_fun = [str_fun,'+0*',vars{1}];
                    end
                    obj.fHandle = eval(str_fun);
    
                    % 
                    % 
                    % fUnits = arrayfun(@(k)unitfstr(obj,k),obj.ftag,'uni',0);
                    % fstr = [];
                    % for i=1:obj.nstates
                    %     fstr = [fstr,strrep(fUnits{i},'x',vars{i})];
                    %     if i<obj.nstates
                    %         fstr = [fstr,'.*'];
                    %     end
                    % end
                    % obj.fHandle = eval(['@(',expr,')',num2str(obj.coeff),'*',fstr]);
                end
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

        function Xcell = dat2cell(obj,X)
            n = find(size(X)==obj.nstates);
            if ~isempty(n)
                if n(1)==1
                    Xcell = mat2cell(X,ones(obj.nstates,1),size(X,2));
                elseif n(1)==2
                    Xcell = mat2cell(X,size(X,1),ones(1,obj.nstates))';
                end
            else
                Xcell = {X};
            end
        end

        function m = get_scale(obj,scales)
            m = 1;
            if ~isempty(scales)
                if isequal(class(obj.ftag),'double')
                    if all(isreal(obj.ftag))
                        for k=1:size(obj.ftag,1)
                            m=m*prod(scales(1:obj.nstates).^-obj.ftag(k,:));
                        end
                    end
                end
                if ~isempty(obj.linOp)
                    if isequal(class(obj.linOp),'diffOp')
                        m=m*prod(scales(obj.nstates+1:end).^obj.linOp.difftags);
                    end
                end
            end
        end

        function lap = get_lap(obj)
            if isempty(obj.gradterms)
                obj.get_grads;
            end
            for n=1:obj.nstates
                obj.gradterms(n).get_grads;
            end
            lap = arrayfun(@(n)obj.gradterms(n).gradterms(n),1:obj.nstates);
            % for n=2:obj.nstates
            %     lap = addterm(lap,obj.gradterms(n).gradterms(n));
            % end
        end
        function hess = get_hess(obj)
            if isempty(obj.gradterms)
                obj.get_grads;
            end
            for n=1:obj.nstates
                obj.gradterms(n).get_grads;
            end
            hess = arrayfun(@(n)obj.gradterms(n).gradterms(:),(1:obj.nstates)','uni',0);
            hess = vertcat(hess{:});
            % hess = cell2mat(hess);
            % for n=2:obj.nstates
            %     lap = addterm(lap,obj.gradterms(n).gradterms(n));
            % end
        end
    end

end
