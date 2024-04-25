function t3 = prodterms(t1,t2,gradon)
    try
        s1 = func2str(t1.fHandle);
        [vars,expr] = get_vars_str(s1);
        s2 = func2str(t2.fHandle);
        ftemp = eval([expr,'(',strrep(s1,expr,''),').*(',fcn_rep_vars(s2,vars),')']);
        testp = repmat({1},1,length(vars));
        ftemp(testp{:});
    catch
        ftemp = @(varargin)t1.fHandle(varargin{:}).*t2.fHandle(varargin{:});
    end
    t3 = term('fHandle',ftemp,'gradon',0);
    t3.nstates = t1.nstates;
    if isequal(class(t1.ftag),class(t2.ftag))
        t3.ftag = [t1.ftag;t2.ftag];
    else
        t3.ftag = {t1.ftag,t2.ftag};
    end
    if gradon
        t3.gradon = 1;
        t3.gradterms = repmat(term('gradon',0),1,t3.nstates);
        for j=1:t3.nstates
            t3.gradterms(j) = addterms(prodterms(t1,t2.gradterms(j),0),prodterms(t1.gradterms(j),t2,0),0);
        end
    end
end
