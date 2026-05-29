classdef custom_term < term
    properties
        evaltermfun
        getstrfun
    end

    methods
        function obj = custom_term(evaltermfun,getstrfun)
            obj = obj@term('gradon',0);
            obj.evaltermfun = evaltermfun;
            obj.getstrfun = getstrfun;
        end

    end

    methods

        function Y = evalterm(obj,dat,tf)
            Y = obj.evaltermfun(dat,tf);
        end

        function str = get_str(obj)
            str=obj.getstrfun;
        end
    end

end
