function [Ai_cell,bi_cell] = stable_diffusion_constraints(WS,tol)

    Ai_cell = cell(WS.numeq,1);
    bi_cell = cell(WS.numeq,1);
    for n = 1:WS.numeq
        ts = WS.lib(n).terms;
        T = length(ts);
        Ai = zeros(T);
        for t=1:T
            LO = ts{t}.linOp;
            if isequal(class(LO),'diffOp')
                if sum(LO.difftags)==2
                    Ai(t,t) = -1;
                end
            end
        end
        Ai = Ai( logical(sum(Ai ~= 0,2)), :);
        bi = -tol*ones(size(Ai,1),1);
        Ai_cell{n} = Ai;
        bi_cell{n} = bi;
    end

    if isequal(WS.catm,'blkdiag')
        Ai_cell = {blkdiag(Ai_cell{:})};
        bi_cell = {cell2mat(bi_cell)};
    end

end
