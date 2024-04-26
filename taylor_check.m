syms q p
format long
pt = mean(cell2mat(arrayfun(@(j)arrayfun(@(i)mean(Uobj(j).Uobs{i}),1:nstates),(1:ntraj)','uni',0)),1);
rr = max(cell2mat(arrayfun(@(j)arrayfun(@(i)range(Uobj(j).Uobs{i}),1:nstates),(1:ntraj)','uni',0)),1);

t1 = sym2poly(taylor(Hfun(q,0)-C,q,pt(1),'order',10))';
t2 = sym2poly(taylor(Hfun_true(q,0)-C_true,q,pt(1),'order',length(t1)))';
t1'
t2'
h = matlabFunction(Hfun_true(q,0)-C_true);

s = 2;
figure(ntraj+1)
subplot(2,1,1)
r = rr(1)/2*s;
x = linspace(pt(1)-r,pt(1)+r,200);
plot(x,h(x),'b')
hold on
plot(x,Hfun(x,0)-C,'g--')
hold off
title(norm(h(x)-Hfun(x,0)+C)/norm(h(x)))

subplot(2,1,2)
t1 = sym2poly(taylor(Hfun(0,p)-C,p,pt(2),'order',10))';
t2 = sym2poly(taylor(Hfun_true(0,p)-C_true,p,pt(2),'order',length(t1)))';
t2=[zeros(length(t1)-length(t2),1);t2];
h = matlabFunction(Hfun_true(0,p)-C_true);
r = rr(2)/2*s;
x = linspace(pt(2)-r,pt(2)+r,200);
plot(x,h(x),'b')
hold on
plot(x,Hfun(0,x)-C,'g--')
hold off
title(norm(h(x)-Hfun(0,x)+C)/norm(h(x)))
