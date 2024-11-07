errs = vecnorm(w_its - w_true,2,1)/max(norm(w_true),eps);
err_windy=errs(end);

%%% compute conf intervals
w_jac = w_plot(w_plot~=0);
xflip = [1:length(w_jac) length(w_jac):-1:1];
c = 0.05; % <(100)c % chance of not containing true val
stdW = sqrt(diag(CovW));
conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
confbounds = [w_jac-conf_int;flipud(w_jac+conf_int)];

%%% print results
disp(['err WENDy:',num2str(err_windy)])

%%% wendy iterates
subplot(3,2,1)
semilogy(1:length(errs),errs,'bo-')
legend({'err(wsindyJac)'},'location','best')
ylabel('||w-w^*||_2/||w^*||_2')
title(['err(1)=',num2str(errs(1)),'; ',...
    'err(end)=',num2str(errs(end))])
ylims = [10^(log10(min(errs)))*0.9 10^(log10(max(errs)))*1.2];
ylim(ylims)
set(gca,'Ytick',exp(linspace(log(ylims(1)),log(ylims(2)),8)))

subplot(3,2,2)
dat = vecnorm(diff(w_its,[],2),2,1)./vecnorm(w_its(:,1:end-1),2,1);
semilogy(dat,'o-')
title('||w^n-w^{n+1}||/||w^n||')
set(gca,'Ytick',10.^(floor(min(log10(dat))):1:ceil(max(log10(dat)))))

%%% confidence intervals
subplot(3,2,3)
h=fill(xflip,confbounds,'red','FaceColor',[1 0.8 0.8],'edgecolor','none'); hold on;
try
    h1=plot(1:length(w_jac),w_jac,'ro',1:length(w_jac),w_true(w_plot~=0),'bx');
    legend(h1(1:2),{'w','w^*'})
catch
    h1=plot(1:length(w_jac),w_jac,'ro');
    legend('w')
end
hold off

%%% p-values
subplot(3,2,4)
pvals = arrayfunvec(res,@(v)outn(@swtest,v,2),1);
plot(pvals,'o-')
title(['p-val=',num2str(pvals(end))])
xlabel('iter')
legend('p-val')


%%% GLS residual: WENDy final to true
subplot(3,2,5)
plot(res(:,end)); 

subplot(3,2,6)
plot(res_0(:,1))
