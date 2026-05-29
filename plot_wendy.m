errs = vecnorm(w_its - w_true,2,1)/max(norm(w_true),eps);
err_windy=errs(end);

%%% compute conf intervals
w_hat = w_plot(w_plot~=0);
w_true_hat = w_true(w_plot~=0);
c0 = 0.05; c_maxmin = false;
[conf_int,c,confbounds,xflip,coverage] = get_conf_ints(w_hat,w_true_hat,CovW,c0,c_maxmin);

%%% print results
disp(['err WENDy:',num2str(err_windy)])
fprintf('coverage=')
fprintf('%1.2e ',coverage);
fprintf('\n')

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
% set(gca,'Ytick',10.^(floor(min(log10(dat))):1:ceil(max(log10(dat)))))

%%% confidence intervals
subplot(3,2,3)
% h=fill(xflip,confbounds,'red','FaceColor',[1 0.8 0.8],'edgecolor','none'); hold on;
% try
%     h1=plot(1:length(w_hat),w_hat,'ro',1:length(w_hat),w_true(w_plot~=0),'bx');
%     legend(h1(1:2),{'w','w^*'})
% catch
%     h1=plot(1:length(w_hat),w_hat,'ro');
%     legend('w')
% end
% hold off
for j=1:length(w_hat)
    gg = 0.3;
    h0 = fill([j-gg j+gg j+gg j-gg],[w_hat(j)-conf_int(j) w_hat(j)-conf_int(j)  w_hat(j)+conf_int(j) w_hat(j)+conf_int(j)],...
        'w','linewidth',1.5);
    hold on
    line([j-gg j+gg],[w_hat(j) w_hat(j)],'color','m','linewidth',3)
end
try
    h1=plot(1:length(w_hat),w_true_hat,'bx','linewidth',3,'markersize',6);
    h2=plot(1:length(w_hat),w_its(:,1),'gx','linewidth',3,'markersize',6);
    legend([h0;h1;h2],{[num2str((1-c)*100),'% CI'],'true val.','OLS'},'location','best')
catch
    legend(h0,{[num2str((1-c)*100),'% CI']},'location','best')
end
hold off
grid on

%%% p-values
if ~isempty(res_0)
    subplot(3,2,4)
    pvals = arrayfunvec(res,@(v)outn(@swtest,v,2),1);
    semilogy(pvals,'o-')
    grid on;
    title(['p-val=',num2str(pvals(end))])
    xlabel('iter')
    legend('p-val')

    %%% GLS residual: WENDy final to true
    subplot(3,2,5)
    plot(res(:,end)); 
    title('C^{-1/2}(b-G*w)')
    
    subplot(3,2,6)
    plot(res_0(:,1))
    title('b-G*w')
end