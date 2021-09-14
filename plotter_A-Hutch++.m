clc
clear

addpath('A-Hutch++')
addpath('other')
addpath('results')

filename = 'results/test';
load(filename);

relative_errors_adap_hpp = sort(abs((outputs_adap_hpp-tr)/tr),2);
relative_errors_hpp = sort(abs((outputs_hpp-tr)/tr),2);

mean_error_adap_hpp = mean(abs((outputs_adap_hpp-tr)/tr),2);
mean_error_hpp = mean(abs((outputs_hpp-tr)/tr),2);
matvecs = mean(total_matvecs_hpp,2);

prctile_10_adap_hpp = relative_errors_adap_hpp(:,10);
prctile_90_adap_hpp = relative_errors_adap_hpp(:,90);
prctile_10_hpp = relative_errors_hpp(:,10);
prctile_90_hpp = relative_errors_hpp(:,90);


loglog(matvecs',prctile_10_hpp','r--')
hold on
loglog(matvecs',prctile_90_hpp','r--')
patch([matvecs' fliplr(matvecs')],[prctile_10_hpp'...
    fliplr(prctile_90_hpp')],'r')


loglog(matvecs',prctile_10_adap_hpp','b--')
loglog(matvecs',prctile_90_adap_hpp','b--')
patch([matvecs' fliplr(matvecs')],[prctile_10_adap_hpp'...
    fliplr(prctile_90_adap_hpp')],'b')


alpha(.1)

h(1)=loglog(matvecs,mean_error_adap_hpp,'b-*','LineWidth',3);
h(2)=loglog(matvecs,mean_error_hpp,'r-*','LineWidth',3);
h(3)=loglog(matvecs,tolerance_list,'k--*','LineWidth',3);

xlabel('Number of matrix-vector multiplies','interpreter','latex')
ylabel('Relative error','interpreter','latex')
legend(h,{'A-Hutch++','Hutch++','Input tolerance'},'interpreter','latex')

set(gca,'FontSize',24)
hold off

T = table(mean(total_matvecs_hpp,2),mean(lowrank_matvecs_hpp,2),mean(trest_matvecs_hpp,2),...
    mean(trest_matvecs_hpp,2)./mean(total_matvecs_hpp,2),...
    'VariableNames',{'Total matvecs','Low rank','Trace est','Ratio:Hutchinson'});
T

print(filename,'-depsc')