function plotter_Nystrompp(filename)

addpath('Nystrom++')
addpath('other')
addpath('results')

load(filename);

relative_errors_npp = sort(abs((trace_outputs(:,:,1)-tr)/tr),2);
relative_errors_hpp = sort(abs((trace_outputs(:,:,2)-tr)/tr),2);
relative_errors_nahpp = sort(abs((trace_outputs(:,:,3)-tr)/tr),2);

mean_error_npp = mean(abs((trace_outputs(:,:,1)-tr)/tr),2);
mean_error_hpp = mean(abs((trace_outputs(:,:,2)-tr)/tr),2);
mean_error_nahpp = mean(abs((trace_outputs(:,:,3)-tr)/tr),2);

matvecs = matvecs_list';

prctile_10_npp = relative_errors_npp(:,10);
prctile_90_npp = relative_errors_npp(:,90);
prctile_10_hpp = relative_errors_hpp(:,10);
prctile_90_hpp = relative_errors_hpp(:,90);
prctile_10_nahpp = relative_errors_nahpp(:,10);
prctile_90_nahpp = relative_errors_nahpp(:,90);


semilogy(matvecs',prctile_10_hpp','r--')
hold on
semilogy(matvecs',prctile_90_hpp','r--')
patch([matvecs' fliplr(matvecs')],[prctile_10_hpp'...
    fliplr(prctile_90_hpp')],'r')

semilogy(matvecs',prctile_10_npp','b--')
semilogy(matvecs',prctile_90_npp','b--')
patch([matvecs' fliplr(matvecs')],[prctile_10_npp'...
    fliplr(prctile_90_npp')],'b')


alpha(.1)

h(1)=semilogy(matvecs,mean_error_npp,'b-*','LineWidth',2);
h(2)=semilogy(matvecs,mean_error_hpp,'r-*','LineWidth',2);
h(3)=semilogy(matvecs,mean_error_nahpp,'g-*','LineWidth',2);

xlabel('Number of matrix-vector products','interpreter','latex')
ylabel('Relative error','interpreter','latex')
legend(h,{'Nystrom++','Hutch++','Single Pass Hutch++'},'interpreter','latex')

set(gca,'FontSize',18)
hold off

print(filename,'-depsc')