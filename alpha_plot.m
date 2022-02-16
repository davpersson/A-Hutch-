%------------- DESCRIPTION -------------

%This script reproduces Figure 2.  

%---------------------------------------

clc
clear

addpath('A-Hutch++')
addpath('results')
filename = 'results/alpha';

k_list = 1:1000;
delta_list = [0.1,0.05,0.01];
alpha_list = zeros(length(delta_list),length(k_list));
legends = {'$\delta = 0.1$','$\delta = 0.05$','$\delta = 0.01$'};
lines = {'b-','r-','g-'};

for k = k_list
    
    inner_iteration = 0;
    for delta = delta_list
        
        inner_iteration = inner_iteration + 1;
        alpha_list(inner_iteration,k) = supfind(k,delta);
        
    end
    
end

%plot(k_list,alpha_list(1,:),lines{1},'LineWidth',3)
figure('Renderer', 'painters', 'Position', [10 10 900 400])
semilogx(k_list,alpha_list(1,:),lines{1},'LineWidth',3)
hold on
for index = 2:length(delta_list)
    
    %plot(k_list,alpha_list(index,:),lines{index},'LineWidth',3)
    semilogx(k_list,alpha_list(index,:),lines{index},'LineWidth',3)
    
end

legend(legends,'interpreter','latex','Location','best')
xlabel('$k$','interpreter','latex')
ylabel('$\alpha_k$','interpreter','latex')
set(gca,'FontSize',20)

print(filename,'-depsc')
hold off
