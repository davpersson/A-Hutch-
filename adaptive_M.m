%------------- DESCRIPTION -------------

%This script reproduces Figure 2.  

%---------------------------------------

clc
clear

addpath('A-Hutch++')
addpath('other')
addpath('results')

filename = 'results/adaptive_M';

rng(1)

%Parameters for the matrix
n = 1000;
c = 1.5;
matrix_size = n;

%Create matrix
[Q,~] = qr(randn(n));
D = diag((1:n).^(-c));
A = Q*D*Q';
Afun = @(X) A*X;

%Experiment parameters
epsilon = 0.01*trace(A);
delta = 0.05;

C = (4*log(2/delta)/epsilon^2); %Constant for number of matvecs

flag = 1;
fnc = [];

%First iteration
y = Afun(randn(matrix_size,1));
q = y/norm(y);
Q = q;

x = Afun(q);
t = q'*x;
c = t^2;
trest1 = t;
iteration = 1;
b = norm(x)^2;
fnc(iteration) = 2*iteration+C*(c - 2*b);

%Remaining iterations
while flag
    
    %Get new column in Q
    y = Afun(randn(matrix_size,1));
    qt = y-Q*(Q'*y); 
    qt = qt/norm(qt); qt = qt - Q*(Q'*qt); q = qt/norm(qt);
    Q = [Q q];
    x = Afun(q);
    
    %Update recursion
    b = b + norm(x)^2;
    t = q'*x;
    trest1 = trest1 + t;
    c = c + 2*norm(Q(:,1:iteration)'*x)^2 + t^2;
    
    %Update function
    iteration = iteration + 1;
    fnc(iteration) = 2*iteration+C*(c - 2*b);
    
    %Check if iteration should be stopped
    if (iteration > 2) && (fnc(iteration-1) < fnc(iteration)) && (fnc(iteration-2) < fnc(iteration-1))
        
        flag = 0;
        
    end
    
end

lowrank_matvecs = 2*iteration;

%Combine Hutchinson and Frobenius norm estimation
iteration = 0;
flag = 1;
flag2 = 1;
t = 0;
trest2 = [];
M_list = [];


while flag||flag2
    
    psi = randn(matrix_size,1);
    y = psi - Q*(Q'*psi);
    y = Afun(y);
    y = y - Q*(Q'*y);
    trest2 = [trest2 psi'*y];
    t = t + y'*y;
    estFrob = t/(iteration+1);
    alpha = supfind(iteration+1,delta); 
    %if iteration > 99
    %    alpha = 1;
    %end
    M = C*estFrob/alpha;
    
    if (iteration + 1 > M)&&flag
        
        flag = 0;
        
        iteration_stop = iteration + 1;
        
    end
    
    if (~flag)&&(iteration > 2*iteration_stop)
        
        flag2 = 0;
        
    end
    
    iteration = iteration + 1;
    
    M_list(iteration) = M;
    
end

P = eye(n)-Q*Q';

figure('Renderer', 'painters', 'Position', [10 10 900 300])
plot(4:iteration,M_list(4:end),'r','LineWidth',2)
%BreakPlot(2:iteration,M_list(2:end),200,400,'r','LineWidth',2)
hold on
plot(4:iteration,4:iteration,'b','LineWidth',2)
%BreakPlot(1:iteration,1:iteration,200,400,'b','LineWidth',2)
plot(4:iteration,C*(norm(P*A*P,'fro')^2)*ones(1,iteration-3),'k','LineWidth',2)
plot([iteration_stop iteration_stop],[0 iteration_stop],'k--','LineWidth',1)
%BreakPlot(1:iteration,C*(norm(P*A*P,'fro')^2)*ones(1,iteration),200,400,'k','LineWidth',2)
legend({'$$M_k$$','$$k$$',...
    '$$C(\varepsilon,\delta)\|\textbf{\emph A}_{\mathrm{rest}}\|_F^2$$'},...
    'interpreter','latex')
xlabel('$$k$$','interpreter','latex')
set(gca,'FontSize',18)
hold off

print(filename,'-depsc')