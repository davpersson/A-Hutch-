%------------- DESCRIPTION -------------

%This script reproduces Figure 1.  

%---------------------------------------

clc
clear

addpath('other')
addpath('results')

filename = 'results/minima';

rng(0)

n = 1000;
matrix_size = n;
c = 2;
n_max = 50;

D = diag((1:n).^(-c));
[Q,~] = qr(randn(n));
A = Q*D*Q';
fronorm2 = norm(A,'fro')^2;


epsilon = 0.05*trace(A);
delta = 0.01;
C = 4*log(2/delta)/epsilon^2;

Afun = @(X) A*X;

%First iteration
y = Afun(randn(matrix_size,1));
q = y/norm(y);
Q = q;

x = Afun(q);
t = q'*x;
c = t^2;
iteration = 1;
b = norm(x)^2;
fnc(iteration) = 2*iteration+C*(fronorm2+c - 2*b);

%Remaining iterations
while iteration < n_max
    
    %Get new column in Q
    y = Afun(randn(matrix_size,1));
    qt = y-Q*(Q'*y); 
    qt = qt/norm(qt); qt = qt - Q*(Q'*qt); q = qt/norm(qt);
    Q = [Q q];
    x = Afun(q);
    
    %Update recursion
    b = b + norm(x)^2;
    t = q'*x;
    c = c + 2*norm(Q(:,1:iteration)'*x)^2 + t^2;
    
    %Update function
    iteration = iteration + 1;
    fnc(iteration) = 2*iteration+C*(fronorm2+c - 2*b);
    
end

figure('Renderer', 'painters', 'Position', [10 10 900 300])
plot(fnc,'k-*','LineWidth',2)
ylabel('$$m(r)$$','interpreter','latex')
xlabel('$$r$$','interpreter','latex')
set(gca,'FontSize',18)

print(filename,'-depsc')