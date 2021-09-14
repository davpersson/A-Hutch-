%------------- DESCRIPTION -------------

%To reproduce the numerical experiments from Section 3.4, first uncomment
%the corresponding block of code under "MATRIX CHOICE". Then run the
%script. The script will store the results in a file "test". To see the
%plots run script "plotter_Nystrom++.m". 

%---------------------------------------

clc
clear

addpath('Nystrom++')
addpath('other')
addpath('results')

filename = 'results/test';
rng(0)

%------------- MATRIX CHOICE -------------

%Uncomment one of the following blocks.

%----------------------------------------

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 8a)
% n = 5000;
% c = 0.1;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 8b)
% n = 5000;
% c = 0.5;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 8c)
 n = 5000;
 c = 2;
 D = (1:n).^(-c);
 [Q,~] = qr(randn(n));
 A = Q*diag(D)*Q';
 tr = trace(A);
 Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 8d)
% n = 5000;
% c = 5;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Estrada index (Figure 9)
% B = sparse(create_roget_mat()); 
% n_it = 35;
% n = size(B,1);
% f = @(X) expm(X);
% tr = trace(f(B));
% Bfun = @(X) B*X;
% Afun = @(X) matmat(n,Bfun,X,f,n_it);

%Inverse of tridiag(-1,4,-1) (Figure 10a)
% A = sparse(4*eye(10000) - diag(ones(9999,1),1) - diag(ones(9999,1),-1));
% tr = trace(inv(A));
% n = size(A,1);
% Afun = @(X) A\X;

%Inverse of discretization of Poisson's equation (Figure 10b)
% A = gallery('poisson',100);
% tr = trace(inv(A));
% n = size(A,1);
% Afun = @(X) A\X;

%Synthetic matrix with exponentially decaying eigenvalues/singular values.
%(Figure 11a)
% n = 1000;
% s = 10;
% D = diag(exp(-(1:n)/s));
% [Q,~] = qr(randn(n));
% A = Q*D*Q';
% tr = trace(A);
% Afun = @(x) A*x;

%Synthetic matrix with exponentially decaying eigenvalues/singular values.
%(Figure 11b)
% n = 1000;
% s = 100;
% D = diag(exp(-(1:n)/s));
% [Q,~] = qr(randn(n));
% A = Q*D*Q';
% tr = trace(A);
% Afun = @(x) A*x;

%------------- NUMERICAL TESTS -------------

%-------------------------------------------

matvecs_list = 12:48:1000; %List of number of matvecs
inner_repeats = 100; %Number of repeats to get a sample mean of the error
iteration = 0; %Counter

%Initiate vectors to save results
trace_outputs = zeros(length(matvecs_list),inner_repeats,3);

tic
for matvecs = matvecs_list
    
    %Update counter
    iteration = iteration + 1;
    disp(length(matvecs_list)-iteration)
    
    %Repeat experiment repeats times
    for k = 1:inner_repeats
        
        %Run experiment
        trest_npp = nystrompp(n,Afun,matvecs);
        trest_hpp = hutchpp(n,Afun,matvecs);
        trest_nahpp = nahpp(n,Afun,matvecs);
        
        %Save results
        trace_outputs(iteration,k,1) = trest_npp;
        trace_outputs(iteration,k,2) = trest_hpp;
        trace_outputs(iteration,k,3) = trest_nahpp;
        
    end
    
    
end
toc

save(filename,'tr','matvecs_list','trace_outputs');

beep
