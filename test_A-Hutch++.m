%------------- DESCRIPTION -------------

%To reproduce the numerical experiments from Section 2.3, first uncomment
%the corresponding block of code under "MATRIX CHOICE". Then run the
%script. The script will store the results in a file "test". To see the
%plots run script "plotter_A-Hutch++.m". 

%---------------------------------------

clc
clear

addpath('A-Hutch++')
addpath('other')
addpath('results')

filename = 'results/test';
rng(0)

%------------- MATRIX CHOICE -------------

%Uncomment one of the following blocks.

%----------------------------------------

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 3a)
n = 5000;
c = 0.1;
D = (1:n).^(-c);
[Q,~] = qr(randn(n));
A = Q*diag(D)*Q';
tr = trace(A);
Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 3b)
% n = 5000;
% c = 0.5;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 3c)
% n = 5000;
% c = 2;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Synthetic matrix with algebraically decaying eigenvalues/singular values.
%(Figure 3d)
% n = 5000;
% c = 5;
% D = (1:n).^(-c);
% [Q,~] = qr(randn(n));
% A = Q*diag(D)*Q';
% tr = trace(A);
% Afun = @(X) A*X;

%Triangle counting with Wikipedia vote network (Figure 4a)
% A = readWikiNet(1);
% tr = trace(A^3);
% n = size(A,1);
% Afun = @(X) A*(A*(A*X));

%Triangle counting with arXiv GR-QC network (Figure 4b)
% A = load('other/ca-GrQc.mat'); A = A.Problem; A = A.A;
% tr = trace(A^3);
% n = size(A,1);
% Afun = @(X) A*(A*(A*X));

%Estrada index (Figure 5)
% B = sparse(create_roget_mat()); 
% n_it = 35;
% n = size(B,1);
% f = @(X) expm(X);
% tr = trace(f(B));
% Bfun = @(X) B*X;
% Afun = @(X) matmat(n,Bfun,X,f,n_it);

%Log-determinant of matrix with eigenvalue gap (Figure 6a)
% n = 5000;
% n_it = 25;
% h = 10; 
% l = 1;
% X = sprand(n, 300, 0.025);
% D = spdiags([h*ones(40, 1)./(((1:40).^2)'); l * ones(260, 1)./(((41:300).^2)')], 0, 300, 300);
% Bfun = @(x) x + X * (D * (X' * x));
% f = @(Y) logm(Y);
% Afun = @(Y) matmat(n,Bfun,Y,f,n_it);
% tr = trace(logm(eye(300) + full(D * (X' * X))));

%Log-determinant of matrix Thermomech_TC (Figure 6b)
% A = load('thermomech_TC.mat'); A = A.Problem; A = A.A;
% n = size(A,1);
% n_it = 35;
% Bfun = @(x) A*x;
% f = @(Y) logm(Y);
% Afun = @(Y) matmat(n,Bfun,Y,f,n_it);
% tr = -546786.561681857;

%Inverse of tridiag(-1,4,-1) (Figure 7a)
% A = sparse(4*eye(10000) - diag(ones(9999,1),1) - diag(ones(9999,1),-1));
% tr = trace(inv(A));
% n = size(A,1);
% Afun = @(X) A\X;

%Inverse of discretization of Poisson's equation (Figure 7b)
% A = gallery('poisson',100);
% tr = trace(inv(A));
% n = size(A,1);
% Afun = @(X) A\X;

%------------- NUMERICAL TESTS -------------

%-------------------------------------------

%Parameters
repeats = 100;
tolerance_list = 2.^(-(1:9));
delta = 0.05;

%Allocate space for experiments
outputs_adap_hpp = zeros(length(tolerance_list),repeats);
outputs_hpp = zeros(length(tolerance_list),repeats);
outputs_adap_npp = zeros(length(tolerance_list),repeats);
outputs_npp = zeros(length(tolerance_list),repeats);

total_matvecs_hpp = zeros(length(tolerance_list),repeats);
lowrank_matvecs_hpp = zeros(length(tolerance_list),repeats);
trest_matvecs_hpp = zeros(length(tolerance_list),repeats);
total_matvecs_npp = zeros(length(tolerance_list),repeats);
lowrank_matvecs_npp = zeros(length(tolerance_list),repeats);
trest_matvecs_npp = zeros(length(tolerance_list),repeats);

outputs_hpp = zeros(length(tolerance_list),repeats);
outputs_npp = zeros(length(tolerance_list),repeats);

%Running tests and save results
outer_iteration = 0;
tic
for tolerance = tolerance_list
    
    outer_iteration = outer_iteration + 1
    
    for inner_iteration = 1:repeats
        %inner_iteration
        [outputs_adap_hpp(outer_iteration,inner_iteration),...
            total_matvecs_hpp(outer_iteration,inner_iteration),...
            lowrank_matvecs_hpp(outer_iteration,inner_iteration),...
            trest_matvecs_hpp(outer_iteration,inner_iteration)] = adap_hpp(n,Afun,tolerance*abs(tr),delta);
        
        outputs_hpp(outer_iteration,inner_iteration) = hutchpp(n,Afun,...
            total_matvecs_hpp(outer_iteration,inner_iteration));
        
    end
    
end
toc

save(filename,'tr','tolerance_list','delta',...
    'outputs_adap_hpp','outputs_hpp',...
    'total_matvecs_hpp','lowrank_matvecs_hpp',...
    'trest_matvecs_hpp')

beep