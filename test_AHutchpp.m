%------------- DESCRIPTION -------------

% This script reproduces the results in Figures 4 -- 8

%---------------------------------------

clc
clear

addpath('A-Hutch++')
addpath('other')
addpath('results')

rng(0)

for matrix = 1:12

    %------------- MATRIX CHOICE -------------
    
    if (matrix == 1)
        
        %Synthetic matrix with algebraically decaying eigenvalues/singular values.
        %(Figure 4a)
        n = 5000;
        c = 0.1;
        D = (1:n).^(-c);
        [Q,~] = qr(randn(n));
        A = Q*diag(D)*Q';
        tr = trace(A);
        Afun = @(X) A*X;
        filename = 'results/algebraic_c=01';
        hutchinson = true;
        
    elseif (matrix == 2)

        %Synthetic matrix with algebraically decaying eigenvalues/singular values.
        %(Figure 4b)
        n = 5000;
        c = 0.5;
        D = (1:n).^(-c);
        [Q,~] = qr(randn(n));
        A = Q*diag(D)*Q';
        tr = trace(A);
        Afun = @(X) A*X;
        filename = 'results/algebraic_c=05';
        hutchinson = true;
        
    elseif (matrix == 3)

        % Synthetic matrix with algebraically decaying eigenvalues/singular values.
        % (Figure 4c)
        n = 5000;
        c = 1;
        D = (1:n).^(-c);
        [Q,~] = qr(randn(n));
        A = Q*diag(D)*Q';
        tr = trace(A);
        Afun = @(X) A*X;
        filename = 'results/algebraic_c=1';
        hutchinson = false;
        
    elseif (matrix == 4)

        % Synthetic matrix with algebraically decaying eigenvalues/singular values.
        % (Not included in the paper)
        n = 5000;
        c = 2;
        D = (1:n).^(-c);
        [Q,~] = qr(randn(n));
        A = Q*diag(D)*Q';
        tr = trace(A);
        Afun = @(X) A*X;
        filename = 'results/algebraic_c=2';
        hutchinson = false;
        
    elseif (matrix == 5)

        %Synthetic matrix with algebraically decaying eigenvalues/singular values.
        %(Figure 4d)
        n = 5000;
        c = 3;
        D = (1:n).^(-c);
        [Q,~] = qr(randn(n));
        A = Q*diag(D)*Q';
        tr = trace(A);
        Afun = @(X) A*X;
        filename = 'results/algebraic_c=3';
        hutchinson = false;
        
    elseif (matrix == 6)

        %Triangle counting with Wikipedia vote network (Figure 5a)
        A = readWikiNet(1);
        tr = trace(A^3);
        n = size(A,1);
        Afun = @(X) A*(A*(A*X));
        filename = 'results/cube_wiki';
        hutchinson = false;
        
    elseif (matrix == 7)

        %Triangle counting with arXiv GR-QC network (Figure 5b)
        A = load('other/ca-GrQc.mat'); A = A.Problem; A = A.A;
        tr = trace(A^3);
        n = size(A,1);
        Afun = @(X) A*(A*(A*X));
        filename = 'results/cube_arxiv';
        hutchinson = false;
        
    elseif (matrix == 8)

        %Estrada index (Figure 6)
        B = sparse(create_roget_mat()); 
        n_it = 35;
        n = size(B,1);
        f = @(X) expm(X);
        tr = trace(f(B));
        Bfun = @(X) B*X;
        Afun = @(X) matmat(n,Bfun,X,f,n_it);
        filename = 'results/estrada';
        hutchinson = false;
        
    elseif (matrix == 9)

        %Log-determinant of matrix with eigenvalue gap (Figure 7a)
        n = 5000;
        n_it = 25;
        h = 10; 
        l = 1;
        X = sprand(n, 300, 0.025);
        D = spdiags([h*ones(40, 1)./(((1:40).^2)'); l * ones(260, 1)./(((41:300).^2)')], 0, 300, 300);
        Bfun = @(x) x + X * (D * (X' * x));
        f = @(Y) logm(Y);
        Afun = @(Y) matmat(n,Bfun,Y,f,n_it);
        tr = trace(logm(eye(300) + full(D * (X' * X))));
        filename = 'results/saibaba_ipsen';
        hutchinson = false;
        
    elseif (matrix == 10)

        %Log-determinant of matrix Thermomech_TC (Figure 7b)
        A = load('thermomech_TC.mat'); A = A.Problem; A = A.A;
        n = size(A,1);
        n_it = 35;
        Bfun = @(x) A*x;
        f = @(Y) logm(Y);
        Afun = @(Y) matmat(n,Bfun,Y,f,n_it);
        tr = -546786.561681857;
        filename = 'results/thermomec';
        hutchinson = true;
        
    elseif (matrix == 11)

        %Inverse of tridiag(-1,4,-1) (Figure 8a)
        A = sparse(4*eye(10000) - diag(ones(9999,1),1) - diag(ones(9999,1),-1));
        tr = trace(inv(A));
        n = size(A,1);
        Afun = @(X) A\X;
        filename = 'results/tridiaginv';
        hutchinson = true;
        
    elseif (matrix == 12)

        %Inverse of discretization of Poisson's equation (Figure 8b)
        A = gallery('poisson',100);
        tr = trace(inv(A));
        n = size(A,1);
        Afun = @(X) A\X;
        filename = 'results/poissoninv';
        hutchinson = false;
        
    end

    %------------- NUMERICAL TESTS -------------

    %-------------------------------------------

    %Parameters
    repeats = 100;
    %tolerance_list = 2.^(-(1:9));
    tolerance_list = 2.^(-(2:10));
    
    if matrix == 10
        
        tolerance_list = 2.^(-(3:11));
        
    end
    
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
    outputs_hutch = zeros(length(tolerance_list),repeats);

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
            
            if hutchinson == true
                
                outputs_hutch(outer_iteration,inner_iteration) = hutch(n,Afun,total_matvecs_hpp(outer_iteration,inner_iteration));
                
            end

        end

    end
    toc

    save(filename,'tr','tolerance_list','delta','hutchinson',...
        'outputs_adap_hpp','outputs_hpp','outputs_hutch',...
        'total_matvecs_hpp','lowrank_matvecs_hpp',...
        'trest_matvecs_hpp')
    
    plotter_AHutchpp(filename);
end

% beep