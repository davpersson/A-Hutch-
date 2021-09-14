function [T, Q] = lanczos(n, Afun, n_it, u)
%Code from Alice
    %n = size(A, 1);
    if (n_it > n)
        n_it = n;
    end
   

    x = u/norm(u);
    Q = [];
    xold = zeros(n,1);
    beta = 0;

    for j=1:n_it
        y = Afun(x);
        alpha = x'*y;
        Q = [Q, x];
        r = y - alpha*x - beta*xold;
        % Reorthogonalization:
        r = r - Q*(Q'*r);

        % Gauss quadrature rule
        if (j == 1)
            T = alpha;
        else
            T(j,j) = alpha;
            T(j-1,j) = beta;
            T(j,j-1) = beta;
        end

        % Update stuff
        beta = norm(r);
        xold = x;

        if (beta < 1e-10)
            return;
        end

        x = r/beta;
    end
end