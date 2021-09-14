function [trest] = hutchpp(n,Afun,matvecs)
%function [trest,approximation_error] = hutchpp(n,Afun,matvecs)

if mod(matvecs,3) == 1
    
    matvecs = matvecs-1;
    
elseif mod(matvecs,3) == 2
    
    matvecs = matvecs + 1;
    
end

Q = orth(Afun(randn(n,matvecs/3)));
Psi = randn(n,matvecs/3);
Psi = Psi - Q*(Q'*Psi);

tr1 = product_trace(Q',Afun(Q));
tr2 = 3*(product_trace(Psi',Afun(Psi)))/matvecs;
trest = tr1 + tr2;

%To check behaviour of trace error. Uncomment if necessary. 
%AQ = Afun(Q);
%approximation_error = 2*norm(AQ,'fro')^2-norm(Q'*AQ,'fro')^2;

end

function t = product_trace(A,B)
%Computes trace(A*B) without computing all entries of A*B
t = sum(sum(A.*B',2));
end