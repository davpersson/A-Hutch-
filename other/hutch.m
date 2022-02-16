function trest = hutch(n,Afun,matvecs)

Omega = randn(n,matvecs);
trest = product_trace(Omega',Afun(Omega))/matvecs;

end

function t = product_trace(A,B)
%Computes trace(A*B) without computing all entries of A*B
t = sum(sum(A.*B',2));
end