function y = matvec(n,Afun,x,f,n_it)

[T,U] = lanczos(n,Afun,n_it,x);
F = f(T);
y = norm(x)*U*F(:,1);

end