function C = matmat(n,Afun,B,f,n_it)
%Computes C = f(A)*B by n_it iterations of Lanczos

%m = size(A,1); 
m = size(B,2);
C = zeros(n,m);

for col = 1:m
    
    C(:,col) = matvec(n,Afun,B(:,col),f,n_it);
    
end

end