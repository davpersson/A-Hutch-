function [trest,total_matvecs,lowrank_matvecs,trest_matvecs] = block_adap_hpp(matrix_size,Afun,epsilon,delta,block_size)
%INPUTS
%A: Square matrix
%delta: Failure probability. Algorithm succeeds wpral 1-2*delta
%epsilon: Error tolerance

%OUTPUTS:
%trest: Approximation to trace(A)
%matvecs: Total number of matvecs used
%rsvd: Total number of matvecs used for randomized svd
%M: Total number of matvecs used for Hutchinson estimation

%Set up parameters
C = (4*log(2/delta)/epsilon^2); %Constant for number of matvecs

flag = 1;
fnc = [];

%First iteration
AOmega = Afun(randn(matrix_size,block_size));
Q = orth(AOmega);

AQ = Afun(Q);
QtAQ = Q'*AQ;
t = trace(QtAQ); 
c = norm(QtAQ,'fro'); %t^2;
trest1 = t;
iteration = block_size;
b = norm(AQ,'fro')^2;
fnc(iteration) = 2*iteration+C*(c - 2*b);

%Remaining iterations
while flag
    
    %Get new column in Q
    AOmega = Afun(randn(matrix_size,block_size));
    qt = AOmega-Q*(Q'*AOmega); 
    
    if norm(qt,'fro')/sqrt(block_size) < 1e-10
        
        lowrank_matvecs = 2*iteration;
        trest_matvecs = 0;
        total_matvecs = lowrank_matvecs + trest_matvecs;
        trest = trest1;
        return
        
        
    end
    
    qt = orth(qt); qt = qt - Q*(Q'*qt); q = orth(qt);
    Q = [Q q];
    Aq = Afun(q);
    
    %Update recursion
    b = b + norm(Aq,'fro')^2;
    qtAq = q'*Aq;
    t = trace(qtAq);
    trest1 = trest1 + t;
    c = c + 2*norm(Q(:,1:iteration)'*Aq,'fro')^2 + norm(qtAq,'fro')^2;
    
    %Update function
    iteration = iteration + block_size;
    fnc(iteration) = 2*iteration+C*(c - 2*b);
    
    %Check if iteration should be stopped
    if (iteration > block_size) && (fnc(iteration-block_size) < fnc(iteration))
        
        flag = 0;
        
    end
    
end

lowrank_matvecs = 2*iteration;

%Combine Hutchinson and Frobenius norm estimation
iteration = 0;
flag = 1;
t = 0;
trest2 = [];

while flag
    
    psi = randn(matrix_size,block_size);
    y = psi - Q*(Q'*psi);
    y = Afun(y);
    y = y - Q*(Q'*y);
    trest2 = [trest2 trace(psi'*y)/block_size];
    t = t + trace(y'*y);
    estFrob = t/(iteration+block_size);
    alpha = supfind(iteration+block_size,delta);
    M = ceil(C*estFrob/alpha);
    
    if iteration + block_size > M
        
        flag = 0;
        
    end
    
    iteration = iteration + block_size;
    
end

%Perform Hutchinson estimation
trest2 = mean(trest2);
trest_matvecs = iteration;


%Return outputs
trest = trest1 + trest2;
total_matvecs = lowrank_matvecs + trest_matvecs; 
end