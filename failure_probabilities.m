%------------- DESCRIPTION -------------

%This script reproduces Table 1. 

%---------------------------------------

clc
clear

addpath('A-Hutch++')
addpath('other')
addpath('results')

filename = 'results/failure_probabilities';

rng(0)

n = 5000;
repeats = 100000;

%Initialize tensor ton save results
table_tensor = zeros(3,3,5);

%Initialize tube index for table_tensor
tube = 0;
tic
for c = [0.1,0.5,1,2,3]
    
    tube = tube + 1;
    
    %Initialize column index for table_tensor
    column = 0;
    
    %Create matrix
    D = (1:n).^(-c);
    D = sparse(diag(D));
    tr = trace(D);
    Afun = @(X) D*X;
    
    for delta = [0.1,0.05,0.01]
        
        column = column + 1;
        
        %Initialize row index for table_tensor
        row = 0;
        
        for epsilon = [0.1,0.01,0.005]
            
            row = row + 1;
            
            %Create list to save outputs from each iteration
            outputs = zeros(1,repeats);
            
            %Run adap_hpp repeats number times
            
            for iteration = 1:repeats

                if mod(iteration,10000)==0

                    disp(iteration)

                end

                outputs(iteration) = adap_hpp(n,Afun,epsilon*tr,delta);

            end
            
            
            table_tensor(row,column,tube) = mean((abs(outputs-tr)/tr) > epsilon);
            fprintf("c = %f, epsilon/tr(A) = %f, delta = %f, Estimated failure probability = %f\n\n"...
                ,c,epsilon,delta,table_tensor(row,column,tube))
            
        end
        
    end
    
end
toc

save(filename,'table_tensor')