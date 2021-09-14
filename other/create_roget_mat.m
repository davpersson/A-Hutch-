function A = create_roget_mat()
%based on https://www.mathworks.com/matlabcentral/answers/169557-read-a-text-file-with-varying-number-of-colums

%Save rows in a cell
fid = fopen( 'roget.net' );
data_cell = textscan( fid, '%s%s%s%s%s%s%s%s%s', 'CollectOutput' ...
    , true, 'Delimiter', ':;' );
[~] = fclose( fid );
data_cell = data_cell{1};

%Create matrix
A = zeros(1022,1022);
rows = size(data_cell); rows = rows(1);

%Append rows. Last row is empty
for row = 1:rows
    
    indices = str2num(data_cell{row,1});
    A(indices(1),indices(2:end)) = 1;
    
end

A = A + A';
A(A == 2) = 1;


%if norm(A-A','fro') > 1e-10
    
%    errorStruct.message = 'Matrix is not symmetric';
%    errorStruct.identifier = 'create_roget_mat:noSymmetry';
%    error(errorStruct)
    
%end


end