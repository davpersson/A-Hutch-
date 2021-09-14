function A = readWikiNet(sym)
% Creates the adjacency matrix of the graph

fileID = fopen('wikiNet.txt', 'r');
v = fscanf(fileID, '%d');

u = unique(v);
map = sparse(max(v), 1);
map(u) = 1:7115;

A = sparse(7115, 7115);

for i = 1:2:length(v)-1
    A(map(v(i)), map(v(i+1))) = 1;
    if (sym == 1)
        A(map(v(i+1)), map(v(i))) = 1;
    end
end