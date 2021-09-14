function alpha = supfind(i,d)

f = @(a) gammainc(a*i/2,i/2);
alphalist = flip(0:0.01:1);

fa = 1;
index = 1;

while fa > d
    
    fa = f(alphalist(index));
    index = index + 1;
    
end

alpha = alphalist(index-1);

end