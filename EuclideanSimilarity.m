function ESim = EuclideanSimilarity(Data1,Data2)
if isvector(Data1)==0 || isvector(Data2)==0
    error('x and y have to be vectors!')
end
if length(Data1)~=length(Data2)
    error('x and y have to be same length!')
end
ESim=sqrt(sum((Data1 - Data2) .^ 2));
if isnan(ESim) || isinf(ESim)
   ESim=1; 
end
end