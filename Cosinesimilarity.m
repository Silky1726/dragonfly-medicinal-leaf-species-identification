function CosSim = CosineSimilarity(Data1,Data2)
if isvector(Data1)==0 || isvector(Data2)==0
    error('x and y have to be vectors!')
end
if length(Data1)~=length(Data2)
    error('x and y have to be same length!')
end
CosSim=dot(Data1,Data2)/(norm(Data1)*norm(Data2));
if isnan(CosSim) || isinf(CosSim)
   CosSim=1; 
end
end
