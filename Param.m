function ErrorRate = Param(A,B)
% PSNR (Peak Signal to noise ratio)
try
A=imresize(A,[256 256]);
B=imresize(B,[256 256]);
end
try
A=rgb2gray(A);
end
try
B=rgb2gray(B);
end
if (size(A) ~= size(B))
   error('The size of the 2 matrix are unequal')
   psnr_Value = NaN;
   return; 
elseif (A == B)
   psnr_Value = Inf;
   return;   
else
    mseImage = (double(A) - double(B)) .^ 2;
    [rows, columns] = size(A);
    ErrorRate = (sum(mseImage(:)) / (rows * columns))/1e3;
end
end 