function [keypoints]=SIFT(I,J)
I=imresize(I,[256 256]);
J=imresize(J,[256 256]);
I=im2double(I);
octavemtrx=[];
k2=0;
[m,n]=size(I);
I(m:m+4,n:n+4)=0;
for k1=0:3                   
k=sqrt(2);
sigma=(k^(k1+(2*k2)));
for x=-2:2 
for y=-2:2
h(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma))); 
end
end
for i=1:m
for j=1:n
tempimg=I(i:i+4,j:j+4)'.*h;
convmtrx(i,j)=sum(sum(tempimg));
end
end
octavemtrx=[octavemtrx convmtrx];
end
pause(round(5*rand));
diffmtrx1=octavemtrx(1:m,m+1:2*m)-octavemtrx(1:m,((2*m)+1):3*m);
diffmtrx2=octavemtrx(1:m,((2*m)+1):3*m)-octavemtrx(1:m,((3*m)+1):4*m);
key=[];
x1=0;
y1=0;
for i=2:m-1
for j=2:m-1
if (((diffmtrx1(i,j)>diffmtrx1(i-1,j))&&(diffmtrx1(i,j)>diffmtrx1(i+1,j))....
&&(diffmtrx1(i,j)>diffmtrx1(i,j-1))&&(diffmtrx1(i,j)>diffmtrx1(i+1,j+1))....
&&(diffmtrx1(i,j)>diffmtrx1(i-1,j-1))&&(diffmtrx1(i,j)>diffmtrx1(i-1,j+1))....
&&(diffmtrx1(i,j)>diffmtrx1(i+1,j-1))&&(diffmtrx1(i,j)>diffmtrx1(i,j+1))))
x1=x1+1;
else
continue;
end
if x1>0
if((diffmtrx1(i,j)>diffmtrx2(i,j))&&(diffmtrx1(i,j)>diffmtrx2(i-1,j))....
&&(diffmtrx1(i,j)>diffmtrx2(i+1,j))&&(diffmtrx1(i,j)>diffmtrx2(i,j-1))....
&&(diffmtrx1(i,j)>diffmtrx2(i+1,j+1))&&(diffmtrx1(i,j)>diffmtrx2(i-1,j-1))....
&&(diffmtrx1(i,j)>diffmtrx2(i-1,j+1))&&(diffmtrx1(i,j)>diffmtrx2(i+1,j-1))&&(diffmtrx1(i,j)>diffmtrx2(i,j+1)))
y1=y1+1;
else
continue;
end 
end
key(i,j)=diffmtrx1(i,j);
end
end
keypoints=key*255-(8*mean2(J));