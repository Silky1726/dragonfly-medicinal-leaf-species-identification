clc;
clear all;
close all;
load allfeatures_updated_184_2_classes.mat
load gtvalues.mat
%labels=labels';
%[labels,k_cent]=kmeans(Allfeatures,2);

total=numel(labels);
 %GTv=zeros(2,total);
 GTv=zeros(2,total);
 for j=1:numel(labels)
    GTv(labels(j),j)=labels(j); 
 end
[rows,cols]=size(Allfeatures);
totalindexes=310;
af=[];
for i=1:rows
   for j=1:totalindexes
      af(i,j)=Allfeatures(i,j); 
   end
end
Allfeatures=[];
Allfeatures=af;
TrainingData=Allfeatures;
nttrain=GTv;
net=patternnet(10);
%net=newff(TrainingData',nttrain,41);
net.trainParam.epochs=100;
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;                                                                                                  
net.divideParam.testRatio = 15/100;
net=train(net,TrainingData',nttrain);
%score=abs(round(sim(net,TrainingData')));
%result=abs(score);
%resultplot=[];
%trainplot=[];
%[rows,cols]=size(nttrain);
%result(1,1)=2;


% Assuming net is already trained
score = sim(net, TrainingData');
[~, result] = max(score, [], 1); % result now contains class labels 1 or 2

total1=numel(result);
 %GTv=zeros(2,total);
 PTv=zeros(2,total1);
 for j=1:numel(result)
    PTv(result(j),j)=result(j); 
 end


resultplot=[];
trainplot=[];
[rows,cols]=size(nttrain);


Classification_Score(nttrain,PTv);
 
%  plotroc(nttrain,result);
 
 classes={'Apple','Blueberry'}; % add class name 
 figure(22)
 aucv=0;
 for i=1:2 % change here 2 to 3
 [X,Y,T,AUC] = perfcurve(nttrain(i,:),PTv(i,:),i); 
 aucv=aucv+AUC;
 hold on;
 plot(X,Y,'x-')
xlabel('False positive rate') 
ylabel('True positive rate')

 end
 legend(classes);
title(['AUC Value:',num2str(aucv/2)]); 

figure,
plotconfusion(nttrain,PTv);
save('neuralvalue','net');

