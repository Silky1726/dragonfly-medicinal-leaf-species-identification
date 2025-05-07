function [Con_Matrix, Result, RefereceResult] = ClassificationScore(Actual, Predict)
[rows,cols]=size(Actual);
Actual=Actual(:);
Predict=Predict(:);

%[X,Y,T,AUC] = perfcurve(Actual,Predict,[1 2 3]);
if length(Actual)~=length(Predict)
   error('Input have different lengths')
end
UActual=unique(Actual);
UPredict=unique(Predict);
UActual=1:2;
UPredict=1:2;
condition=length(UActual)==length(UPredict);
if ~condition
   error('Class List is not same in given inputs')
end
condition=(sum(UActual==UPredict)==length(UActual));
if ~condition
   error('Class List in given inputs are different')
end

Class_List=UActual;
disp('Total Class in Data...')
disp(Class_List)
fprintf('Total Instance = %d\n',cols);

N_Class=length(UActual);
Con_Matrix=zeros(N_Class);
PredictClass=cell(1,N_Class);
Class_Ref=cell(N_Class,1);
Row_name=cell(1,N_Class);
%Calculate conufsion for all classes
for i=1:N_Class
    Class_Ref{i,1}=['Class ',num2str(i),' = ',num2str(Class_List(i))];
    for j=1:N_Class
        val=(Actual==Class_List(i)) & (Predict==Class_List(j));
        Con_Matrix(i,j)=sum(val);
        PredictClass{i,j}=sum(val);
    end
    Row_name{i}=['ActualClass',num2str(i)];
    disp(Class_Ref{i})
end

C_Matrix_Table=cell2table(PredictClass);
C_Matrix_Table.Properties.RowNames=Row_name;
disp('Confusion Matrix')
disp(C_Matrix_Table)
[Result,RefereceResult]=CalParam(Con_Matrix);
RefereceResult.Class=Class_Ref;
Param=struct2table(RefereceResult);
disp(Param)
% disp('Overall Parameters')
%disp(Result)
end

function  [Result,RefereceResult] = CalParam(Con_Matrix)
[row,~]=size(Con_Matrix);
N_Class=row;
TP=zeros(1,N_Class);
FN=zeros(1,N_Class);
FP=zeros(1,N_Class);
TN=zeros(1,N_Class);
for i=1:N_Class
    TP(i)=Con_Matrix(i,i);
    FN(i)=sum(Con_Matrix(i,:))-Con_Matrix(i,i);
    FP(i)=sum(Con_Matrix(:,i))-Con_Matrix(i,i);
    TN(i)=sum(Con_Matrix(:))-TP(i)-FP(i)-FN(i);
end
P=TP+FN;
N=FP+TN;
Accuracy=(TP+TN)./(P+N);
Recall=TP./(TP+TN);
Precision=TP./(TP+FP);
Fmeasure=(2*Precision.*Recall)./(Precision+Recall);

Result.Accuracy=mean(Accuracy)*100;
Result.Precision=mean(Precision);
Result.Recall=mean(Recall);
Result.Fmeasure=mean(Fmeasure);

RefereceResult.Accuracy=100*Accuracy';
RefereceResult.Precision=Precision';
RefereceResult.Recall=Recall';
RefereceResult.Fmeasure=Fmeasure';
RefereceResult.TP=TP';
RefereceResult.FP=FP';
RefereceResult.FN=FN';
RefereceResult.TN=TN';
end
