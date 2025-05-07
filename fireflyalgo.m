function [BestSelection,Brighterdata]=fireflyalgo(Data,nOfSelection,nOfFireFlies,itrationMax,fun)
alphaval=0.5;% step size ,
gamaval=1;% light absorbance 
BrightnessInit=zeros(1,nOfFireFlies); %  Light intensity 
upperBound=size(Data,2);
Boundary=[1,upperBound];
FireflyPositions=zeros(nOfFireFlies,nOfSelection);
for i=1:nOfFireFlies
    tempVar=randperm(Boundary(2));
    FireflyPositions(i,:)=tempVar(1:nOfSelection);
end
%% Optimization
for itr=1:itrationMax
    distance_r_=dist(FireflyPositions,FireflyPositions');
    distance_r_=distance_r_/mean(distance_r_(:));
    % Calculate Brightness Io
    for k=1:nOfFireFlies
        BrightnessInit(k)=fun(Data(:,FireflyPositions(k,:)));
    end
    % Get Best 
    [Brighterdata,bestBrigthnessfilefly]=min(BrightnessInit);
    bestBrigthnessfilefly=bestBrigthnessfilefly(1);
    Brighterdata=abs(Brighterdata(1)/1000);
    % Update X(FireflyPositions)
    for i=1:nOfFireFlies
        if i~=bestBrigthnessfilefly
            A=(BrightnessInit(i) * exp(-gamaval*(distance_r_(bestBrigthnessfilefly,i))^2) * (FireflyPositions(i,:)-FireflyPositions(bestBrigthnessfilefly,:)) )...
                + (alphaval*(rand(1,nOfSelection)-0.5));
            FireflyPositions(i,:)=abs( FireflyPositions(i,:) +(A) );
        end
    end
    % Ranging X(FireflyPositions)
    FireflyPositions=round(FireflyPositions);
    FireflyPositions(FireflyPositions<Boundary(1))=Boundary(1);
    FireflyPositions(FireflyPositions>Boundary(2))=Boundary(2);
    for i=1:nOfFireFlies
        temp=(FireflyPositions(i,:));
        if length(unique(temp))~=nOfSelection
            tempVar=randperm(Boundary(2));
            FireflyPositions(i,:)=tempVar(1:nOfSelection);
        end
    end
    fprintf('itration %d best Atrraction value is %f\n',itr,Brighterdata)
end
BestSelection=mean2(FireflyPositions(bestBrigthnessfilefly,:))/(numel(FireflyPositions(bestBrigthnessfilefly,:)));
BestSelection=abs((BestSelection)/1e3);