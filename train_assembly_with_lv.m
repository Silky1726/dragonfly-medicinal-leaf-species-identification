clc;
clear all;
close all;
total_images=500
Allfeatures=[];
errorfound=0;
image_info=dir('Test data/*.JPG')
imgc=0;
image_index=[];
for img=1:total_images
    pathname=[image_info(img).name];
    fullpath=image_info.folder;
    flpath=[fullpath,'/',pathname];
    Timg=imread(flpath);
    figure(2)
    imshow(Timg);title(['Training Image',num2str(img),],'color','k','FontSize',12);
    Eimg=imadjust(Timg,stretchlim(Timg));

    imshow(Eimg)
    title('Pre-processed Image','color','k','FontSize',12);
    CTstruct=makecform('srgb2lab');
CTimg=applycform(Eimg,CTstruct);
ab=double(CTimg(:,:,2:3));
nrows=size(ab,1);
ncols=size(ab,2);
ab=reshape(ab,nrows*ncols,2);
nPart=2;
[Spart,Spartcen]=kmeans(ab,nPart,'distance','sqEuclidean','Replicates',2);
pixel_labels=reshape(Spart,nrows,ncols);
Simg=cell(1,2);
RGBlabel=repmat(pixel_labels,[1,1,3]);
for k=1:nPart
    colors=Eimg;
    colors(RGBlabel ~= k)=0;
    Simg{k}=colors;
end

    posimg=[mean2(Simg{1});mean2(Simg{2})];
    pos=find(posimg==min(posimg));
    PartSimg=Simg{pos};
    nOfSelection=5;
    DataPos=find(PartSimg>0);
    Data=PartSimg(DataPos);
    Data=Data(1:100,1)';
    nOfFireFlies=50;
    itrationMax=100;
    fitfun=@meandata;
    [FFAThresh,~]=fireflyalgo(Data,nOfSelection,nOfFireFlies,itrationMax,fitfun);
    Bimg=im2bw(rgb2gray(PartSimg),FFAThresh);
    Bimg=imfill(Bimg,'holes');
    figure(2);cla
    imshow(Bimg);title('Mask Image','color','k','FontSize',12);
    sts=regionprops(Bimg,'Eccentricity','Area');
    eccn=[sts.Eccentricity];
    ios=find(eccn);
    figure(3)
    imshow(Eimg);
    title('Selected Region','color','k','FontSize',12);
    hold on;
    boundaries=bwboundaries(Bimg);
    Marea=0;
    for k=1:length(ios)
    if sts(k).Area>Marea
        Marea=sts(k).Area;
        krec=k;
    end  
    end
    b = boundaries{krec};
    plot(b(:,2),b(:,1),'r','LineWidth',3);
    [rows, columns, pln] = size(Eimg);
    Segreg = poly2mask(b(:,2),b(:,1), rows, columns);
    Redimg=Eimg(:,:,1);
    Greenimg=Eimg(:,:,2);
    Blueimg=Eimg(:,:,3);
    FinalSimgR = double(Redimg).*double(imresize(Segreg,[size(Redimg,1) size(Redimg,2)]));
    FinalSimgG = double(Greenimg).*double(imresize(Segreg,[size(Greenimg,1) size(Greenimg,2)]));
    FinalSimgB = double(Blueimg).*double(imresize(Segreg,[size(Blueimg,1) size(Blueimg,2)]));
    FinalSimg = uint8(cat(3,FinalSimgR,FinalSimgG,FinalSimgB));
    figure(3)
    imshow(FinalSimg,[]);
    title('Segmented Image','color','k','FontSize',12);  
    %image_info1=dir('databases/GT Data/*.JPG');
    %pathname1=[image_info1(img).name];
    %fullpath1=image_info1.folder;
    %flpath1=[fullpath1,'/',pathname1];

    %TotalRegion=imread(flpath1);
    %TotalRegion=imresize(TotalRegion,[size(PartSimg,1) size(PartSimg,2)]);
    %if size(TotalRegion,3)>1
     %  A=im2bw(rgb2gray(TotalRegion));
    %else
     %  A=im2bw(TotalRegion); 
    %end
    %if size(FinalSimg,3)>1
     %  B=im2bw(rgb2gray(FinalSimg));
    %else
     %  B=im2bw(FinalSimg); 
    %end

    if size(FinalSimg,3)>1
    GRegimg=rgb2gray(FinalSimg);
    else
    GRegimg=FinalSimg;
    end
    Fpoints = detectBRISKFeatures(GRegimg);
    [features, validPoints] = extractFeatures(GRegimg,Fpoints);
    figure(3);imshow(uint8(Eimg)); hold on;
    plot(Fpoints.Location(:,1),Fpoints.Location(:,2),'ro','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','y');
    title('BRISK Points','color','k','FontSize',12);
   
   
    %ERBRISK=Param(features.Features,Eimg);
    allevaluations=[];
    errorfound=0;
    try
selected_indexes=[];
sind=0;
dragon_set=[];
dragon_scale=double(Fpoints.Scale);
dragon_points=double(Fpoints.Metric);
[rows,cols]=size(dragon_points);
dragon_set(1:rows,1)=dragon_points(:,1);
dragon_set(1:rows,2)=dragon_scale(:,1);
[dragon_index,dragon_global_best]=kmeans(dragon_set,2); % dividing data into two regions, outer and inner
for dr=1:numel(dragon_index)
   current_dragon_at=dragon_set(dr,:);
   dragon_home_index=dragon_index(dr);
   
   dragon_population= find(dragon_index==dragon_home_index); % these dragons belongs to the same group 
   
   %%%%% performing exploitation 
   dragons_to_group=round(numel(dragon_population)*.30); %taking 20 percent margin for grouping
   dragon_home_all_at=dragon_set(dragon_population,:);
   levy_flights=10;
   Lambda=1.5;
   acceptance_rejection_probablity=zeros(1,levy_flights); % probablity of being accepted 
   acceptance_rejection_order=zeros(1,levy_flights); % boolean , either true or false
   for lv=1:levy_flights
      dragon_grouping_indexes=round(rand(1,dragons_to_group)*dragons_to_group);
      dragon_grouping_indexes=findandreplace(dragon_grouping_indexes,0,dragon_population(1));
      dragon_grouped_at=dragon_home_all_at(dragon_grouping_indexes,:);
      
      %modified Lv
      for i = 1:size(dragon_grouped_at, 1)
        step = levy_flight(Lambda);
        dragon_grouped_at(i, :) = dragon_grouped_at(i, :) + step; % Update position with Levy step
      end
      
      
      dragon_grouped_at(dragons_to_group+1,1)=current_dragon_at(1,1);
      dragon_grouped_at(dragons_to_group+1,1)=current_dragon_at(1,2);
      global_best=dragon_global_best(dragon_home_index,:);
      cohesion_d1=0;% cosine similarity 
      cohesion_d2=0; % euclidean distance 
      
      for coh=1:size(dragon_grouped_at,1)
      cohesion_d1=cohesion_d1+Cosinesimilarity(dragon_grouped_at(coh,:),global_best);
      cohesion_d2=cohesion_d2+EuclideanSimilarity(dragon_grouped_at(coh,:),global_best);
      
      end
    
      %modified Lv
      [rows, ~] = size(dragon_grouped_at);
     cohesion_d1 = cohesion_d1 / rows;
     cohesion_d2 = cohesion_d2 / rows;
     cohesion = cohesion_d1 / cohesion_d2;


      %earlier withoutlv
      %[rows,cols]=size(dragon_grouped_at);
      %cohesion_d1=cohesion_d1/rows;
      %cohesion_d2=cohesion_d2/rows;
      %cohesion=cohesion_d1/cohesion_d2;
      allignment_d1=0;
      allignment_d2=0;
       for al=1:rows-1 % we are calculating the allignment of current dragon to other dragons
      allignment_d1=allignment_d1+Cosinesimilarity(dragon_grouped_at(al,:),dragon_grouped_at(rows,:));
      allignment_d2=allignment_d2+EuclideanSimilarity(dragon_grouped_at(al,:),dragon_grouped_at(rows,:));
      
      end
      
      allignment_d1=allignment_d1/rows;
      allignment_d2=allignment_d2/rows;
      allignment=allignment_d1/allignment_d2;
      
      seperation=0;
      if dragon_home_index==1
         seperation_cluster=2;
      else
          seperation_cluster=1;
      end
   
      seperation_best=dragon_global_best(seperation_cluster,:);
      seperation=0;
      seperation_d1=0;
      seperation_d2=0;
      
      for sep=1:size(dragon_grouped_at,1)
      seperation_d1=seperation_d1+Cosinesimilarity(dragon_grouped_at(sep,:),seperation_best);
      seperation_d2=seperation_d2+EuclideanSimilarity(dragon_grouped_at(sep,:),seperation_best);
      
      end
      [rows,cols]=size(dragon_grouped_at);
      seperation_d1=seperation_d1/rows;
      seperation_d2=seperation_d2/rows;
      seperation=seperation_d1/seperation_d2;
      
      [fitness_value,acceptance_probablity]=dragon_fitness(allignment,cohesion,seperation);
      
      acceptance_rejection_probablity(lv)=acceptance_probablity; % probablity of being accepted 
      acceptance_rejection_order(lv)=fitness_value; % boolean , either true or false
  
   

      
      
     
      
      
      
      
   end
   allevaluations(dr,1:levy_flights)=acceptance_rejection_probablity(1,:);
   f1=numel(find(acceptance_rejection_order==1));
   f0=numel(find(acceptance_rejection_order==0));
   
   selected=0;
   resultt='Rejected';
%    if f1>f0
%        selected=1;
%        sind=sind+1;
%        selected_indexes(sind,1)=dr;
%        resultt='Accepted';
%    end
%    
 mean_probablity=mean(acceptance_rejection_probablity);
   if mean_probablity>.30
        selected=1;
        sind=sind+1;
        selected_indexes(sind,1)=dr;
        resultt='Accepted';
       
   end
   figure(11);
   bar(acceptance_rejection_probablity);
   xlabel('Number of Levy Flights');
   ylabel('Acceptance Probablity');
  
   title(['Brisk index',num2str(dr),'Mean Problity:',num2str(mean_probablity),'Result',resultt])
   %%%%%%%
end
    catch
        errorfound=1;
    end
    if errorfound==0
        imgc=imgc+1;
image_index{imgc}=pathname;
for kkk=1:numel(selected_indexes)
Allfeatures(imgc,kkk)=dragon_points(selected_indexes(kkk));
end
disp(['Image:',num2str(img)])
    end
end
save allfeatures_updated_Test_data_5classes_with_lv