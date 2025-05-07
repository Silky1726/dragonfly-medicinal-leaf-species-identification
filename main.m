function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 13-Jan-2023 14:09:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Timg filename
set(handles.dip,'string','Please wait...');pause(0.25);
[filename,pathname]= uigetfile({'*.png'},'Select an Image');
fullpath=fullfile(pathname, filename);
Timg=imread(fullpath);
axes(handles.axes1);
imshow(Timg);title('Test Image','color','k','FontSize',12);
set(handles.dip,'string','Test Image Loaded !!!');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Timg Eimg
set(handles.dip,'string','Please wait...');pause(0.25);
Eimg=imadjust(Timg,stretchlim(Timg));
axes(handles.axes2);cla
imshow(Eimg)
title('Pre-processed Image','color','k','FontSize',12);
set(handles.dip,'string','Pre-processing Done!!!');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
clear
close all
main;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
clear
close all


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg Simg filename
set(handles.dip,'string','Please wait...');pause(0.25);
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
Bimg=im2bw(rgb2gray(PartSimg),0.5);
Bimg=imfill(Bimg,'holes');
se=strel('disk',50);
Bimg=imclose(Bimg,se);
BGpos=find(posimg==max(posimg));
axes(handles.axes3);cla
imshow(Simg{BGpos});title('BG Image','color','k','FontSize',12);
sts=regionprops(Bimg,'Eccentricity','Area');
eccn=[sts.Eccentricity];
ios=find(eccn);
boundaries=bwboundaries(Bimg);
Marea=0;
for k=1:length(ios)
if sts(k).Area>Marea
    Marea=sts(k).Area;
    krec=k;
end  
end
b = boundaries{krec};
[rows, columns, pln] = size(Eimg);
Segreg = poly2mask(b(:,2),b(:,1), rows, columns);
Redimg=Eimg(:,:,1);
Greenimg=Eimg(:,:,2);
Blueimg=Eimg(:,:,3);
FinalSimgR = double(Redimg).*double(imresize(Segreg,[size(Redimg,1) size(Redimg,2)]));
FinalSimgG = double(Greenimg).*double(imresize(Segreg,[size(Greenimg,1) size(Greenimg,2)]));
FinalSimgB = double(Blueimg).*double(imresize(Segreg,[size(Blueimg,1) size(Blueimg,2)]));
FinalSimg = uint8(cat(3,FinalSimgR,FinalSimgG,FinalSimgB));
axes(handles.axes4);cla
imshow(FinalSimg);title('Segmentated Image (K-means)','color','k','FontSize',12);

    image_info1=dir('GT Data/*.jpg');
    pathname1=filename;
    
    fullpath1=image_info1.folder;
    flpath1=[fullpath1,'/',pathname1];


TotalRegion=imread(flpath1);
TotalRegion=imresize(TotalRegion,[size(PartSimg,1) size(PartSimg,2)]);
if size(TotalRegion,3)>1
   A=im2bw(rgb2gray(TotalRegion));
else
   A=im2bw(TotalRegion); 
end
if size(FinalSimg,3)>1
   B=im2bw(rgb2gray(FinalSimg));
else
   B=im2bw(FinalSimg); 
end
[Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = SegmentationScores(A, B);
res=sprintf('1. Accuracy :%4f %%\n2. Sensitivity :%4f\n3. F-measure :%4f\n4. Precision :%.2f\n5. MCC :%.2f\n6. Dice :%.2f\n7. Jaccard :%.2f\n8. Specitivity :%.2f\n',Accuracy*100, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity);
msgbox(res,'Results (K-means)')
set(handles.dip,'string','Segmentation Done (K-means)!!!');

% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Timg filename
set(handles.dip,'string','Please wait...');pause(0.25);
[filename,pathname]= uigetfile({'*.jpg'},'Select an Image');
fullpath=fullfile(pathname, filename);
Timg=imread(fullpath);
axes(handles.axes1);
imshow(Timg);title('Test Image','color','k','FontSize',12);
set(handles.dip,'string','Test Image Loaded !!!');


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Timg Eimg
set(handles.dip,'string','Please wait...');pause(0.25);
Eimg=imadjust(Timg,stretchlim(Timg));
axes(handles.axes2);cla
imshow(Eimg)
title('Pre-processed Image','color','k','FontSize',12);
set(handles.dip,'string','Pre-processing Done!!!');

% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg Simg FinalSimg filename
set(handles.dip,'string','Please wait...');pause(0.25);
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
axes(handles.axes5);cla
imshow(Bimg);title('Mask Image','color','k','FontSize',12);
sts=regionprops(Bimg,'Eccentricity','Area');
eccn=[sts.Eccentricity];
ios=find(eccn);
axes(handles.axes6);cla;
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
axes(handles.axes7);cla;
imshow(FinalSimg,[]);
title('Segmented Image','color','k','FontSize',12);  
set(handles.dip,'string','Segmentation Done (FFA) !!!');


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg Simg FinalSimg filename
set(handles.dip,'string','Please wait...');pause(0.25);
posimg=[mean2(Simg{1});mean2(Simg{2})];
pos=find(posimg==min(posimg));
PartSimg=Simg{pos};
PSOThresh=[];
bit_value=PartSimg;
[r,c,p]=size(bit_value);
r=round(r/(32));
c=round(c/(32));
options=optimoptions('particleswarm','SwarmSize',50);
for si=1:r
for s=1:c
Fs=bit_value(si,s);
Ft=mean2(bit_value);
FitnessFunction=@(e)fitnessfunctionPSO(e,Fs,Ft);
numberOfVariables=1;
[~, fval]=particleswarm(FitnessFunction,numberOfVariables,[],[],options);
pause(0.025);
PSOThresh=abs(fval);
end
end
Bimg=im2bw(rgb2gray(PartSimg),PSOThresh);
Bimg=imfill(Bimg,'holes');
se=strel('disk',25);
Bimg=imclose(Bimg,se);
axes(handles.axes5);cla
imshow(Bimg);title('Mask Image','color','k','FontSize',12);
sts=regionprops(Bimg,'Eccentricity','Area');
eccn=[sts.Eccentricity];
ios=find(eccn);
axes(handles.axes6);cla;
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
[rows, columns, pln] = size(Bimg);
Segreg = poly2mask(b(:,2),b(:,1), rows, columns);
Redimg=Eimg(:,:,1);
Greenimg=Eimg(:,:,2);
Blueimg=Eimg(:,:,3);
FinalSimgR = double(Redimg).*double(imresize(Bimg,[size(Redimg,1) size(Redimg,2)]));
FinalSimgG = double(Greenimg).*double(imresize(Bimg,[size(Greenimg,1) size(Greenimg,2)]));
FinalSimgB = double(Blueimg).*double(imresize(Bimg,[size(Blueimg,1) size(Blueimg,2)]));
FinalSimg = uint8(cat(3,FinalSimgR,FinalSimgG,FinalSimgB));
axes(handles.axes7);cla;
imshow(FinalSimg,[]);
title('Segmented Image','color','k','FontSize',12);  
TotalRegion=imread(strcat('GT Data\',filename));
TotalRegion=imresize(TotalRegion,[size(PartSimg,1) size(PartSimg,2)]);
if size(TotalRegion,3)>1
   A=im2bw(rgb2gray(TotalRegion));
else
   A=im2bw(TotalRegion); 
end
if size(FinalSimg,3)>1
   B=im2bw(rgb2gray(FinalSimg));
else
   B=im2bw(FinalSimg); 
end
[Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = SegmentationScores(A, B);
res=sprintf('1. Accuracy :%4f %%\n2. Sensitivity :%4f\n3. F-measure :%4f\n4. Precision :%.2f\n5. MCC :%.2f\n6. Dice :%.2f\n7. Jaccard :%.2f\n8. Specitivity :%.2f\n',Accuracy*100, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity);
msgbox(res,'Results (K-means+PSO)')

set(handles.dip,'string','Segmentation Done (PSO) !!!');


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg Simg FinalSimg filename
set(handles.dip,'string','Please wait...');pause(0.25);
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
axes(handles.axes5);cla
imshow(Bimg);title('Mask Image','color','k','FontSize',12);
sts=regionprops(Bimg,'Eccentricity','Area');
eccn=[sts.Eccentricity];
ios=find(eccn);
axes(handles.axes6);cla;
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
axes(handles.axes7);cla;
imshow(FinalSimg,[]);
title('Segmented Image','color','k','FontSize',12);  

TotalRegion=imread(strcat('GT Data\',filename));
TotalRegion=imresize(TotalRegion,[size(PartSimg,1) size(PartSimg,2)]);
if size(TotalRegion,3)>1
   A=im2bw(rgb2gray(TotalRegion));
else
   A=im2bw(TotalRegion); 
end
if size(FinalSimg,3)>1
   B=im2bw(rgb2gray(FinalSimg));
else
   B=im2bw(FinalSimg); 
end
[Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = SegmentationScores(A, B);
res=sprintf('1. Accuracy :%4f %%\n2. Sensitivity :%4f\n3. F-measure :%4f\n4. Precision :%.2f\n5. MCC :%.2f\n6. Dice :%.2f\n7. Jaccard :%.2f\n8. Specitivity :%.2f\n',Accuracy*100, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity);
msgbox(res,'Results (K-means+FFA)')

set(handles.dip,'string','Segmentation Done (FFA) !!!');


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FinalSimg
set(handles.dip,'string','Please wait...');pause(0.25);
tic;
if size(FinalSimg,3)>1
GRegimg=rgb2gray(FinalSimg);
else
GRegimg=FinalSimg;
end
Keypoints=SIFT(GRegimg,FinalSimg);
axes(handles.axes8);imshow(Keypoints,[]); 
title('SIFT Points','color','k','FontSize',12);
save SIFT_Feature Keypoints
EtimeSIFT=toc;
ERSIFT=Param(Keypoints,FinalSimg);
save MSESIFT ERSIFT
save EtimeSIFT EtimeSIFT
%res=sprintf('1. Error :%4f %%\n2. Time :%4f Sec\n',ERSIFT/100, EtimeSIFT);
%msgbox(res,'Feature Results (SIFT)')
set(handles.dip,'string','Feature Extraction Done (SIFT)');

% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg FinalSimg
set(handles.dip,'string','Please wait...');pause(0.25);
tic;
if size(FinalSimg,3)>1
GRegimg=rgb2gray(FinalSimg);
else
GRegimg=FinalSimg;
end
Fpoints = detectSURFFeatures(GRegimg);
[features, validPoints] = extractFeatures(GRegimg,Fpoints);
axes(handles.axes8);imshow(uint8(Eimg)); hold on;
plot(Fpoints.Location(:,1),Fpoints.Location(:,2),'ro','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','y');
title('SURF Points','color','k','FontSize',12);
save SURF_Feature features validPoints
EtimeSURF=toc;
ERSURF=Param(features,Eimg);
save MSESURF ERSURF
save EtimeSURF EtimeSURF
%res=sprintf('1. Error :%4f %%\n2. Time :%4f Sec\n',ERSURF, EtimeSURF);
%msgbox(res,'Feature Results (SURF)')
set(handles.dip,'string','Feature Extraction Done (SURF)');

% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Eimg FinalSimg
set(handles.dip,'string','Please wait...');pause(0.25);
tic;
if size(FinalSimg,3)>1
GRegimg=rgb2gray(FinalSimg);
else
GRegimg=FinalSimg;
end
Fpoints = detectBRISKFeatures(GRegimg);
[features, validPoints] = extractFeatures(GRegimg,Fpoints);
axes(handles.axes8);imshow(uint8(Eimg)); hold on;
plot(Fpoints.Location(:,1),Fpoints.Location(:,2),'ro','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','y');
title('BRISK Points','color','k','FontSize',12);
save BRISK_Feature features validPoints
EtimeBRISK=toc;
ERBRISK=Param(features.Features,Eimg);
save('Briskpoints','Fpoints');
save MSEBRISK ERBRISK
save EtimeBRISK EtimeBRISK
%res=sprintf('1. Error :%4f %%\n2. Time :%4f Sec\n',ERBRISK, EtimeBRISK);
%msgbox(res,'Feature Results (BRISK)')
set(handles.dip,'string','Feature Extraction Done (BRISK)');


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load Briskpoints
allevaluations=[];
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
   acceptance_rejection_probablity=zeros(1,levy_flights); % probablity of being accepted 
   acceptance_rejection_order=zeros(1,levy_flights); % boolean , either true or false
   for lv=1:levy_flights
      dragon_grouping_indexes=round(rand(1,dragons_to_group)*dragons_to_group);
      dragon_grouping_indexes=findandreplace(dragon_grouping_indexes,0,dragon_population(1));
      dragon_grouped_at=dragon_home_all_at(dragon_grouping_indexes,:);
      dragon_grouped_at(dragons_to_group+1,1)=current_dragon_at(1,1);
      dragon_grouped_at(dragons_to_group+1,1)=current_dragon_at(1,2);
      global_best=dragon_global_best(dragon_home_index,:);
      cohesion_d1=0;% cosine similarity 
      cohesion_d2=0; % euclidean distance 
      
      for coh=1:size(dragon_grouped_at,1)
      cohesion_d1=cohesion_d1+Cosinesimilarity(dragon_grouped_at(coh,:),global_best);
      cohesion_d2=cohesion_d2+EuclideanSimilarity(dragon_grouped_at(coh,:),global_best);
      
      end
      [rows,cols]=size(dragon_grouped_at);
      cohesion_d1=cohesion_d1/rows;
      cohesion_d2=cohesion_d2/rows;
      cohesion=cohesion_d1/cohesion_d2;
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
      
      acceptance_rejection_probablity(lv)=acceptance_probablity % probablity of being accepted 
      acceptance_rejection_order(lv)=fitness_value % boolean , either true or false
  
   

      
      
     
      
      
      
      
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
   if mean_probablity>.20
        selected=1;
        sind=sind+1;
        selected_indexes(sind,1)=dr;
        resultt='Accepted';
       
   end
   figure(11);
   bar(acceptance_rejection_probablity);
   xlabel('Number of Levy Flights');
   ylabel('Acceptance Probablity');
  
   title(['Brisk index',num2str(dr),'Mean Problity:',num2str(mean_probablity),'Result',resultt]);
   
   %%%%%%%
end
save('indexesbydragon','selected_indexes');
global Eimg FinalSimg
%set(handles.dip,'string','Please wait...');pause(0.25);
tic;
load indexesbydragon
if size(FinalSimg,3)>1
GRegimg=rgb2gray(FinalSimg);
else
GRegimg=FinalSimg;
end
Fpoints1=[];


[features, validPoints] = extractFeatures(GRegimg,Fpoints);
figure(23);imshow(uint8(Eimg)); hold on;
plot(Fpoints.Location(:,1),Fpoints.Location(:,2),'ro','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','y');
title('BRISK Points','color','k','FontSize',12);
save BRISK_Feature features validPoints
EtimeBRISK=toc;
ERBRISK=Param(features.Features,Eimg);
save('Briskpoints','Fpoints');
save MSEBRISK ERBRISK
save EtimeBRISK EtimeBRISK
%res=sprintf('1. Error :%4f %%\n2. Time :%4f Sec\n',ERBRISK, EtimeBRISK);
%msgbox(res,'Feature Results (BRISK)')
set(handles.dip,'string','Feature Extraction Done (BRISK)');

Scale=Fpoints.Scale(selected_indexes);
Orientation=Fpoints.Orientation(selected_indexes);
Location=Fpoints.Location(selected_indexes);
Metric=Fpoints.Metric(selected_indexes);
% % Fpoints.Scale=Scale;
% % Fpoints.Orientation=Orientation;
% % Fpoints.Location=Location;
% % Fpoints.Metric=Metric;

% [features, validPoints] = extractFeatures(GRegimg,Fpoints);
figure(23); hold on;
plot(Fpoints.Location(selected_indexes,1),Fpoints.Location(selected_indexes,2),'bo','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','y');
title('BRISK Points','color','k','FontSize',12);
save BRISK_Feature features validPoints
save dragonfeatures Metric;
EtimeBRISK=toc;
ERBRISK=Param(features.Features,Eimg);
save('Briskpoints','Fpoints');
save MSEBRISK ERBRISK
save EtimeBRISK EtimeBRISK
%res=sprintf('1. Error :%4f %%\n2. Time :%4f Sec\n',ERBRISK, EtimeBRISK);
%msgbox(res,'Feature Results (BRISK)')
set(handles.dip,'string','Feature Extraction Done (BRISK)');


% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load Briskpoints

PSO_scale=double(Fpoints.Scale);
PSO_points=double(Fpoints.Metric);
bit_value=PSO_points;
[r,c]=size(PSO_points);
options=optimoptions('particleswarm','SwarmSize',50);
for si=1:r
Fs=bit_value(si);
Ft=mean2(bit_value(1:si));
FitnessFunction=@(e)fitnessfunctionPSO(e,Fs,Ft);
numberOfVariables=1;
[~, fval]=particleswarm(FitnessFunction,numberOfVariables,[],[],options);
pause(0.025);
PSOThresh=abs(fval);
PSO_points(si)=PSOThresh;

end
test_feature=PSO_points
save PSOFeatures_test test_feature
msgbox('PSO Feature List Updated');


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load dragonfeatures.mat
load neuralvalue.mat
%load decisionTreeModel.mat
indexes=310;
testset= double(Metric(1:indexes));
score = sim(net, testset);
[~, result] = max(score, [], 1);
%msgbox(['Class Result Species',num2str(result(1))]);
%score = sim(net, testset);
%[~, result] = max(score, [], 1);
%result = predict(decisionTreeModel, testset);
%msgbox(['Class Result Species',num2str(result(1))]);
if result(1) == 1
        speciesName = 'Apple';
elseif result(1) == 2
        speciesName = 'Blueberry';
else
        speciesName = 'Unknown species';
end
    
    % Display the result in a message box
    msgbox(['Class Result: Species ', speciesName]);
