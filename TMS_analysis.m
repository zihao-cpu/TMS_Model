
raw={};
raw1={};

for i=1:18
    
    [quanpin,firstword] = FunPinyin(raw{i,1});
    name=[];
    for j=1:length(quanpin)
        name=[name,char(quanpin(j))];
        
    end

    name1{i,1}=lower(name);
    
    
end

for i=1:14
    
    [quanpin,firstword] = FunPinyin(raw1{i,1});
    name=[];
    for j=1:length(quanpin)
        name=[name,char(quanpin(j))];
        
    end

    name2{i,1}=lower(name);
    
end
%%
load('subject_excel.mat');
files=dir('I:\D1D2-33个0318to梓豪\D2');
files([1,2])=[];

% for i=1:length(files)
%     movefile([files(i).folder,'\',files(i).name],[files(i).folder,'\',lower(files(i).name(14:end))]);
% end




suball=cellfun(@(name) lower([name{2},'_',name{3},'_',name{4}]),...
    cellfun(@(names) strsplit(names,'_'),...
    {files.name},'UniformOutput',false),...
    'UniformOutput',false);

suball=cellfun(@(name) lower([name{1}]),...
    cellfun(@(names) strsplit(names,'_'),...
    {files.name},'UniformOutput',false),...
    'UniformOutput',false);

path1='I:\TMS\Nii\d2\fun_preprocess\normalise'
files=dir(path1);
files(1:2)=[];
d1all=cellfun(@(name) lower([name{1}]),...
    cellfun(@(names) strsplit(names,'_'),...
    {files.name},'UniformOutput',false),...
    'UniformOutput',false);

[a1,b1]=ismember(d1all,suball);
suball(find(a1==0));






%%

name1d1=name1;
name1d2=name1;
pathd1='I:\TMS\Nii\d1\fun_preprocess\normalise';

files=dir(pathd1);
files(1:2)=[];
d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name1d1(:,1),d1sub);

pathd2='I:\TMS\Nii\d2\fun_preprocess\normalise';
files=dir(pathd2);
files(1:2)=[];
d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];

nSub=length(b1);
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nDimTimePoints,1);
    SubjectRegressors(i:nSub:nDimTimePoints,i) = 1;
end
Regressors = [Regressors,SubjectRegressors];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;

    

    
%% 



load('I:\TMS\Nii\new_results\sub.mat');
name1d1=name1;
name1d2=name1;
name1d1{6}='chenkaiju';
name1d2{6}='chenkaiju';
name1d1{14}='fang_mei';
name1d2{11}='jiachuanhua';

sex=raw(a1,2);
sex=cellfun(@(temp) strcmp(temp,'女'),sex);
sex=double(sex);
sex(sex==0)=2;

age=cell2mat(raw(a1,3));
edu=cell2mat(raw(a1,4));

load('I:\TMS\Nii\new_results\behave.mat')
[a1,b1]=ismember(raw1(:,1),behave_D1(:,1));

[a2,b2]=ismember(raw1(:,1),behave_D2(:,1));


behave_D1=cell2mat(behave_D1(b1,6:15));
behave_D2=cell2mat(behave_D2(b2,6:15));
names={};

for i=1:10
  [h,p(i),c,d]=ttest(behave_D2(:,i),behave_D1(:,i));
  t(i)=d.tstat;
  close all;
  plot_specificity_box(behave_D1(:,i),behave_D2(:,i));
 print(gcf,['I:\TMS\Nii\behave_movie\',names{i},'.png'],'-dpng','-r600');
end


for i=1:10
  [h,p(i),c,d]=ttest(behave_D2(:,i),behave_D1(:,i));
  t(i)=d.tstat;
  close all;
  plot_specificity_box(behave_D1(:,i),behave_D2(:,i));
 print(gcf,['I:\TMS\Nii\behave\',names{i},'.png'],'-dpng','-r600');
end





dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\ReHo';
pathd2='I:\TMS\Nii\d2\fun_preprocess\ReHo';

files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
% 
%     dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
    
        dataall=[dataall;data2(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
% 
%     dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
    
          dataall=[dataall;data1(graymask~=0)'];
      
end



%%
pathd1='I:\TMS\Nii\d1\fun_preprocess\alff';
pathd2='I:\TMS\Nii\d2\fun_preprocess\alff';

files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];

dataall=[];
MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    
%     Vsourcesmoothed./mean(Vsourcesmoothed(graymask~=0));
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
%      Vsourcesmoothed./mean(Vsourcesmoothed(graymask~=0));
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end


%%
dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\FCD\Correlation_0.6\longRange_FCD';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FCD\Correlation_0.6\longRange_FCD';
name='longRange_FCD';
files1d1=dir(pathd1);
files1d1(1:2)=[];

files1d2=dir(pathd2);
files1d2(1:2)=[];


files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];



for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
%     data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
%     data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end




%%
dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\FOCA\mean';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FOCA\mean';
name='FOCA'
files1d1=dir(pathd1);
files1d1(1:2)=[];

files1d2=dir(pathd2);
files1d2(1:2)=[];


files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];



for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
%     data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
%     data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end
%%

dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\FC\140\z*';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FC\140\z*';
name='140'
MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));
files1d1=dir(pathd1);

files1d2=dir(pathd2);


d1sub=cellfun(@(name) name{1}(4:end),cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);

d2sub=cellfun(@(name) name{1}(4:end),cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


pathd1='I:\TMS\Nii\d1\fun_preprocess\FC\140';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FC\140';
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
%     data2=data2./mean(data2(graymask~=0));  
%     smooth_factor=6;
%     Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;data2(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
%     data1=data1./mean(data1(graymask~=0));  
%     smooth_factor=6;
%     Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;data1(graymask~=0)'];
      
end
%%
dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\Spreading';
pathd2='I:\TMS\Nii\d2\fun_preprocess\Spreading';
% name='FOCA'

files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


num=0;
for i=1:length(b2)
    num=num+1;
    data2=load([pathd2,'\',files1d2(b2(i)).name]);

   S_mci_nc(num,:,:)=squeeze(data2.T(:,:,7));


      
end

for i=1:length(b1)
       num=num+1;
    data1=load([pathd1,'\',files1d1(b1(i)).name]);

   S_mci_nc(num,:,:)=squeeze(data1.T(:,:,7));

      
end
pathd1='I:\TMS\Nii\d1\fun_preprocess\realignParameters';
pathd2='I:\TMS\Nii\d2\fun_preprocess\realignParameters';
files1d1=dir(pathd1);
files1d1(1:2)=[];

files1d2=dir(pathd2);
files1d2(1:2)=[];
num=0;
for i=1:length(b2)
    num=num+1;
    head=dir([pathd2,'\',files1d2(b2(i)).name]);
    rp1=importdata([head(3).folder,'\',head(3).name]);
    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];
    
     meanFD(num,1)=mean(FD);

 
end

for i=1:length(b1)


    num=num+1;
    head=dir([pathd1,'\',files1d1(b1(i)).name]);
    rp1=importdata([head(3).folder,'\',head(3).name]);
    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];
    
     meanFD(num,1)=mean(FD);


      
end





S_mci_nc=S_mci_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==32)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end

new=S_mci_nc_temp;
new(isnan(new))=0;
T2map=[];
Regressors=[Regressors,meanFD];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;

for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
%     network=zscore(network')';
    for i=1:416
        DependentVariable=network(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors,Contrast,'T'); 
    end
        Ttemp=TF_ForContrast';
           T2map=[T2map;TF_ForContrast];
    df=size(Regressors,1)-...
        size(Regressors,2);
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
%     
%     if ~isempty(P)
%     Thresholded(find(PMap<=P))=1;
%     end
    

        P=0.01;
    Thresholded(find(PMap<=P))=1;

    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
T2map=[];
for net=1:18
    network=new(:,:,yeoindex17==net);

    network=mean(network,3);
    network=squeeze(network);
%     network=zscore(network')';
    for i=1:416
        DependentVariable=network(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors,Contrast,'T'); 
    end
        Ttemp=TF_ForContrast';
           T2map=[T2map;TF_ForContrast];
    df=size(Regressors,1)-...
        size(Regressors,2);
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
%     
%     if ~isempty(P)
%     Thresholded(find(PMap<=P))=1;
%     end
    

        P=0.01;
    Thresholded(find(PMap<=P))=1;

    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=data_mat_mean_all;
for i=1:18
for j=1:18
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors,Contrast,'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;

end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors,1)-...
    size(Regressors,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;


TF_ForContrast_brain_mean3=[];
for i=1:416
for j=1:416
DependentVariable=new(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors,Contrast,'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;

end
end
TF_ForContrast_brain_mean3(find(eye(416)))=0;
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors,1)-...
    size(Regressors,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(416,416);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;







new=S_mci_nc_temp(:,yeo_17labels_17_to_7,yeo_17labels_17_to_7);
new(isnan(new))=0;
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:8];
    mask=yeoindex;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=data_mat_mean_all;
for i=1:8
for j=1:8
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors,Contrast,'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors,1)-...
    size(Regressors,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(8,8);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;


%%
nSub=length(b1);
nDimTimePoints=2*nSub;
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nDimTimePoints,1);
    SubjectRegressors(i:nSub:nDimTimePoints,i) = 1;
end
Regressors = [Regressors,SubjectRegressors];

% Regressors = [Regressors,SubjectRegressors,[age;age],[sex;sex],[edu;edu]];

% Regressors = [Regressors,SubjectRegressors];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;


for i=1:size(dataall,2);
    DependentVariable=dataall(:,i);
  [~,r,SSE,SSR, ~, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,[Regressors],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
        
end
TF_ForContrast(isnan(TF_ForContrast))=0;
Df1=size(Regressors,1)-size(Regressors,2);
T=zeros(91,109,91);
T(graymask~=0)=TF_ForContrast;
vgraymask=spm_vol(MaskFile);


vgraymask.fname=['I:\TMS\Nii\new_results\',name,'_T.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);

T=TF_ForContrast';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;
Ttemp=T(mask~=0);
df=Df1;
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if isempty(P)
    P=0.05
    Thresholded(find(PMap<=P))=1;
else
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,1);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T.*Thresholdedmask;


T_pos=T_fdr;
T_pos(T_pos<=0)=0;
T=zeros(91,109,91);
T(graymask~=0)=T_pos;
vgraymask=spm_vol(MaskFile);
vgraymask.fname=['I:\TMS\Nii\new_results\',name,'_T_pos.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);


T_neg=T_fdr;
T_neg(T_neg>=0)=0;
T=zeros(91,109,91);
T(graymask~=0)=T_neg;
vgraymask=spm_vol(MaskFile);
vgraymask.fname=['I:\TMS\Nii\new_results\',name,'_T_neg.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);






% IsTwoTailed=0;
% StatsImgFile='I:\TMS\Nii\new_results\alff_T.nii';
% name='alff';
% VoxelPThreshold=0.05;
% ClusterPThreshold=0.05;
% [BrainVolume, VoxelSize, Header]=y_ReadRPI(StatsImgFile);
% Header_DF = w_ReadDF(Header);
% Flag = Header_DF.TestFlag;
% Df1 = Header_DF.Df;
% Df2 = Header_DF.Df2;
% [nDim1,nDim2,nDim3]=size(BrainVolume);
% Outpath=fileparts(StatsImgFile);
% OutputName=[Outpath,'\',name,'_GRF_Vp_',num2str(VoxelPThreshold),'_Cp_',num2str(ClusterPThreshold),'_.nii'];
% MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
%  [Data_Corrected, ClusterSize, Header,BrainVolumeRawStats]=y_GRF_Threshold_zzh(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,...
%     OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,[],[]);









%%


load('I:\TMS\Nii\new_results\sub.mat');
name2d1=name2;
name2d2=name2;

name2d1{14}='zhaoxiaocui';

sex=raw1(a1,2);
sex=cellfun(@(temp) strcmp(temp,'女'),sex);
sex=double(sex);
sex(sex==0)=2;

age=cell2mat(raw1(a1,3));
edu=cell2mat(raw1(a1,4));



pathd1='I:\TMS\Nii\d1\fun_preprocess\alff';
pathd2='I:\TMS\Nii\d2\fun_preprocess\alff';
name='alff';
files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name2d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

dataall=[];
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end
%%


pathd1='I:\TMS\Nii\d1\fun_preprocess\ReHo';
pathd2='I:\TMS\Nii\d2\fun_preprocess\ReHo';
name='ReHo';
files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name2d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

dataall=[];
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=3;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

%%
pathd1='I:\TMS\Nii\d1\fun_preprocess\FCD\Correlation_0.6\global_FCD';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FCD\Correlation_0.6\global_FCD';
name='global';
files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name2d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

dataall=[];
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end
%%
pathd1='I:\TMS\Nii\d1\fun_preprocess\FCD\Correlation_0.6\longRange_FCD';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FCD\Correlation_0.6\longRange_FCD';
name='longRange';
files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name2d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

dataall=[];
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end
%%
pathd1='I:\TMS\Nii\d1\fun_preprocess\FOCA\mean';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FOCA\mean';
name='FOCA';
files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
[a1,b1]=ismember(name2d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{2},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

dataall=[];
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
    data2=data2./mean(data2(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
    data1=data1./mean(data1(graymask~=0));  
    smooth_factor=6;
    Vsourcesmoothed=zeros(91,109,91);
    spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);

    dataall=[dataall;Vsourcesmoothed(graymask~=0)'];
      
end

%%
dataall=[];
pathd1='I:\TMS\Nii\d1\fun_preprocess\FC\140\z*';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FC\140\z*';
name='140'

MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));
files1d1=dir(pathd1);

files1d2=dir(pathd2);


d1sub=cellfun(@(name) name{1}(4:end),cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);

[a1,b1]=ismember(name2d1(:,1),d1sub);


d2sub=cellfun(@(name) name{1}(4:end),cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name2d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


pathd1='I:\TMS\Nii\d1\fun_preprocess\FC\140';
pathd2='I:\TMS\Nii\d2\fun_preprocess\FC\140';
for i=1:length(b2)
    data2=spm_read_vols(spm_vol([pathd2,'\',files1d2(b2(i)).name]));
%     data2=data2./mean(data2(graymask~=0));  
%     smooth_factor=6;
%     Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data2, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;data2(graymask~=0)'];
      
end

for i=1:length(b1)
    data1=spm_read_vols(spm_vol([pathd1,'\',files1d1(b1(i)).name]));
%     data1=data1./mean(data1(graymask~=0));  
%     smooth_factor=6;
%     Vsourcesmoothed=zeros(91,109,91);
%     spm_smooth(data1, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);
    dataall=[dataall;data1(graymask~=0)'];
      
end



%%
nSub=length(b1);
nSub=14;

nDimTimePoints=2*nSub;
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nDimTimePoints,1);
    SubjectRegressors(i:nSub:nDimTimePoints,i) = 1;
end
Regressors = [Regressors,SubjectRegressors];

% Regressors = [Regressors,SubjectRegressors,[age;age],[sex;sex],[edu;edu]];

% Regressors = [Regressors,SubjectRegressors];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;


for i=1:size(dataall,2);
    DependentVariable=dataall(:,i);
  [~,r,SSE,SSR, ~, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,[Regressors],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
        
end
TF_ForContrast(isnan(TF_ForContrast))=0;

Df1=size(Regressors,1)-size(Regressors,2);
T=zeros(91,109,91);
T(graymask~=0)=TF_ForContrast;
vgraymask=spm_vol(MaskFile);


vgraymask.fname=['I:\TMS\Nii\new_results_movie\',name,'_T.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);

T=TF_ForContrast';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;
Ttemp=T(mask~=0);
df=Df1;
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if isempty(P)
    P=0.05
    Thresholded(find(PMap<=P))=1;
else
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,1);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T.*Thresholdedmask;


T_pos=T_fdr;
T_pos(T_pos<=0)=0;
T=zeros(91,109,91);
T(graymask~=0)=T_pos;
vgraymask=spm_vol(MaskFile);
vgraymask.fname=['I:\TMS\Nii\new_results_movie\',name,'_T_pos.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);


T_neg=T_fdr;
T_neg(T_neg>=0)=0;
T=zeros(91,109,91);
T(graymask~=0)=T_neg;
vgraymask=spm_vol(MaskFile);
vgraymask.fname=['I:\TMS\Nii\new_results_movie\',name,'_T_neg.nii'];
vgraymask.dt=[16,0];
vgraymask.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',Df1,0,0,0,0);
spm_write_vol(vgraymask,T);






%%

dataD1=load_data_TMS_D1('ReHo');

dataD2=load_data_TMS_D2('ReHo');

dataD1_1=dataD1(size(dataD1,1)/2+1:end,:); %第一种基线

dataD1_2=dataD1(1:size(dataD1,1)/2,:);%第一种治疗


dataD2_1=dataD2(size(dataD2,1)/2+1:end,:);

dataD2_2=dataD2(1:size(dataD2,1)/2,:);
for name={'global_FCD','local_FCD','longRange_FCD','FOCA'}
    dataD1=load_data_TMS_D1(name);

dataD2=load_data_TMS_D2(name);

dataD1_1=dataD1(size(dataD1,1)/2+1:end,:); %第一种基线

dataD1_2=dataD1(1:size(dataD1,1)/2,:);%第一种治疗


dataD2_1=dataD2(size(dataD2,1)/2+1:end,:);

dataD2_2=dataD2(1:size(dataD2,1)/2,:);
save(['Anova_',char(name)],'dataD1_2','dataD1_1','dataD2_2','dataD2_1');
    
    
end






sub=repmat([1:size(dataD1,1)/2+size(dataD2,1)/2]',2,1);
c1=[ones((size(dataD1,1)/2+size(dataD2,1)/2),1);...
    ones((size(dataD1,1)/2+size(dataD2,1)/2),1)*2];
c2=repmat([ones(size(dataD1,1)/2,1);ones(size(dataD2,1)/2,1)*2],2,1);





MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));


c1=repmat([1;1;2;2],(size(dataD1,1)+size(dataD2,1))/4,1)

c2=repmat([1;2],31,1);
for i=1:size(dataD1,1)/2+size(dataD2,1)/2
    sub(2*i-1)=i;
     sub(2*i)=i;
    
end




parfor i=1:size(dataD1,2)
    temp=zeros(length(sub),1);
    temp(c1==1 & c2==1)=dataD1_1(:,i); %第一种基线
    
    temp(c1==2 & c2==1)=dataD2_1(:,i);%第2种基线
    
    temp(c1==1 & c2==2)=dataD1_2(:,i); %第1种治疗
    
    temp(c1==2 & c2==2)=dataD2_2(:,i);%第二种治疗
    
     [F,P,DF] = RM_2WayAnova([temp,c1,c2,sub]);
        Fc1(i)=F(1);
        Fc2(i)=F(2);
        Fc3(i)=F(3);
        Fc4(i)=F(4);
    
    
    
end









for i=1:size(dataD1,2)
   temp=[[ones(size(dataD1_1,1),1);ones(size(dataD2_1,1),1)*2],[[dataD1_1(:,i),dataD1_2(:,i)];[dataD2_1(:,i),dataD2_2(:,i)]]];
   model = mixed_anova(temp, 0); 
   Between_F(i)=model.between.F;
   Within_F(i)=model.within.F;
   Inter_F(i)=model.interaction.F; 
end


design=temp;
groups      = unique(design(:,1));
ngroups     = length(groups);
gind        = cell(size(groups));
group_n     = zeros(1,length(gind));
nlevels     = size(design,2) - 1;
% compute grand mean, cell means and marginal totals
data        = design(:,2:end);
for ii = 1:length(gind);
    gind{ii}                 = find(design(:,1) == groups(ii));
    group_data{ii}           = data(gind{ii},:);
    group_n(ii)              = length(find(design(:,1) == groups(ii)));
    cell_n(ii,:)             = repmat(group_n(ii),1,nlevels);
    cell_means(ii,1:nlevels) = mean(data(gind{ii},:));
    subj_subtotals(ii)       = sum(sum(data(gind{ii},:)));
    subj_ratios(ii)          = subj_subtotals(ii)^2/(nlevels*group_n(ii));
    group_subtot(ii,:)       = sum(group_data{ii});
    group_ratios(ii)         = sum(group_subtot(ii,:).^2)/group_n(ii);
end  
df_A=length(unique(temp(:,1)))-1;
df_SA=sum(group_n-1);
df_B = nlevels - 1 ;
df_SBA = df_B * df_SA;
df_AB = df_A * df_B;


MaskFile='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
vgraymask=spm_vol(MaskFile);
graymask=spm_read_vols(spm_vol(vgraymask));

Fc3=zeros(91,109,91);
Fc3(graymask~=0)=Between_F;
vimage=spm_vol(MaskFile);
vimage.descrip=sprintf('DPABI{F_[%.1f,%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',df_A,df_SA,0,0,0,0);
vimage.dt=[16,0];
vimage.fname='I:\TMS\Nii\Anova\reho_BF.nii';
spm_write_vol(vimage,Fc3);


Fc3=zeros(91,109,91);
Fc3(graymask~=0)=Within_F;
vimage=spm_vol(MaskFile);
vimage.descrip=sprintf('DPABI{F_[%.1f,%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',df_B,df_SBA,0,0,0,0);
vimage.dt=[16,0];
vimage.fname='I:\TMS\Nii\Anova\reho_WF.nii';
spm_write_vol(vimage,Fc3);

Fc3=zeros(91,109,91);
Fc3(graymask~=0)=Inter_F;
vimage=spm_vol(MaskFile);
vimage.descrip=sprintf('DPABI{F_[%.1f,%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',df_AB,df_SBA,0,0,0,0);
vimage.dt=[16,0];
vimage.fname='I:\TMS\Nii\Anova\reho_IF.nii';
spm_write_vol(vimage,Fc3);



files=dir('I:\TMS\Nii\Anova\*.nii');
for i=1:length(files)
    
    [BrainVolume, VoxelSize, Header]=y_ReadRPI([files(i).folder,'\',files(i).name]);
    Header_DF = w_ReadDF(Header);
    
    P = 1-fcdf(BrainVolume,  Header_DF.Df,  Header_DF.Df2);
    P(graymask==0)=0;
    
    BrainVolume(P>=0.01)=0;
    
    vimage=spm_vol(MaskFile);
    vimage.descrip=sprintf('DPABI{F_[%.1f,%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',df_AB,df_SBA,0,0,0,0);
    vimage.dt=[16,0];
    vimage.fname=['I:\TMS\Nii\Anova\',files(i).name(1:end-4),'lt_0.01.nii'];
    spm_write_vol(vimage,BrainVolume);

    
end


%%
clear all
load('I:\TMS\Nii\new_results\sub.mat');
name1d1=name1;
name1d2=name1;
name1d1{6}='chenkaiju';
name1d2{6}='chenkaiju';
name1d1{14}='fang_mei';
name1d2{11}='jiachuanhua';


path1='I:\TMS\Nii\d1\ASL_new\uninter';
uninter=dir(path1);
uninter(1:2)=[];
path2='I:\TMS\Nii\d1\ASL_new\inter';
inter=dir(path2);
inter(1:2)=[];

uninter={uninter.name};
inter={inter.name};
alldatapath='I:\TMS\Nii\d1\ASL';
alldata=dir(alldatapath);
alldata(1:2)=[];
alldata={alldata.name};

cortex=spm_read_vols(spm_vol('E:\WorkSpace\MRI_MASK\r2schaefer400MNI.nii'));
cortex(isnan(cortex))=0;
subcortex=spm_read_vols(spm_vol('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii'));


[a1,b1]=ismember(uninter,alldata);
[a2,b2]=ismember(inter,alldata);
path1='I:\TMS\Nii\d1\ASL_new\uninter';
cd(path1);
subfile=dir(pwd);
M_youngcortex=[];M_youngsubcortex=[];MRIdata_in=[];
for i=3:length(subfile)
    
    asl=spm_read_vols(spm_vol([subfile(i).name,'\FUNC\wmeanCBF_0_sASLflt_rASL.nii']));
    asl(asl<=0)=0;
%     asl=(asl./mean(asl))-1;
%     data(:,i-2)=asl(graymask~=0);


   image=asl;
    for roi=1:400
    M_youngcortex(roi,i-2)=mean(image(cortex==roi));
    end
    
        for roi=1:16
        M_youngsubcortex(roi,i-2)=mean(image(subcortex==roi));
    end
    
end
MRIdata_un=[ M_youngcortex;M_youngsubcortex];

path1='I:\TMS\Nii\d1\ASL_new\inter';
cd(path1);
subfile=dir(pwd);
 M_youngcortex=[];M_youngsubcortex=[];MRIdata_in=[];
for i=3:length(subfile)
    
    asl=spm_read_vols(spm_vol([subfile(i).name,'\FUNC\wmeanCBF_0_sASLflt_rASL.nii']));
     asl(asl<=0)=0;
%     asl=(asl./mean(asl))-1;
%     data(:,i-2)=asl(graymask~=0);


   image=asl;
    for roi=1:400
    M_youngcortex(roi,i-2)=mean(image(cortex==roi));
    end
    
        for roi=1:16
        M_youngsubcortex(roi,i-2)=mean(image(subcortex==roi));
    end
    
end
MRIdata_in=[M_youngcortex;M_youngsubcortex];

ASL=zeros(32,416);
ASL(b1,:)=MRIdata_un';
ASL(b2,:)=MRIdata_in';

Boldpath='I:\TMS\Nii\d1\fun_preprocess\BOLD_smooth';

Boldpath='I:\TMS\Nii\d1\fun_preprocess\BOLD';

cd(Boldpath)
boldfile=dir('*.mat');
for i=[1:32]
    load(boldfile(i).name);
    BOLD=theROITimeCoursesTotal;
    clear beta1
    for time=1:size(BOLD,1)
        
%         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(time,:))',[zscore(ASL(i,:)'),ones(416,1)],[1,0],'T');
        beta1(time)=b(1);
           
    end
    
     for roi=1:416
        
%         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(:,roi)),[zscore(beta1'),ones(size(BOLD,1),1)],[1,0],'T');
        beta2(roi,i)=b(1);
           
    end
    
    
    
end
TMS_d1=beta2;


path1='I:\TMS\Nii\d2\ASL_new';
uninter=dir(path1);
uninter(1:2)=[];
uninter={uninter.name};
  
uninter{11}='jiachuanhua_20240108_d2';

alldatapath='I:\TMS\Nii\d2\fun_preprocess\BOLD_smooth';

alldatapath='I:\TMS\Nii\d2\fun_preprocess\BOLD';
alldata=dir(alldatapath);
alldata(1:2)=[];
alldata=cellfun(@(name) name(1:end-21),{alldata.name},'UniformOutput',false);


[a1,b1]=ismember(uninter,alldata);

path1='I:\TMS\Nii\d2\ASL_NEW';
cd(path1);
subfile=dir(pwd);
M_youngcortex=[];M_youngsubcortex=[];MRIdata_in=[];
for i=3:length(subfile)
    
    asl=spm_read_vols(spm_vol([subfile(i).name,'\FUNC\wmeanCBF_0_sASLflt_rASL.nii']));
    asl(asl<=0)=0;
%     asl=(asl./mean(asl))-1;
%     data(:,i-2)=asl(graymask~=0);


   image=asl;
    for roi=1:400
    M_youngcortex(roi,i-2)=mean(image(cortex==roi));
    end
    
        for roi=1:16
        M_youngsubcortex(roi,i-2)=mean(image(subcortex==roi));
    end
    
end
MRIdata_un=[ M_youngcortex;M_youngsubcortex];



ASL=zeros(31,416);
ASL(b1,:)=MRIdata_un';


Boldpath='I:\TMS\Nii\d2\fun_preprocess\BOLD_smooth';

Boldpath='I:\TMS\Nii\d2\fun_preprocess\BOLD';

cd(Boldpath)
boldfile=dir('*.mat');
clear beta2;
for i=[1:31]
    load(boldfile(i).name);
    BOLD=theROITimeCoursesTotal;
    clear beta1
    for time=1:size(BOLD,1)
        
%         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(time,:))',[zscore(ASL(i,:)'),ones(416,1)],[1,0],'T');
        beta1(time)=b(1);
           
    end
    
     for roi=1:416
        
%         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(:,roi)),[zscore(beta1'),ones(size(BOLD,1),1)],[1,0],'T');
        beta2(roi,i)=b(1);
           
    end
    
    
    
end
TMS_d2=beta2;



pathd1='I:\TMS\Nii\d1\fun_preprocess\BOLD';
pathd2='I:\TMS\Nii\d2\fun_preprocess\BOLD';


pathd1='I:\TMS\Nii\d1\fun_preprocess\BOLD_smooth';
pathd2='I:\TMS\Nii\d2\fun_preprocess\BOLD_smooth';

files1d1=dir(pathd1);
files1d1(1:2)=[];
d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';
[a1,b1]=ismember(name1d1(:,1),d1sub);


files1d2=dir(pathd2);
files1d2(1:2)=[];
d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{files1d2.name},'UniformOutput',false),...
    'UniformOutput',false);
[a2,b2]=ismember(name1d2(:,1),d2sub);

b1(b1==0)=[];
b2(b2==0)=[];


dataall=[TMS_d2(:,b2),TMS_d1(:,b1)];

nSub=length(b1);
nDimTimePoints=2*nSub;
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nDimTimePoints,1);
    SubjectRegressors(i:nSub:nDimTimePoints,i) = 1;
end
Regressors = [Regressors,SubjectRegressors];

% Regressors = [Regressors,SubjectRegressors,[age;age],[sex;sex],[edu;edu]];

% Regressors = [Regressors,SubjectRegressors];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;

for i=1:size(dataall,1);
    DependentVariable=dataall(i,:)';
  [~,r,SSE,SSR, ~, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,[Regressors],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
        
end

Df1=size(Regressors,1)-size(Regressors,2);
T=TF_ForContrast';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;
Ttemp=T(mask~=0);
df=Df1;
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if isempty(P)
    P=0.01
    Thresholded(find(PMap<=P))=1;
else
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,1);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T.*Thresholdedmask;


TF2_ForContrast_brain_H4=[];
for j=1:8
   
     
              DependentVariable=mean(dataall(yeoindex==j,:))';
               [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,[Regressors],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
            
              TF2_ForContrast_brain_H4(j)=TF_ForContrast;

end
TF2_ForContrast_brain_H4(isnan(TF2_ForContrast_brain_H4))=0;



dataall=dataall(yeo_17labels,:);

TF2_ForContrast_brain_H4=[];
for j=1:18
   
     
              DependentVariable=mean(dataall(yeoindex17==j,:))';
               [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,[Regressors],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
            
              TF2_ForContrast_brain_H4(j)=TF_ForContrast;

end
TF2_ForContrast_brain_H4(isnan(TF2_ForContrast_brain_H4))=0;


load('H:\TSA_testCamcan\FC\covariance\groupnew.mat')
roiindex=[];
yeoindex=[groupnew;ones(16,1)*8];
yeoindex_new=[];
for i=1:8
[ind1,ind2]=find(yeoindex==i);
roiindex=[roiindex;ind1];
yeoindex_new=[yeoindex_new;ones(size(ind1,1),1)*i];
end

new=TF_ForContrast';
mean_grad=T_mci';

scatter(new,mean_grad,60,yeoindex,'filled');
yeoindex_netcolor=[(othercolor('Spectral8',8))];

colorbar('Ticks',[])
col=yeoindex_netcolor;
colormap(col);
set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold')
set(gca,'LineWidth',2);
hold on;
x=new';
y=mean_grad';
a=[x;y];
a1=a(1,:);
[a1,pos]=sort(a1);
a2=a(2,pos);
[p,s]=polyfit(a1,a2,1);
y1=polyval(p,a1);
% 95% confidence interval 计算：
[yfit,dy] = polyconf(p,a1,s,'predopt','curve');
% [yfit,dy] = polyconf(p,a1,s,'alpha',0.05);
hold on
%置信区域绘制
fill([a1,fliplr(a1)],[yfit-dy,fliplr(yfit+dy)],[0.75,0.75,0.75],'EdgeColor','none','FaceAlpha',0.5);
p=polyfit(new,mean_grad,1);
yfit=polyval(p,new);
plot(new,yfit,'k-','Markersize',3.5,'LineWidth',4)
set(gcf,'color','w')



%%
clear all
load('I:\TMS\Nii\new_results\sub.mat');
name1d1=name1;
name1d2=name1;
name1d1{6}='chenkaiju';
name1d2{6}='chenkaiju';
name1d1{14}='fang_mei';
name1d2{11}='jiachuanhua';

name2d1=name2;
name2d2=name2;
name2d1{14}='zhaoxiaocui';

BOLDd1=dir('I:\TMS\Nii\d1\fun_preprocess\BOLD_smooth100\*.mat');
BOLDd2=dir('I:\TMS\Nii\d2\fun_preprocess\BOLD_smooth100\*.mat');

d1sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{BOLDd1.name},'UniformOutput',false),...
    'UniformOutput',false);
d1sub{8}='fang_mei';

d2sub=cellfun(@(name) name{1},cellfun(@(names) strsplit(names,'_'),{BOLDd2.name},'UniformOutput',false),...
    'UniformOutput',false);


[a11,b11]=ismember(name1d1,d1sub);
[a12,b12]=ismember(name1d2,d2sub);

[a21,b21]=ismember(name2d1,d1sub);
[a22,b22]=ismember(name2d2,d2sub);

b11(b11==0)=[];
b12(b12==0)=[];


BOLDs={};
num=0;

for i=1:length(b11)
    load([BOLDd1(b11(i)).folder,'\',BOLDd1(b11(i)).name]);
    
    num=num+1;
    
    BOLDs{num}=theROITimeCoursesTotal;
end
for i=1:length(b21)
    load([BOLDd1(b21(i)).folder,'\',BOLDd1(b21(i)).name]);
    
    num=num+1;
    
    BOLDs{num}=theROITimeCoursesTotal;
end

for i=1:length(b12)
    load([BOLDd2(b12(i)).folder,'\',BOLDd2(b12(i)).name]);
    
    num=num+1;
    
    BOLDs{num}=theROITimeCoursesTotal;
end
for i=1:length(b22)
    load([BOLDd2(b22(i)).folder,'\',BOLDd2(b22(i)).name]);
    
    num=num+1;
    
    BOLDs{num}=theROITimeCoursesTotal;
end

poscm = colormap_tor([0.96 0.41 0], [1,1,0],'n', 100);  % warm
negcm = colormap_tor([.23 1 1], [0.11,0.46,1],'n', 100);  % cools
poscolor=poscm;   negcolor=negcm;

poscm = colormap_tor([1,0,0], [1 1 0],'n', 100);  % warm
negcm = colormap_tor([0,0,1], [0 1 128/255],'n', 100);  % cools
poscolor=poscm;   negcolor=negcm;

% poscm = colormap_tor([0,255,253]./255, [0,61,252]./255,'n', 100);  % warm
% negcm = colormap_tor([255,248,13]./255, [0 1 128/255],'n', 100);  % cools
% poscolor=poscm;   negcolor=negcm;
% cmp_pos=bw(1:32,:);
% cmp_neg=bw(33:end,:);
% poscm = interp1(linspace(0,1,size(cmp_pos,1)),cmp_pos,linspace(0,1,100),'cubic');
% negcm = interp1(linspace(0,1,size(cmp_neg,1)),cmp_neg,linspace(0,1,100),'cubic');
% poscolor=poscm;   negcolor=negcm;
poscolor=ccs_mkcolormap('hotcolors.tif',100);
negcolor=ccs_mkcolormap('coldcolors.tif',100);

for j=1:11
    temp=zeros(416,1);
for i=1:18
    
    temp(yeoindex17==i)=allw(i,j);
    
end
plot_hemisphere_zhengnew(temp,{left_surface,right_surface},brainmask,1,poscolor,negcolor);
end

allw=squeeze(temp(4,:,:));
load('H:\light_therapysyy\baselinall\controal\asl_enriched\spin.mat');
for j=1:11
%    close all;
    plot_hemisphere_zhengnew(allw(:,j),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
%     print(gcf,['I:\ADNI_DMN_state1_',num2str(j),'.png'],'-dpng','-r600');

% real_data=allw(:,j);
% pos_real_data=real_data;
% pos_real_data(pos_real_data<0)=0;
% pos_real_data(pos_real_data>0)=1;
% 
% neg_real_data=real_data;
% neg_real_data(neg_real_data>0)=0;
% neg_real_data(neg_real_data<0)=1;
% 
% 
% for i=1:18
%    pos_real_dic(i)=sum(yeoindex17(find(pos_real_data))==i)./sum(yeoindex17==i);
%     
%    neg_real_dic(i)=sum(yeoindex17(find(neg_real_data))==i)./sum(yeoindex17==i);
% end
% 
% for per=1:5000
%  perm_data=real_data(perm_id(:,per));
% pos_perm_data= perm_data;
% pos_perm_data(pos_perm_data<0)=0;
% pos_perm_data(pos_perm_data>0)=1;
% 
% neg_perm_data=perm_data;
% neg_perm_data(neg_perm_data>0)=0;
% neg_perm_data(neg_perm_data<0)=1;
%     
% for i=1:18
%    pos_perm_dic(i,per)=sum(yeoindex17(find(pos_perm_data))==i)./sum(yeoindex17==i);
%     
%    neg_perm_dic(i,per)=sum(yeoindex17(find(neg_perm_data))==i)./sum(yeoindex17==i);
% end
% end
% 
% for i=1:18
%  up=prctile(pos_perm_dic(i,:),95);
% 
% 
% if pos_real_dic(i)>=up
%       pos_p_percent(i)=1-length(find(pos_real_dic(i)>=pos_perm_dic(i,:)))/5000;
% 
% else
%    pos_p_percent(i)=0.06;
% end
% end
% 
% for i=1:18
%  up=prctile(neg_perm_dic(i,:),95);
% 
% 
% if neg_real_dic(i)>=up
%       neg_p_percent(i)=1-length(find(neg_real_dic(i)>=neg_perm_dic(i,:)))/5000;
% 
% else
%     neg_p_percent(i)=0.06;
% end
% end
% find(pos_p_percent<0.05)
% p_percent_new=pval_adjust(pos_p_percent,'BH');
% find(p_percent_new<0.05)
% 
% find(neg_p_percent<0.05)
% p_percent_new=pval_adjust(neg_p_percent,'BH');
% find(p_percent_new<0.05)



end


for j=1:11
clear pos_real_dic pos_perm_dic pos_p_percent neg_real_dic neg_perm_dic neg_p_percent
real_data=allw(:,j);
real_data=real_data(yeo_17labels)
pos_real_data=real_data;
pos_real_data(pos_real_data<0)=0;
pos_real_data(pos_real_data>0)=1;

neg_real_data=real_data;
neg_real_data(neg_real_data>0)=0;
neg_real_data(neg_real_data<0)=1;

pos_real_data=real_data;
pos_real_data(pos_real_data<0)=0;
pos_real_data(pos_real_data>0)=1;

neg_real_data=real_data;
neg_real_data(neg_real_data>0)=0;
neg_real_data(neg_real_data<0)=1;


for i=1:8
   pos_real_dic(i)=sum(yeoindex(find(pos_real_data))==i)./sum(yeoindex==i);
    
   neg_real_dic(i)=sum(yeoindex(find(neg_real_data))==i)./sum(yeoindex==i);
end

for per=1:5000
 perm_data=real_data(perm_id(:,per));
pos_perm_data= perm_data;
pos_perm_data(pos_perm_data<0)=0;
pos_perm_data(pos_perm_data>0)=1;

neg_perm_data=perm_data;
neg_perm_data(neg_perm_data>0)=0;
neg_perm_data(neg_perm_data<0)=1;
    
for i=1:8
   pos_perm_dic(i,per)=sum(yeoindex(find(pos_perm_data))==i)./sum(yeoindex==i);
    
   neg_perm_dic(i,per)=sum(yeoindex(find(neg_perm_data))==i)./sum(yeoindex==i);
end
end

for i=1:8
 up=prctile(pos_perm_dic(i,:),95);


if pos_real_dic(i)>=up
      pos_p_percent(i)=1-length(find(pos_real_dic(i)>=pos_perm_dic(i,:)))/5000;

else
   pos_p_percent(i)=0.06;
end
end

for i=1:8
 up=prctile(neg_perm_dic(i,:),95);


if neg_real_dic(i)>=up
      neg_p_percent(i)=1-length(find(neg_real_dic(i)>=neg_perm_dic(i,:)))/5000;

else
    neg_p_percent(i)=0.06;
end
end
find(pos_p_percent<0.05)
p_percent_new=pval_adjust(pos_p_percent,'BH');
find(p_percent_new<0.05)

find(neg_p_percent<0.05)
p_percent_new=pval_adjust(neg_p_percent,'BH');
find(p_percent_new<0.05)


end

for j=1:11
clear pos_real_dic pos_perm_dic pos_p_percent neg_real_dic neg_perm_dic neg_p_percent
real_data=allw(:,j);
real_data=real_data;
pos_real_data=real_data;
pos_real_data(pos_real_data<0)=0;
pos_real_data(pos_real_data>0)=1;

neg_real_data=real_data;
neg_real_data(neg_real_data>0)=0;
neg_real_data(neg_real_data<0)=1;

pos_real_data=real_data;
pos_real_data(pos_real_data<0)=0;
pos_real_data(pos_real_data>0)=1;

neg_real_data=real_data;
neg_real_data(neg_real_data>0)=0;
neg_real_data(neg_real_data<0)=1;


for i=1:18
   pos_real_dic(i)=sum(yeoindex17(find(pos_real_data))==i)./sum(yeoindex17==i);
    
   neg_real_dic(i)=sum(yeoindex17(find(neg_real_data))==i)./sum(yeoindex17==i);
end

for per=1:5000
 perm_data=real_data(perm_id(:,per));
pos_perm_data= perm_data;
pos_perm_data(pos_perm_data<0)=0;
pos_perm_data(pos_perm_data>0)=1;

neg_perm_data=perm_data;
neg_perm_data(neg_perm_data>0)=0;
neg_perm_data(neg_perm_data<0)=1;
    
for i=1:18
   pos_perm_dic(i,per)=sum(yeoindex17(find(pos_perm_data))==i)./sum(yeoindex17==i);
    
   neg_perm_dic(i,per)=sum(yeoindex17(find(neg_perm_data))==i)./sum(yeoindex17==i);
end
end

for i=1:18
 up=prctile(pos_perm_dic(i,:),95);


if pos_real_dic(i)>=up
      pos_p_percent(i)=1-length(find(pos_real_dic(i)>=pos_perm_dic(i,:)))/5000;

else
   pos_p_percent(i)=0.06;
end
end

for i=1:18
 up=prctile(neg_perm_dic(i,:),95);


if neg_real_dic(i)>=up
      neg_p_percent(i)=1-length(find(neg_real_dic(i)>=neg_perm_dic(i,:)))/5000;

else
    neg_p_percent(i)=0.06;
end
end


find(pos_p_percent<0.05)
p_percent_new=pval_adjust(pos_p_percent,'BH');
find(p_percent_new<0.05)

find(neg_p_percent<0.05)
p_percent_new=pval_adjust(neg_p_percent,'BH');
find(p_percent_new<0.05)


end



















clear;clc

times=[-15:3:15];
nFrames=32;
mov(1:nFrames) = struct('cdata', [],'colormap', []);
namell=dir([ '*' '.png']);%读取需要合并的图片名字    '*' '.png'表示读取文件夹中任意字段+.png格式的图片

for i=1:11
   Img=imread(['I:\ADNI_DMN_state1_',num2str(i),'.png']);
     imshow(Img,[]);
     title(num2str(times(i)),'FontSize',20);
     set(gcf,'color','w');
    frame=getframe(gcf);
    im=frame2im(frame);%制作gif文件，图像必须是index索引图像
    [I,map]=rgb2ind(im,256);
    mov(i)=getframe(gcf);

    if(i==1)
        imwrite(I,map,'DMN_s1.gif','DelayTime',0.8,'LoopCount',Inf)
    else
        imwrite(I,map,'DMN_s1.gif','WriteMode','append','DelayTime',0.8)    
    end
end
[Yeo,Header]=y_Read('D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii');
data=zeros(size(Yeo,1),size(Yeo,2),size(Yeo,3));
for j=1:11
    data=zeros(size(Yeo,1),size(Yeo,2),size(Yeo,3));
for i=1:400
    data(Yeo==i)=allw(i,j);
    
end
y_Write(data,Header,['C:\Users\dell\Desktop\test\',num2str(j),'.nii'])
end

unique(mask(:));
%%
%plot
left_surface=gifti('D:\BrainSpace-master\shared\surfaces\conte69_32k_left_hemisphere.gii');
right_surface=gifti('D:\BrainSpace-master\shared\surfaces\conte69_32k_right_hemisphere.gii');
surface={left_surface,right_surface};

for net=1:18
    load(['I:\ADNI\lag\Summary_measures_yeo17new_',num2str(net),'.mat'])
    
    for state=1:4
        close all;
        h=[];
        h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.3,1]);
       allw=squeeze(temp(state,:,:)); 
            for j=1:11

                   T=allw(:,j);
                   T=T(yeo_17labels);
                    T(isnan(T))=0;
                    T_pos=T;
                    T_pos(T<=0)=0;
                    T_neg=T;
                    T_neg(T>=0)=0;

                    resultMap=zeros(size(brainmask,1),1);
                    for i=1:400
                        resultMap(brainmask==i)=T(i);
                    end
                    resultsubMap=zeros(16,1);
                    for i=401:416
                        resultsubMap(i-400)=T(i);
                    end

                    maxpos=max(T_pos);
                    minpos=min(T_pos(T_pos~=0));
                    minneg=min(T_neg);
                    maxneg=max(T_neg(T_neg~=0));
                    resultcol=zeros(size(resultMap,1),3);
                    if  maxpos==minpos
                        if find(T==maxpos)<=400
                            flagpos=1;
                        else
                            flagpos=2;
                        end
                        idx_pos=100;
                        pos_step=1;
            %                         poscolor=(othercolor('YlOrRd9',100)); 
                        resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

                    else
                        flagpos=3;
                        pos_step=((maxpos)-(minpos))/100;
                        idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                        idx_pos(idx_pos<1)=1;
                        idx_pos(idx_pos>100)=100;
            %                         poscolor=(othercolor('YlOrRd9',100));  
                        resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
                    end

                   if  maxneg==minneg
                        if find(T==maxneg)<=400
                            flagneg=1;
                        else
                            flagneg=2;
                        end
                        idx_neg=1;
                        neg_step=1;
                         %negcolor=(othercolor('YlGnBu9',100));
                        resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

                    else
                        flagneg=3;
                        neg_step=((-minneg)-(-maxneg))/100;
                        idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                        idx_neg(idx_neg<1)=1;
                        idx_neg(idx_neg>100)=100;
                         %negcolor=(othercolor('YlGnBu9',100));
                        resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
                   end
                    resultcol(resultMap==0,:)=repmat([0.9 0.9 0.9],sum(resultMap==0),1);


                    S = convert_surface(surface);
                    data=resultcol;
                    D={};
                    for ii = 1:2
                    D{1,ii} = data(1:size(S{1}.coord,2),:);
                    if numel(S) == 2
                        D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
                    end
                    end
                    jj=1;
                        for ii = 1:numel(S)*2
                            idx = ceil(ii/2);
                            h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-.09*j .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                            h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                        ,'FaceColor', 'interp','EdgeColor', 'none');

                            material dull; lighting phong;
                        end
                        set(h.axes(:,jj)                    , ...
                        'Visible'           , 'off'         , ...
                        'DataAspectRatio'   , [1 1 1]       , ...
                        'PlotBoxAspectRatio', [1 1 1]      )
                    set(h.axes(1,:),'View',[-90 0]);
                    set(h.axes(2,:),'View',[90 0]);
                    set(h.axes(3,:),'View',[-90 0]);
                    set(h.axes(4,:),'View',[90 0]);


                    for ii = 1:4
                        axes(h.axes(ii));
                        h.camlight(ii) = camlight();
                    end

                   idx_subcluster = find(resultsubMap ~= 0);
                    idx_subcluster_out=setdiff([1:16],idx_subcluster);

                    for ii=5:8
                        for jj=1
                            h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-.09*j .1 .1]);

                            for i=1:length(idx_subcluster)
                                thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                                temp1=thala.dat;
                                temp1(temp1~=idx_subcluster(i))=0;
                                thala.dat=temp1;
                                thalaregion=region(thala);
                                if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                                    subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                                    subidx_pos(subidx_pos<1)=1;
                                    subidx_pos(subidx_pos>100)=100;
                                    p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                                elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                                    p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                                elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                                    subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                                    subidx_neg(subidx_neg<1)=1;
                                    subidx_neg(subidx_neg>100)=100;
                                    p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                                elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                                    p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                                end


                            end

                          material dull; lighting phong;
                        set(h.axes(:,jj)               , ...
                        'Visible'           , 'off'         , ...
                        'DataAspectRatio'   , [1 1 1]       , ...
                        'PlotBoxAspectRatio', [1 1 1]  ); 
                        end
                    end
                    set(h.axes(5,1),'View',[-90 0]);
                    set(h.axes(6,1),'View',[90 0]);     
                    set(h.axes(7,1),'View',[45 0]);    
                    set(h.axes(8,1),'View',[-45 0]);  
                    material dull; lighting phong;
                        set(gca                 , ...
                        'Visible'           , 'off'         , ...
                        'DataAspectRatio'   , [1 1 1]       , ...
                        'PlotBoxAspectRatio', [1 1 1]  );
                    set(gca,'color','w');
                    set(gcf,'color','w');


                
            end
            print(gcf,['I:\ADNI\lag\yeonet_',num2str(net),'state',num2str(state),'.png'],'-dpng','-r600');

    end
end

cent=csvread('D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Centroid_coordinates\Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
imagepath='D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii';
vimage=spm_vol(imagepath);
image=spm_read_vols(vimage);

Mapimage = inv(vimage.mat);
for i=1:400
x= round(Mapimage(1,1)*cent(i,2)+Mapimage(1,4)); 
y= round(Mapimage(2,2)*cent(i,3)+Mapimage(2,4));
z= round(Mapimage(3,3)*cent(i,4)+Mapimage(3,4)); 
yeo_17labels(i)=image(x,y,z);
end


cent=csvread('D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Centroid_coordinates\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
imagepath='D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii';
vimage=spm_vol(imagepath);
image=spm_read_vols(vimage);

Mapimage = inv(vimage.mat);
for i=1:400
x= round(Mapimage(1,1)*cent(i,2)+Mapimage(1,4)); 
y= round(Mapimage(2,2)*cent(i,3)+Mapimage(2,4));
z= round(Mapimage(3,3)*cent(i,4)+Mapimage(3,4)); 
yeo_17labels(i)=image(x,y,z);
end
yeo_17labels=[yeo_17labels,[401:416]];       

allw=squeeze(temp(4,:,:));
cent=csvread('D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Centroid_coordinates\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
imagepath='D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii';
vimage=spm_vol(imagepath);
image=spm_read_vols(vimage);

Mapimage = inv(vimage.mat);
for i=1:400
x= round(Mapimage(1,1)*cent(i,2)+Mapimage(1,4)); 
y= round(Mapimage(2,2)*cent(i,3)+Mapimage(2,4));
z= round(Mapimage(3,3)*cent(i,4)+Mapimage(3,4)); 
yeo_17labels(i)=image(x,y,z);
end
yeo_17labels=[yeo_17labels,[401:416]]; 

for lags=1:11
T=allw(:,lags);
T=T(yeo_17labels);

resultMap=zeros(size(brainmask,1),1);
for i=1:400
    resultMap(brainmask==i)=T(i);
end
cortex=spm_read_vols(spm_vol('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii'));
cortex(isnan(cortex))=0;
new=zeros(91,109,91);
for i=401:416
    new(cortex==i-400)=T(i);
end
lefts = read_surface('conte69_32k_left_hemisphere.gii');
lefts.coord=lefts.vertices';
lefts.tri=lefts.faces;

rights = read_surface('conte69_32k_right_hemisphere.gii');
rights.coord=rights.vertices';
rights.tri=rights.faces;
RTmap=[];
RTmap.data=new;
RTmap.origin = [-90 -126 -72];
RTmap.vox    = [2 2 2];
RTmapsleft=SurfStatVol2Surf(RTmap,lefts); %vol->surf  vol.region(-90 -126 -72)
RTmapsright=SurfStatVol2Surf(RTmap,rights);
sub=[RTmapsleft;RTmapsright];
resultMap(sub~=0)=sub(sub~=0);

[correlationall(:,lags), featureall(lags,:)] = meta_analytic_decoder_zzh(resultMap, ...
'template', 'fslr32k','data_dir','I:\ADNI\lag\neurosynth');
end

%step1 
for net=1:18
load(['I:\ADNI\lag\Summary_measures_yeo17new_',num2str(net),'.mat'])
data=reshape(temp_all,416,11,81,334);

for sub=1:size(data,4)
    for lag=1:size(data,2)
       
        data_sub=squeeze(data(:,lag,:,sub));
        states_lag=squeeze(temp(:,:,lag)); 
        [~, state_assignment] = min(pdist2(data_sub', states_lag), [], 2);
        
       
            S=size(temp,1);
            sub_st=state_assignment;
            
       
      
            for s=1:S
                fo(sub,s,lag)=100*nnz(state_assignment==s)/size(state_assignment,1);
            end
   

            for k=1:S
                loc=(sub_st==k);
                try
                    cnt=diff([0;find(diff(loc));numel(loc)]);
                catch
                    keyboard;
                end
                if loc(1)==1
                    life_time{sub,k}=cnt(1:2:end);  
                else
                    life_time{sub,k}=cnt(2:2:end);
                end
                avg_life(sub,k,lag)=mean(life_time{sub,k});
            end
            
              emp=zeros(S,S);
              stat=state_assignment;
               for j=1:S
                    for k=1:S
                        for q=2:length(stat)
                            if stat(q)==k && stat(q-1)==j
                                emp(j,k)=emp(j,k)+1;
                            end
                        end
                    end
                end
            
                for r=1:S
                    emp(r,:)=emp(r,:)/sum(emp(r,:)); % to make sure every row adds to 1
                    % i.e., total probability = 1
                end
                emp(isnan(emp))=0; 
                row_has_all_zeros = ~any(emp, 2);
                indices = find(row_has_all_zeros);
                for k=1:length(indices)
                    emp(indices(k),indices(k))=1;
                end
                subj_emp(lag,:,:,sub)=emp;
                  
        
    end
    
end

avg_life(isnan(avg_life))=0;
fo(isnan(fo))=0;

% save(['Summary_measures_yeo17new_',num2str(net),'.mat'],'temp','temp_all','avg_life','fo','idx');

[d1,d2,d3]=size(avg_life);
for i=1:d1
    for j=1:d3
        temp_life(i,:,j)=zscore(avg_life(i,:,j));
        
        
    end
end
avg_life=temp_life;

[d1,d2,d3]=size(fo);
for i=1:d1
    for j=1:d3
        temp_fo(i,:,j)=zscore(fo(i,:,j));
        
        
    end
end
fo=temp_fo;



Contrast = zeros(1,size(Regress_Simens,2));
Contrast(1) = 1;

[d1,d2,d3]=size(avg_life);
for i=1:d2;
    for j=1:d3
    DependentVariable=squeeze(avg_life(:,i,j));
  [~,r,SSE,SSR, ~, TF_ForContrast(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T'); 
    end
end

TF_ForContrast_all_avg(net,:,:)=TF_ForContrast;


[d1,d2,d3]=size(fo);
for i=1:d2;
    for j=1:d3
    DependentVariable=squeeze(fo(:,i,j));
  [~,r,SSE,SSR, ~, TF_ForContrast(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T'); 
    end
end
TF_ForContrast_all_fo(net,:,:)=TF_ForContrast;

% 
% [d1,d2,d3,d4]=size(subj_emp);
% for k=1:d1
%     for i=1:d2
%         for j=1:d3
%         DependentVariable=squeeze(subj_emp(k,i,j,:));
%       [~,r,SSE,SSR, ~, TF_ForContrast(i,j,k), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T'); 
%         end
%     end
% end
% TF_ForContrast_emp=TF_ForContrast;



% Df1=size(Regress_Simens,1)-size(Regress_Simens,2);
% T=TF_ForContrast;
% 
% 
% T(isnan(T))=0;
% mask=T;
% mask(mask~=0)=1;
% Ttemp=T(mask~=0);
% df=Df1;

% PMap=2*(1-tcdf(abs(Ttemp), df));
% qThreshold=0.05;%设置p值
% SortP=sort(PMap);
% V=length(SortP);
% I=(1:V)';
% cVID = 1;
% cVN  = sum(1./(1:V));
% P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
% Thresholded=zeros(size(PMap));
% if isempty(P)
%     P=0.05
%     Thresholded(find(PMap<=P))=1;
% else
%     Thresholded(find(PMap<=P))=1;
% end
% Thresholdedmask=zeros(4,11);
% Thresholdedmask(mask~=0)=Thresholded;
% T_fdr=T.*Thresholdedmask;
% figure,imagesc(T_fdr);
% colormap(bw);


% PMap=2*(1-tcdf(abs(TF_ForContrast), df));




end
for net=1:18
load(['I:\ADNI\lag\Summary_measures_yeo17new_',num2str(net),'.mat'])
start=[1:81:27054];
clear avg_life fo subj_emp TF_ForContrast
for sub=1:334
S=size(temp,1);
sub_st=idx(start(sub):start(sub)+81-1);
for s=1:S
fo(sub,s)=100*nnz(sub_st==s)/size(sub_st,1);
end
for k=1:S
loc=(sub_st==k);
try
cnt=diff([0;find(diff(loc));numel(loc)]);
catch
keyboard;
end
if loc(1)==1
life_time{sub,k}=cnt(1:2:end);
else
life_time{sub,k}=cnt(2:2:end);
end
avg_life(sub,k)=mean(life_time{sub,k});
end
emp=zeros(S,S);
stat=sub_st;
for j=1:S
for k=1:S
for q=2:length(stat)
if stat(q)==k && stat(q-1)==j
emp(j,k)=emp(j,k)+1;
end
end
end
end
for r=1:S
emp(r,:)=emp(r,:)/sum(emp(r,:)); % to make sure every row adds to 1
% i.e., total probability = 1
end
emp(isnan(emp))=0;
row_has_all_zeros = ~any(emp, 2);
indices = find(row_has_all_zeros);
for k=1:length(indices)
emp(indices(k),indices(k))=1;
end

emp_temp=zscore(emp(:));
emp=reshape(emp_temp,4,4);

subj_emp(:,:,sub)=emp;
end

Contrast = zeros(1,size(Regress_Simens,2));
Contrast(1) = 1;
[d1,d2,d3]=size(avg_life);
avg_life(isnan(avg_life))=0;
fo(isnan(fo))=0;

avg_life=zscore(avg_life')';
fo=zscore(fo')';
for i=1:4
DependentVariable=squeeze(fo(:,i));
[~,r,SSE,SSR, ~, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T');
x{1,1}=squeeze(fo(Regress_Simens(:,1)==1,i));
x{1,2}=squeeze(fo(Regress_Simens(:,1)==-1,i));
boxplot_wani_2016(x,'color',col,'linewidth', 2, 'boxlinewidth', 3, 'violin', 'dots');
print(gcf,['I:\ADNI\lag\yeonet_',num2str(net),'state',num2str(i),'_fo.png'],'-dpng','-r600');
close all;
end


PMap=2*(1-tcdf(abs(TF_ForContrast), df))
for i=1:4
DependentVariable=squeeze(avg_life(:,i));
[~,r,SSE,SSR, ~, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T');
x{1,1}=squeeze(avg_life(Regress_Simens(:,1)==1,i));
x{1,2}=squeeze(avg_life(Regress_Simens(:,1)==-1,i));
boxplot_wani_2016(x,'color',col,'linewidth', 2, 'boxlinewidth', 3, 'violin', 'dots');
print(gcf,['I:\ADNI\lag\yeonet_',num2str(net),'state',num2str(i),'_avglife.png'],'-dpng','-r600');
close all;
end


PMap=2*(1-tcdf(abs(TF_ForContrast), df))


for i=1:4
    for j=1:4
        DependentVariable=squeeze(subj_emp(i,j,:));
[~,r,SSE,SSR, ~, TF_ForContrast(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T');

    end
end
PMap=2*(1-tcdf(abs(TF_ForContrast), df))
end

[d1,d2,d3]=size(avg_life);
for i=1:d1
    for j=1:d3
        temp_life(i,:,j)=zscore(avg_life(i,:,j));
        
        
    end
end
avg_life=temp_life;

[d1,d2,d3]=size(fo);
for i=1:d1
    for j=1:d3
        temp_fo(i,:,j)=zscore(fo(i,:,j));
        
        
    end
end
fo=temp_fo;
[d1,d2,d3]=size(avg_life);
for i=1:d2;
    for j=1:d3
    DependentVariable=squeeze(avg_life(:,i,j));
  [~,r,SSE,SSR, ~, TF_ForContrast(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T'); 
    end
end


[d1,d2,d3]=size(fo);
for i=1:d2;
    for j=1:d3
    DependentVariable=squeeze(fo(:,i,j));
  [~,r,SSE,SSR, ~, TF_ForContrast(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regress_Simens,Contrast,'T'); 
    end
end
TF_ForContrast_all_fo(net,:,:)=TF_ForContrast;



for net=9:18
    load(['I:\ADNI\lag\Summary_measures_yeo17new_',num2str(net),'.mat'])
    
    states(net,:)=reshape(temp,1,[]);
end
distMatrix = pdist(states, 'euclidean'); 
fullDistMatrix = squareform(distMatrix);
Z = linkage(fullDistMatrix, 'ward');
figure;
dendrogram(Z);
statesnew=[];
for net=1:18

    statesnew=[statesnew;reshape(reshape(states(net,:),4,416,11),4,[])];
    
end
distMatrix = pdist(statesnew, 'euclidean'); 
fullDistMatrix = squareform(distMatrix);
Z = linkage(fullDistMatrix, 'ward');
figure;
dendrogram(Z);


x{1,1}=squeeze(avg_life(Regress_Simens(:,1)==1,i));
x{1,2}=squeeze(avg_life(Regress_Simens(:,1)==-1,i));
col=[253,185,106;78,98,171]./255;
al_goodplot(squeeze(avg_life(Regress_Simens(:,1)==1,i)),[],0.5,[253,185,106]./255,'right')
al_goodplot(squeeze(avg_life(Regress_Simens(:,1)==-1,i)),[],0.5,[78,98,171]./255,'left')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.1, 0.5])
set(gcf,'color','w');
set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','bold');


cent=csvread('D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Centroid_coordinates\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
imagepath='D:\CBIG-0.17.0-Fix_Absolute_Path\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\MNI\Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii';
vimage=spm_vol(imagepath);
image=spm_read_vols(vimage);

Mapimage = inv(vimage.mat);
for i=1:400
x= round(Mapimage(1,1)*cent(i,2)+Mapimage(1,4)); 
y= round(Mapimage(2,2)*cent(i,3)+Mapimage(2,4));
z= round(Mapimage(3,3)*cent(i,4)+Mapimage(3,4)); 
yeo_17labels(i)=image(x,y,z);
end
yeo_17labels=[yeo_17labels,[401:416]];

for net=1:18
     load(['I:\ADNI\lag\Summary_measures_yeo17new_',num2str(net),'.mat'])
     for state=1:4
        allw=squeeze(temp(state,:,:));
        for lags=1:11
        T=allw(:,lags);
        T=T(yeo_17labels);

        resultMap=zeros(size(brainmask,1),1);
        for i=1:400
            resultMap(brainmask==i)=T(i);
        end
        cortex=spm_read_vols(spm_vol('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii'));
        cortex(isnan(cortex))=0;
        new=zeros(91,109,91);
        for i=401:416
            new(cortex==i-400)=T(i);
        end
        lefts = read_surface('conte69_32k_left_hemisphere.gii');
        lefts.coord=lefts.vertices';
        lefts.tri=lefts.faces;

        rights = read_surface('conte69_32k_right_hemisphere.gii');
        rights.coord=rights.vertices';
        rights.tri=rights.faces;
        RTmap=[];
        RTmap.data=new;
        RTmap.origin = [-90 -126 -72];
        RTmap.vox    = [2 2 2];
        RTmapsleft=SurfStatVol2Surf(RTmap,lefts); %vol->surf  vol.region(-90 -126 -72)
        RTmapsright=SurfStatVol2Surf(RTmap,rights);
        sub=[RTmapsleft;RTmapsright];
        resultMap(sub~=0)=sub(sub~=0);

        [correlationall(net,state,lags,:), featureall(net,state,lags,:)] = meta_analytic_decoder_zzh(resultMap, ...
        'template', 'fslr32k','data_dir','I:\ADNI\lag\neurosynth');
        end
     end
end


%%
[d1,H1]=y_Read('E:\WorkSpace\MRI_MASK\AAL3v2.nii');
[d2,H2]=y_Read('E:\WorkSpace\MRI_MASK\r2schaefer400MNI.nii');

temp=(d2(d1==5));
temp(temp==0)=[];
nums=unique(temp);
for i=1:length(unique(temp))
nums(i,2)=length(find(temp==nums(i,1)));
end

sc=readNPY('I:\TMS\Nii\d1\Conn_matrix\temp\SC_con.npy');
sc(isnan(sc))=0;


index=B1(19:25,1);
N=416;
pulse_value = 0.9;
pulse_duration = 5;     % 秒
interburst_interval = 15;  % 秒
n_bursts = 30;
T = 900;  % 15分钟
U = zeros(N,1,T);

for target_region=index
t = 1;  % 当前时间点
for burst = 1:n_bursts
    % 给这 burst 的前5秒赋值 0.9
    U(target_region,1,t : t+pulse_duration-1) = pulse_value;

    % 移动时间指针：5秒刺激 + 15秒休息
    t = t + pulse_duration + interburst_interval;
    
    % 防止越界
    if t > T - pulse_duration
        break
    end
end
end
sc = sc ./ (eigs(sc,1, 'largestabs') + 0 .* eigs(sc,1, 'largestabs') );
A =sc-eye(416);
 [ x, u, energy ] = fcn_optimalControlContinuous( A, B, 1, xi, xf, 1 );

 
 
 
 T_bold = 255;
U_new = zeros(N, 1, T_bold);

% 用前几个刺激串构造短版本 TMS 控制
pulse_duration = 5;
interburst_interval = 15;
pulse_value = 0.9;
stim_node =index;

t = 1;
while t + pulse_duration - 1 <= T_bold
    U_new(stim_node,1,t:t+pulse_duration-1) = pulse_value;
    t = t + pulse_duration + interburst_interval;
end
 
 
 
 
 
 

files1=dir('I:\TMS\Nii\d1\BOLD_temp\*.mat');
files2=dir('I:\TMS\Nii\d2\BOLD_temp\*.mat');
B = zeros(N,N);
for i=1:length(index)
    B(index(i),index(i))=1;
end
A = sc / max(eig(sc)) * 0.95-eye(416);  % 保证系统稳定

B=eye(416);

for i=1:length(files1)
    load([files1(i).folder,'\',files1(i).name])
    xi=mean(zscore(theROITimeCoursesTotal'),2);
    [x] = open_loop_control( A, B, xi, U(:,:,1:585));
    for iter=1:19
         xi=mean(zscore(squeeze(x)),2);
        [x] = open_loop_control( A, B, xi, U(:,:,1:585));
    end
    load([files2(i).folder,'\',files2(i).name])
    xf=mean(zscore(theROITimeCoursesTotal'),2);
    xff=mean(zscore(squeeze(x)),2);
    [r(i),p(i)]=corr(xff,xf);
end


for i=1:length(files1)
    load([files1(i).folder,'\',files1(i).name])
    xi=zscore(mean(theROITimeCoursesTotal)');
    [x] = open_loop_control( A, B, xi, U);
    for iter=1:19
         xi=zscore(mean(squeeze(x),2));
        [x] = open_loop_control( A, B, xi, U);
    end
    load([files2(i).folder,'\',files2(i).name])
    xf=zscore(mean(theROITimeCoursesTotal)');
    xff=zscore(mean(squeeze(x),2));
    [r(i),p(i)]=corr(xff,xf);
end




for i=1:length(files1)
    load([files1(i).folder,'\',files1(i).name])
    xi=zscore(mean(theROITimeCoursesTotal)');
    [x] = open_loop_control( A, B, xi, U_new);
    for iter=1:19
         xi=zscore(mean(squeeze(x),2));
        [x] = open_loop_control( A, B, xi, U_new);
    end
    load([files2(i).folder,'\',files2(i).name])
    xf=zscore(mean(theROITimeCoursesTotal)');
    xff=zscore(mean(squeeze(x),2));
    [r(i),p(i)]=corr(xff,xf);
end


