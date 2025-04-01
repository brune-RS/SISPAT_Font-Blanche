clear;
close all;
clc;
annee=2014;
%addpath('/home/brune/SISPAT/Fonctions matlab/')
addpath('D:\SISPAT\Fonctions matlab\')
% %% Chargement des obs

% obs_FBt=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');
% obs_FBt=table2timetable(obs_FBt);
% obs_FB=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');
% prof2m=importdata('/home/brune/SISPAT/tests FB/prof2m.mat');
% prof20m=importdata('/home/brune/SISPAT/tests FB/prof20m.mat');
% path_prof20m='/home/brune/SISPAT/tests FB/prof_maillageFB.csv';
% ep_20m=importdata('/home/brune/SISPAT/tests FB/ep_internoeuds_FB_20m.txt');

obs_FBt=importdata('D:\SISPAT\\tableaux matlab\Font-Blanche\tableau_obs_FB_14_21.mat');
obs_FBt=table2timetable(obs_FBt);
obs_FB=importdata('D:\SISPAT\tableaux matlab\Font-Blanche\tableau_obs_FB_14_21.mat');
prof2m=importdata('D:\SISPAT\tests FB\prof2m.mat');
prof20m=importdata('D:\SISPAT\tests FB\prof20m.mat');
path_prof20m='D:\SISPAT\tests FB\prof_maillageFB.csv';
ep_20m=importdata('D:\SISPAT\tests FB\ep_internoeuds_FB_20m.txt');

% fopen(path_prof20m,'r');
% prof20m = fread(frid,'single');
G_obs=obs_FB.G_5;
%Pcumuljour=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/Pcumuljour.mat');

%% Lecture du (des) fichier(s) atm
% path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6/';
% path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/F/';
% path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/K1d4/';
path_res2='D:\SISPAT\execution\sortie_modele\FB\K1d6\';
path_res3='D:\fev\rooth2\';
path_res4='D:\fev\mix1\';
%[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_17_21_20m_28.out']),dir([path_res3 'atm_FB_17_21_20m_28.out']),dir([path_res4 'atm_FB_17_21_20m_28.out']));
%[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_19_21_2m.out']),dir([path_res3 'atm_FB_19_21_2m.out']),dir([path_res4 'atm_FB_19_21_2m.out']));
[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'sol_FB_14_21_20m.txt']),dir([path_res3 'sol_FB_14_21_20m.txt']),dir([path_res4 'sol_FB_14_21_20m.txt']));

var1=6;
    noeuds1=409;
    file_atm2=files_atm2.name;
    fid = fopen([path_res2 file_atm2],'r');
    sol2 = fread(fid,'single');
nbtime= (length(sol2)/var1)/noeuds1; 

   sol2 = reshape(sol2,var1,noeuds1,[]);
    %sol=sol2;
    fclose(fid);

% var=6;
%     noeuds=409;
% 
%     file_atm3=files_atm3.name;
%     fid = fopen([path_res3 file_atm3],'r');
%     sol3 = fread(fid,'single');
% nbtime= (length(sol3)/var)/noeuds; 
%     %nbline= length(sol3)/86; % 86 =nombre de variables par ligne
%     sol3 = reshape(sol3,var,noeuds,[]);

%  var3=6;
%     noeuds3=409;
%     file_atm4=files_atm4.name;
%     fid = fopen([path_res4 file_atm4],'r');
%     sol4 = fread(fid,'single');
%      sol4 = reshape(sol4,var3,noeuds3,[]);
%     fclose(fid);

%end 

StartDate = datetime('2014-08-01 00:30');
EndDate =datetime('2021-12-31 23:30'); 

dates= linspace(StartDate,EndDate,12);
Date = StartDate:minutes(30):EndDate;
Date.Format = 'dd-MM-yyyy hh:mm';
dates=Date';
dates2=datenum(dates);
%tronquer les premiers mois de simulation (1er aout au 1er janvier)    
    % dates=dates(7344:end);
    % dates3=dates3(7344:end);
    % dates4=dates4(7344:end);
% X=dates(1:365*48);
% Y=linspace(1,192,192);
% Z=squeeze(sol2(2,1:365*48,:));
% figure(2)
% s=surf(Y,X,Z)
% s.EdgeColor = 'none';
 % dynamicDateTicks([],[],'dd/mm');
    % datetick('x','mmm-yy','keepticks')
figure()
plot(dates,squeeze(sol3(2,1,:)))
hold on 
plot(dates,squeeze(sol3(2,100,:)))


% sol2(3,1:15,:)=sol2(3,1:15,:)/0.15;
% sol2(3,16:44,:)=sol2(3,16:44,:)/0.5;
% moist_FB_ci(45:64)=moist_FB_ci(45:64)*0.4;
% moist_FB_ci(65:79)=moist_FB_ci(65:79)*0.35;
% moist_FB_ci(80:94)=moist_FB_ci(80:94)*0.3;
% moist_FB_ci(95:117)=moist_FB_ci(95:117)*0.2;
% moist_FB_ci(118:145)=moist_FB_ci(118:145)*0.15;
% moist_FB_ci(146:end)=moist_FB_ci(146:end)*0.1;

%% heatmap
dp=squeeze(sol2(1,:,:));% potentiel hydrique
dT=squeeze(sol2(2,:,:));%temp
dH=squeeze(sol2(3,:,:)); % Humidité quelles corrections apporter pour %cailloux ??
dex=squeeze(sol2(4,:,:));  % en m2/s dans l'espace internoeud 

dex2=(dex.*ep_20m)*1000*1800;
som=sum(dex2(:,48*30*8+24));% extraction racinaire  en mm/30min si ça faisiat 1 m, pondéré par la taille de la couche 
%heatmap(dates,prof,d)
%f=heatmap(d);
%f.XDisplayLabels=dates;
%f.YDisplayLabels= prof20m;
prof20m=prof20m./100;

%sum_ep_demi_mm=cumsum(ep_demi_mm);



ep_demi_mm=linspace(0,2.3,4600); %jusqu'à 2.3 m coupe en demi mm
 a=1;
 mat=[];
 for i=1:length(ep_demi_mm)
       if prof20m(a)>= ep_demi_mm(i)
           mat(i,:)=dex(a,:);
       else 
           a=a+1;
           mat(i,:)=dex(a,:);
       end
 end 
ep_5cm=linspace(0,2.3,46);
e=1;
 mat2=[];
 mat2(1,1:130079)=0;
 for i=1:length(ep_demi_mm)
       if ep_5cm(e)>= ep_demi_mm(i)
           mat2(e,:)=mat2(e,:)+mat(i,:);
       else 
           e=e+1;
           mat2(e,:)=mat(i,:);
           i
       end
 end 

mat2u=mat2*1000*1800;%pour des mm et 30min 
mat2u=mat2u/2000; %nombre de 0.5mm dans 1m

%%  ____________________________________________

ep_25_mm=linspace(2.3,20,709); %de 2.3 m à 20m chaque 25 mm
 a=198;
 matb=[];
 for i=1:length(ep_25_mm)
       if prof20m(a)>= ep_25_mm(i)
           matb(i,:)=dex(a,:);
       else 
           a=a+1;
           matb(i,:)=dex(a,:);
       end
 end 


ep_5cm=linspace(2.3,20,355);%5cm
e=1;
 matb2=[];
 matb2(1,:)=matb(1,:);
 for i=1:length(ep_25_mm)
       if ep_5cm(e)>= ep_25_mm(i)
          matb2(e,:)=matb2(e,:)+matb(i,:);
       else 
           e=e+1;
           matb2(e,:)=matb(i,:);
           i
       end
 end 

mat2bu=matb2*1000*1800;%pour des mm et 30min 
mat2bu=mat2bu/400; %nombre de 2.5cm dans 1m
ep_demi_mm(10)
prof20m(2)


extrac_20m=[mat2u ;mat2bu];

figure()
imagesc(dates2,prof20m',dT)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Temperature')
colorbar   
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

figure()
imagesc(dates,prof20m,dp)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Potentiel hydrique')
colorbar
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

%_______________________________________________

mat2u=mat2*1000*1800;%pour des mm et 30min 
mat2u=mat2u/2000; %nombre de 0.5mm dans 1m



ep_20m_5cm=linspace(0,20,401);

figure()
imagesc(dates2,ep_5cm,mat2bu)
ylabel('Depth (m)')
xlabel('Date')
zlabel('Extraction racinaire ')
colorbar
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 


%% Humidités dans le sol + pluies 
%start=datenum()
%dates_d=datetime(dates,'ConvertFrom','datenum');
dates=datenum(dates);
figure_handle=figure('Name','teneur en eau ')

    subplot(7,1,1)
    Precip= Pcumuljour;
    Pcumuljourtime=datenum(Pcumuljour.TIMESTAMP);
     ind_t_Ps=find(Pcumuljourtime==floor(dates(1)));
     ind_t_Pe=find(Pcumuljourtime==floor(dates(end)));
     dates_p=Pcumuljourtime(ind_t_Ps:ind_t_Pe);
     bar(dates_p,Precip.P(ind_t_Ps:ind_t_Pe))
     hold on 
     grid on 
dynamicDateTicks(); 
datetick('x','mmm-yy')
     legend('Dayly rainfall')
     ylabel('P (mm') 
    % dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
    ax.YDir = 'reverse';
    ax.FontSize = 10; 
    ax.FontWeight = 'bold'; 
   
  subplot(7,1,2)  %1cm
 plot(dates,squeeze(sol2(3,8,:)))
datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    hold on   
grid on 
 dynamicDateTicks(); 
datetick('x','mmm-yy','keepticks','keeplimits')
    title([' 1cm ' ],'FontSize',15)
     legend('Observed ','Simulated ','location','north')
       dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

 subplot (7,1,3) %5cm
 % plot(obs_FB.TIMESTAMP(ind_t_Ps:ind_t_Pe),obs_FB.SWC_TD_5_1(ind_t_Ps:ind_t_Pe),'k-','linewidth',1.5)
     hold on
    plot(dates,squeeze(sol2(3,22,:)),'b-','linewidth',1.5)
   %ylim([0.079,0.666])
    datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    hold on   
grid on 
 dynamicDateTicks(); 
datetick('x','mmm-yy','keepticks','keeplimits')
    title([' 3cm ' ],'FontSize',15)
     legend('Observed ','Simulated ','location','north')
       dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

    subplot(7,1,4) %15cm
    plot(obs_FB.TIMESTAMP(ind_t_Ps:ind_t_Pe),obs_FB.SWC_TD_15_2(ind_t_Ps:ind_t_Pe),'k-','linewidth',1.5)
    hold on
   plot(dates,squeeze(sol2(3,29,:)),'b-','linewidth',1.5)
   %ylim([0.079,0.666])
    datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    title(['5 cm' ],'FontSize',15)
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    
      subplot(7,1,5) 
    %plot(obs_FB.TIMESTAMP(ind_t_Ps:ind_t_Pe),obs_FB.SWC_TD_15_2(ind_t_Ps:ind_t_Pe),'k-','linewidth',1.5)
    hold on
   plot(dates,squeeze(sol2(3,39,:)),'b-','linewidth',1.5)
   %ylim([0.079,0.666])
    datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    title(['9 cm' ],'FontSize',15)
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

      subplot(7,1,6) 
    %plot(obs_FB.TIMESTAMP(ind_t_Ps:ind_t_Pe),obs_FB.SWC_TD_15_2(ind_t_Ps:ind_t_Pe),'k-','linewidth',1.5)
    hold on
   plot(dates,squeeze(sol2(3,51,:)),'b-','linewidth',1.5)
   %ylim([0.079,0.666])
    datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    title(['12 cm' ],'FontSize',15)
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

      subplot(7,1,7) 
    %plot(obs_FB.TIMESTAMP(ind_t_Ps:ind_t_Pe),obs_FB.SWC_TD_15_2(ind_t_Ps:ind_t_Pe),'k-','linewidth',1.5)
    hold on
   plot(dates,squeeze(sol2(3,72,:)),'b-','linewidth',1.5)
   %ylim([0.079,0.666])
    datetick('x')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    title(['25 cm' ],'FontSize',15)
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );

%%

% % Your x-axis
% x = linspace(0,2*pi);
% % y-axis
% y = linspace(0,2*pi);
% % a mesh because
% [X,Y] = meshgrid(x,y);
% % you need data at every possible x-y combination
% Z = sin(X).*atan(Y);
% % that scales the Z-values and plots a heatmap
% imagesc(x,y,Z)
% % choose a colormap of your liking
% colormap hot