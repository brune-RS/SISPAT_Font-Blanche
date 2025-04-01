clear;
close all;
clc;

annee=2014;
addpath('/home/brune/SISPAT/Fonctions matlab/')
%addpath('D:\SISPAT\Fonctions matlab\')
%% Chargement des obs

%obs_FBt=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');
%obs_FBt=table2timetable(obs_FBt);
obs_FB=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');
prof2m=importdata('/home/brune/SISPAT/tests FB/prof2m.mat');
prof20m=importdata('/home/brune/SISPAT/tests FB/prof20m.mat');
path_prof20m='/home/brune/SISPAT/tests FB/prof_maillageFB.csv';
ep_20m=importdata('/home/brune/SISPAT/tests FB/ep_internoeuds_FB_20m.mat');

% obs_FBt=importdata('D:\SISPAT\\tableaux matlab\Font-Blanche\tableau_obs_FB_14_21.mat');
% obs_FBt=table2timetable(obs_FBt);
% obs_FB=importdata('D:\SISPAT\tableaux matlab\Font-Blanche\tableau_obs_FB_14_21.mat');
% prof2m=importdata('D:\SISPAT\tests FB\prof2m.mat')
%prof20m=importdata('D:\SISPAT\tests FB\prof20m.mat');
 %path_prof20m='D:\SISPAT\tests FB\prof_maillageFB.csv';
% ep_20m=importdata('D:\SISPAT\tests FB\ep_internoeuds_FB_20m.txt');
% 
% fopen(path_prof20m,'r');
%prof20m = fread(frid,'single');
%G_obs=obs_FB.G_5;
%Pcumuljour=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/Pcumuljour.mat');

%% Lecture du (des) fichier(s) atmpath_res2='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6/';
path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6s1/';
path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6/';
path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/tests/a2/';
% path_res2='D:\SISPAT\execution\sortie_modele\FB\K1d6\';
% path_res3='D:\fev\rooth2\';
% path_res4='D:\fev\mix1\';
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
    %sol=sol2;qq
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
% 
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
% figure()
% plot(dates,squeeze(sol3(2,1,:)))
% hold on 
% plot(dates,squeeze(sol3(2,100,:)))
% 

% sol2(3,1:15,:)=sol2(3,1:15,:)/0.15;
% sol2(3,16:44,:)=sol2(3,16:44,:)/0.5;
% moist_FB_ci(45:64)=moist_FB_ci(45:64)*0.4;
% moist_FB_ci(65:79)=moist_FB_ci(65:79)*0.35;
% moist_FB_ci(80:94)=moist_FB_ci(80:94)*0.3;
% moist_FB_ci(95:117)=moist_FB_ci(95:117)*0.2;
% moist_FB_ci(118:145)=moist_FB_ci(118:145)*0.15;
% moist_FB_ci(146:end)=moist_FB_ci(146:end)*0.1;

%% variables 
% Simulations : Construction des matrices z, t, h, fv, fl
        
h=squeeze(sol2(1,:,:));% potentiel hydrique
t=squeeze(sol2(2,:,:));%temp
swc=squeeze(sol2(3,:,:));%humidité
swc_corr=swc;% Humidité dans le volume du sol (cailloux compris) 
ex=squeeze(sol2(4,:,:));%extraction racinaire
fl=squeeze(sol2(5,:,:));% flux liquide
fv=squeeze(sol2(6,:,:));% flux vapeur

% Changement d'unités fl, fv sont en m/s-->mm/j    % unité a changer 
        fl=fl*(10^3*(60*60*24));
        fv=fv*(10^3*(60*60*24));
        er=ex*100*((60*60*24)); % er en 1/s-->%/j

        % Changement d'unités fl, fv sont en m/s-->mm/h
        %fl=fl*(10^3*(60*60);
        %fv=fv*(10^3*(60*60));
        %er=er*100*((60*60)); % er en 1/s-->%/j
        
%% +correction cailloux





% corrections à apporter compter l'humidité dans partie sol fin 

 % Pourcentage de cailloux FB 50,60,65,70,80,85,90,90%
 c= [0.85,0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1];  % pourcentage de sol
  
 for j=1:409 
    if j<=16
    dHcorr(j,:)=swc_corr(j,:)/c(1);
    end 
    if j>16 & j<=16+29
    dHcorr(j,:)=swc_corr(j,:)/c(2);
    end 
    if j>16+29 & j<=16+29+20
    dHcorr(j,:)=swc_corr(j,:)/c(3);
    end 
    if j>16+29+20 & j<=16+29+20+15
    dHcorr(j,:)=swc_corr(j,:)/c(4);
    end 
    if j>16+29+20+15 & j<=16+29+20+15+15
    dHcorr(j,:)=swc_corr(j,:)/c(5);
    end 
    if j>16+29+20+15+15 & j<=16+29+20+15+15+23
    dHcorr(j,:)=swc_corr(j,:)/c(6);
    end 
    if j>16+29+20+15+15+23 & j<=16+29+20+15+15+23+28
    dHcorr(j,:)=swc_corr(j,:)/c(7);
    end 
    if j>16+29+20+15+15+23+28 
    dHcorr(j,:)=swc_corr(j,:)/c(8);
    end 
end 




%dex=squeeze(sol2(4,:,:));  %extraction racinaire en m2/s dans l'espace internoeud 

%% 
figure()
imagesc(dT)
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
surf(dates,prof20m,dT)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Humidité')
colorbar
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 
     ax = gca;
     shading('interp')
view(0,90)

figure()
surf(dates,prof20m,dHcorr)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Humidité ')
colorbar
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
     shading('interp')
     set(gca,'Ydir','reverse')
view(0,90)
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

%% script jordi 
% hori=[0 2 20 70 120 400];
% ygrid=hori(1):0.01:hori(end);
% [~,~,ind]=intersect(hori,ygrid);
% 

ygrid=0:0.05:2000;  % tous les 0.5 mm
ygrid=ygrid';
[~,~,ind]=intersect(round(prof20m*20),round(ygrid*20));  % intersection maillage sol et grille 0.5mm
ygrid_bis=ygrid(ind); % ind   marque la concordance 
Humid=[];
%dex=squeeze(sol2(4,:,:)); 
dex=(ex.*ep_20m)*1000*1800;  % extrac pondéré par lépaisseur internoeud et passé en mm par 30 min
dex=dex(:,1:24863);
extrac=[];
for i=2:409
    d=ind(i-1);
    d2=ind(i);
    m=d2-d;
    %mat1h=repmat(dHcorr(i-1,:),m,1);
    mat1ex=repmat(dex(i-1,:),m,1);
    %Humid(d:d2-1,:)=mat1h;
    extrac(d:d2-1,:)=mat1ex;
    %Humid(d+1,:)=dHcorr(i-1,:);
    %extrac(d2,:)=dex(i-1,:); 
end 

% reproduit les valeurs plusieurs fois    ex: noeud de 2 à 8 cm  valeur x
% --->> noeud 2-2.05 x    2.05-2.1 x  ... 

Humid= cat(2,ygrid,Humid);
Humid=reshape(Humid,40001,2);
extrac= cat(2,extrac,ygrid);
extrac=reshape(extrac,40000,2);



ep_5cm=linspace(0,2000,2001); % tous les 5 cm 
e=1;
 extrac_mat2=[];
 extrac_mat2(201,1:24863)=0;
         

 for i=2:201
     ind1=(i-1)*20-19;
     ind2=ind1+18;
     extrac_mat2(i-1,:)=sum(extrac(ind1:ind2,:));
 end
% somme des valeurs dans les cases tous les 5 cm 



 for i=1:length(extrac)
       if ep_5cm(e)>= ygrid
           extrac_mat2(e,:)=extrac_mat2(e,:)+extrac(i,:);
       else 
           e=e+1;
           extrac_mat2(e,:)=extrac(i,2);
           i
       end
 end 


figure()
surf(dates,ep_5cm,extrac_mat2)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Humidité ')
colorbar
colorbar.label = "Température";
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
     shading('interp')
     set(gca,'Ydir','reverse')
view(0,90)
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

%% test script cecile adapté (contourf)
surf(dates,prof20m,dHcorr)
ylabel('Depth (cm)')
xlabel('Date')
zlabel('Humidité ')
f=figure;
            
% flux liquide   -3:0.05:0.5     marquer le 0   (:,1:24863)(1:24863)
% extrac    -0.02:0.0001:0.02
            %set(gcf,'Units','Normalized','Position',gcfposition);
           % set(gca,'Units','Normalized','Position',gcaposition);
            
            ax(1)=figure()
           
            %[?,?]   contourf(temps,prof, variable)
            mist=[-3:0.05:0.5];
            [c,hp]=contourf(dates2,prof20m,fl,'LevelList', -4:0.01:2 ); % contourf pour colorer entre les contours et contour trace juste les contours
            %ytick=-log10([400 250 200 150 100 50 20 10 5 1 0.2]);
            %yticklabel=[400 250 200 150 100 50 20 10 5 1 0.2];
            set(hp,'LineColor','none');
            %set(gca,'ytick',ytick','yticklabel',yticklabel,'ylim',-log10([z(toto(end),1);z(toto(1),1)]),'FontSize',fontsize);
            ylabel('Depth (cm)');
            %title('extraction racinaire');
            colormap(flipud(jet));
            set(gca,'Ydir','reverse')
            dynamicDateTicks([],[],'dd/mm')
            dates2d=datetime(dates2, 'ConvertFrom', 'datenum');
            %ytick=[dates2d];
            %yticklabel=[dates2d];
            co=colorbar;
            cmap = colormap;       %current map
            %cmap(1,:) = [0 0 0];   %make first color white
            colormap(cmap);   
           % set(co,'fontsize',fontsize,'Ytick',0:0.02:0.24);
            title(co,'flux liquide','FontSize',9);
            set(ax(1),'CLim', [0 0.24])
            
            %figure des observations
            ax(2)=subplot(2,1,2);
            [c,hp]=contourf(tobs(1:6,toto2),-log10(zobs(1:6,toto2)),wobs(1:6,toto2),0:0.01:0.20); % contourf pour colorer entre les contours et contour trace juste les contours
            ytick=-log10([400 250 200 150 100 50 20 10 5 1 0.2]);
            yticklabel=[400 250 200 150 100 50 20 10 5 1 0.2];
            set(hp,'LineColor','none');
            set(gca,'ytick',ytick','yticklabel',yticklabel,'ylim',-log10([z(toto(end),1);z(toto(1),1)]),'FontSize',fontsize);
            ylabel('Depth (cm)');
            xlabel(sprintf('Julian Day - Year %s - %s',an,Name));
            title(sprintf('Observed soil moisture map - %s',Name));
            colormap(flipud(jet));
            co=colorbar;
            set(co,'fontsize',fontsize,'Ytick',0:0.02:0.19);
            title(co,'\theta (m^3.m^-^3)','FontSize',fontsize);
            set(ax(2),'CLim', [0 0.20])
            
            % Lier les axes des absisses entre les deux graphiques
            linkaxes(ax,'xy');
            
            % colormap(jet(20));
            nom_figure=strcat('Figures/carte_humidité_Sim_Obs_',root_parcelle,'_',an,'_',numstr,'.fig');
            saveas(f,nom_figure);
  
        





%% EXTRACTION RACINAIRE
%dex=squeeze(sol3(4,:,:));  % en m2/s dans l'espace internoeud 


som=sum(dex2(:,48*30*8+24));



%prof20m(end)% extraction racinaire  en mm/30min si ça faisiat 1 m, pondéré par la taille de la couche 
%heatmap(dates,prof,d)
%f=heatmap(d);
%f.XDisplayLabels=dates;
%f.YDisplayLabels= prof20m;
%prof20m=prof20m./100;

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

%% ____________________________________________

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