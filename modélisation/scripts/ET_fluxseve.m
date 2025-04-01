%% ET.m
% Ce script lit un (ou plusieurs) fichier atm de sortie de SiSPAT et le compare aux obs de
% février- juin 2023 Brune
% pour données de FONT-BLANCHE
% graphes ET et sapflow 
clear;
close all;
clc;
 
addpath('/home/brune/SISPAT/Fonctions matlab/')

%quelle année considérée?
annee=2014;
% annee=2017;
Lv=2.484*10^6;
dt=1800;
%% Chargement des tableaux des obs
obs_FBt=importdata(['/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat']);
obs_FBt=table2timetable(obs_FBt);
obs_FB=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');

G_obs=obs_FB.G_5;
Pcumuljour=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/Pcumuljour.mat');
%% Lecture du (des) fichier(s) atm
path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/tests/a2/';
path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/tests/j1/';
path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/tests/j2/';
%[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_17_21_20m_28.out']),dir([path_res3 'atm_FB_17_21_20m_28.out']),dir([path_res4 'atm_FB_17_21_20m_28.out']));
%[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'asftm_FB_19_21_2m.out']),dir([path_res3 'atm_FB_19_21_2m.out']),dir([path_res4 'atm_FB_19_21_2m.out']));
[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_14_21_20m.out']),dir([path_res3 'atm_FB_14_21_20m.out']),dir([path_res4 'atm_FB_14_21_20m.out']));

 % Ordre des variables
    % (01)JOUR      (02)HEURE     (03)MINUTE    (04)HEURE     (05)RN SOL    (06)RN VEG
    % (07)RN TOT   F (08)H SOL     (09)H VEG     (10)H TOT     (11)LE SOL    (12)LE FEUIL
    % (13)LE VEG    (14)LE TOT    (15)G         (16)ETR/ETP   (17)TF        (18)HF
    % (19)Trad      (20)TS        (21)TAF       (22)QAF       (23)ETR TOT    (24)TR
    % (25)E SOL     (26)E FEUIL   (27)RUIS SUR  (28)PERCOL    (29)PLUIE     (30)STOCK
    % (31)Rb        (32)Ras       (33)Rsto      (34)ustar     (35)Uaf       (36)TS 5
    % (37)TS 15cm   (38)TS 25cm   (39)TS 35cm   (40)TS 50cm   (41)TS 200cm  (42)TS 500cm
    % (43)TS 1000   (44)TS 2000   (45)W 5       (46)W 15cm    (47)W 25cm    (48)W 35cm
    % (49)W 50cm    (50)W 200cm   (51)W 500cm   (52)W 1000    (53)W 2000    (54)WI 0-20cm
    % (55)WI 20-30  (56)WI 30-40  (57)WI 40-60  (58)WI 60-100 (59)WI 100-200(60)WI200-500   
    %(61)WI total   (62)Stock     (63)Q surf    (64)Albsol    (65)Albcouv   (66)Ra
    % (67)G 5cm     (68)ETP       (69)Flxliq 5                (70)Flxliq 15cm
    % (71)Flxliq25cm              (72)Flxliq35cm              (73) Flxliq50cm
    % (74)Flxliq200cm             (75)Flxliq500cm             (76)Flxliq1000cm
    % (77)Flxliq2000cm            (78)Flxvap5                 (79)Flxvap15cm
    % (80)Flxvap25cm              (81)Flxvap35cm              (82) Flxvap50cm
    % (83)Flxvap200cm             (84)Flxvap500cm             (85)Flxvap1000
    % (86)Flxvap2000cm

for i=1:length(files_atm4)
    
    file_atm2=files_atm2(i).name;
    fid = fopen([path_res2 file_atm2],'r');
    temp2 = fread(fid,'single');
    nbline= length(temp2)/86; % 86 =nombre de variables par ligne
    temp2 = reshape(temp2,[86,nbline])';
    %temp=temp2;
    fclose(fid);
    file_atm3=files_atm3(i).name;
    fid = fopen([path_res3 file_atm3],'r');
    temp3 = fread(fid,'single');
    nbline= length(temp3)/86; % 86 =nombre de variables par ligne
    temp3 = reshape(temp3,[86,nbline])';
    fclose(fid);
    file_atm4=files_atm4(i).name;
    fid = fopen([path_res4 file_atm4],'r');
    temp4 = fread(fid,'single');
    nbline= length(temp4)/86; % 86 =nombre de variables par ligne
    temp4 = reshape(temp4,[86,nbline])';
    fclose(fid);

end 
    dates=datenum(annee,8,1)-1+temp2(:,1)+(1/24)*temp2(:,2)+(1/(24*60))*temp2(:,3);
    dates3=datenum(annee,8,1)-1+temp3(:,1)+(1/24)*temp3(:,2)+(1/(24*60))*temp3(:,3);
    dates4=datenum(annee,8,1)-1+temp4(:,1)+(1/24)*temp4(:,2)+(1/(24*60))*temp4(:,3);
%tronquer les premiers mois de simulation (1er aout au 1er janvier)    
    dates=dates(7344:end);
    dates3=dates3(7344:end);
    dates4=dates4(7344:end);
   
    %% SÃ©lection des dates communes + processing eddyflux
    obs_FB.TIMESTAMP=datenum(obs_FB.TIMESTAMP);
    
    [ind_t,~,ind_t4]=deal (  find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end)),find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end))  , find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end))   );
    % Variable Ã  12h pour albedo notamment
    ind_t_12=find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
    ind_t_123=find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
    ind_t_124=find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
   
    % Intersection pour les scatterplot
    [~,i_obs,i_sim]=intersect(round(obs_FB.TIMESTAMP*48),round(dates4*48));
    
    % retirer les premiers mois de la simulation (erreurs d'initialisation)
    % variables cumulatives comme ETR cumulé (ETR TOT (23))
% %     % on enlève aout au 1er janvier
    temp2=temp2(7344:end,1:end);
    temp3=temp3(7344:end,1:end);
    temp4=temp4(7344:end,1:end);
    %i_obs=i_obs1(8600:end);
    %i_sim=i_sim1(8600:end);
    %dates_b=dates(8600:end);

    % processing H et LE
    %hard filter (range filter+rainfilter)
    % on rajoute + spike flag activé + Test SS_ITC 
    for i=1:height(obs_FB)
        if obs_FB.LE_SPIKE_FLAG(i)==1
            obs_FB.LE_hf(i)=NaN;
        end 
        if obs_FB.H_SPIKE_FLAG(i)==1
            obs_FB.H_hf(i)=NaN;
        end 
         if obs_FB.LE_SSITC_TEST(i)==2 % 2=mauvais 
            obs_FB.LE_hf(i)=NaN;
         end 
          if obs_FB.H_SSITC_TEST(i)==2
            obs_FB.H_hf(i)=NaN;
        end 

    end 




    % Filtrage des valeurs aberrantes
    obs_FB.SWC_TD_5_1(obs_FB.SWC_TD_5_1<0 | obs_FB.SWC_TD_5_1>0.5)=NaN;     % ????????
    obs_FB.ALB(obs_FB.ALB<0 | obs_FB.ALB>1)=NaN;
    %obs_FB.H_hf(obs_FB.H<-50)=NaN;
    %obs_FB.LE_hf(obs_FB.LE_hf<-50)=NaN;
    
    % CORRECTION DU BILAN D'ENERGIE
    % Calcul du déphasage et du rappport d'amplitude à 5cm partir de la simulation
    G0j_sim=zeros(48,1);% Composite journalier du G simulé sur l'année
    G5j_sim=zeros(48,1);
    for j=1:48
        G0j_sim(j)=mean(temp2(j:48:end,15),1);
        G5j_sim(j)=mean(temp2(j:48:end,67),1);
       
    end
    [A0,A5]= deal ( max(G0j_sim),max(G5j_sim)) ;
    [ia0,ia5]=deal(  find(G0j_sim==A0),find(G5j_sim==A5));
     G0_obs=(G_obs(ind_t)+ia5-ia0)*A0/A5;
    %G0_obs=G_obs(i_obs+ia5-ia0)*A0/A5;
     
    % Contrôle graphique
%         figure
%         plot(obs_FB.TIMESTAMP(ind_t),G0_obs)
%         hold on
%         plot(obs_FB.TIMESTAMP(ind_t),G_obs(ind_t))
%         datetick('x')
%         legend('surface','5cm')
  
    % Calcul de la fermeture du bilan
    H_obs=obs_FB.H_hf(ind_t);
    LE_obs=obs_FB.LE_hf(ind_t);
    Rn_obs=obs_FB.NETRAD(ind_t);
    Residu=Rn_obs-G0_obs-H_obs-LE_obs;
    Residu(abs(Residu)>250)=NaN; % filtrage des résidus abhérants
    
    % Bilan annuel
    bilan=nanmean(Residu(1:8400));
    bilan_std=nanstd(Residu(1:8400));
    
    % Attribution des résidus à LE et H selon leurs amplitudes relatives
    S=abs(H_obs)+abs(LE_obs);
    prc_H=abs(H_obs)./S;
    prc_LE=abs(LE_obs)./S;
    H_corr=H_obs+(Residu.*prc_H);
    LE_corr=LE_obs+(Residu.*prc_LE);

    bis_LE=LE_corr-LE_obs;
    
       %% CORRECTION PRISE EN COMPTE PROPORTION CAILLOUX
    %les teneur en eau mesurées sont pour un sol 100% fin et les teneurs
    %simulés sont pour un sol avec des cailloux
    % pour pouvoir comparer rescaler les teneurs en eau simulés  
    
    % Pourcentage de cailloux FB 50,60,65,70,80,85,90,90%
    % humidités à 
    p=[2,5,15,25,35,50,75,100,150,200];
    col=[45,46,47,48,49,50,51,52,53]; % colonnes de temp correspondantes
    cailloux=[0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1,0.1]; % pourcentage de sol
     c= [0.85,0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1];  
    for s=1:9
    %temp(:,col(s))=temp(:,col(s))/cailloux(s);
    temp2(:,col(s))=temp2(:,col(s))/cailloux(s);
    temp3(:,col(s))=temp3(:,col(s))/cailloux(s);
    temp4(:,col(s))=temp4(:,col(s))/cailloux(s);
    end 


%% sapflow import data
sapflow=importdata("/home/brune/SISPAT/Stage 2023/Données d'entrée/flux de sève/sapflow_2015_2021.mat");
sapflowj=importdata("/home/brune/SISPAT/Stage 2023/Données d'entrée/flux de sève/sapflow_jour_2015_2021.mat");
sapflowj=table2timetable(sapflowj);
%% sapflow jour
figure_handle=figure('Name',['ET sapflow jour' ])
subplot(4,1,1)
plot(datenum(sapflowj.TIMESTAMP),sapflowj.SAP_FLOW_TD_P)
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on
ylabel('Sapflow (dm3/j)')
%ylim([0 0.15])
legend('sapflow PIN')

subplot(4,1,2)
plot(datenum(sapflowj.TIMESTAMP),sapflowj.SAP_FLOW_TD_C)
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on
ylabel('Sapflow (dm3)')
legend('sapflow Chene')

TR_sim=timeavgCK(dates,temp2(:,13),'day',@nansum)*dt/Lv;
subplot(4,1,3)
plot(datenum(sapflowj.TIMESTAMP),TR_sim)
ylabel('Transpi (mm)')
grid on 
legend('Transpi sim')
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')

           subplot(4,1,4)
plot(datenum(sapflowj.TIMESTAMP), ETR_day_obs,'Color','b')
ylabel('ET (mm)')
legend('ET obs')
grid on 
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')

all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

%% sapflow 30 min
figure_handle=figure('Name',['ET sapflow' ])
subplot(4,1,1)
plot(datenum(sapflow.TIMESTAMP),sapflow.SAP_FLOW_TD_P)
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on
ylabel('Sapflow (dm3)')
ylim([0 0.15])
legend('sapflow PIN')
subplot(4,1,2)
plot(datenum(sapflow.TIMESTAMP),sapflow.SAP_FLOW_TD_C)
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on
ylabel('Sapflow (dm3)')
legend('sapflow Chene')

subplot(4,1,3)
plot(dates,temp3(:,13)*dt/Lv)
ylabel('Transpi (mm)')
grid on 
legend('Transpi sim')
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')

           subplot(4,1,4)
plot(dates, obs_FBt.LE_hf(ind_t)*dt/Lv,'Color','b')
ylabel('ET (mm)')
legend('ET obs')
grid on 
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')

all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

%% REPARTITION PAR PERIODE DES COMPOSANTES DE L'ET
TTobs = retime(obs_FBt(ind_t,:),'regular','sum','TimeStep',calmonths(1));
[ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),'day',@nansum);
ETR_day_obs=ETR_day_obs*dt/Lv;
TTobsET=array2table([days_obs ETR_day_obs]);
TTobsET.Var1=datetime(TTobsET.Var1,'ConvertFrom','datenum');
TTobsET= table2timetable(TTobsET);
TTobsET=retime(TTobsET,'regular','sum','TimeStep',calmonths(1));
%TTobsET= groupsummary(TTobsET, 'Var1', 'monthofyear', 'mean');

   ind_LEs=isnan(obs_FB.LE_sf(ind_t));
   ind_LEh=isnan(obs_FB.LE_hf(ind_t));
   ind_LE=isnan(obs_FB.LE(ind_t));
   
TTsim1=1;
TTsim2=1;
TTsim3=1;
ETobs_sim=[];
for i=1:3
       temp_=eval(['temp',num2str(i+1)]);
   
    LEf(:,i)=temp_(:,12)*dt/Lv;% evap feuille (eau interceptée)
    LEtr(:,i)=temp_(:,13)*dt/Lv;%transpi
    LEs(:,i)=temp_(:,11)*dt/Lv;% evap sol
    LEtot(:,i)=temp_(:,14)*dt/Lv;% total

      for j=1:length(obs_FB.LE_hf(ind_t))
        if ind_LEh(j) ==1
             LEf(j,i)=NaN;
             LEs(j,i)=NaN;
             LEtr(j,i)=NaN;
             LEtot(j,i)=NaN;
        end 
      end 
% TTt=timeavgCK(dates,LEtot,'day',@nansum);
% TTtr=timeavgCK(dates,LEtr,'day',@nansum);
% TTes=timeavgCK(dates,LEs,'day',@nansum);
% TTef=timeavgCK(dates,LEf,'day',@nansum);
TTsim=array2table([dates LEtr(:,i) LEs(:,i) LEf(:,i) LEtot(:,i) ind_LE ind_LEs ind_LEh]);
TTsim.Var1=datetime(TTsim.Var1,'ConvertFrom','datenum');
TTsim= table2timetable(TTsim);
TTsimj=retime(TTsim,'regular','sum','TimeStep',days(1));
TTsim=retime(TTsim,'regular','sum','TimeStep',calmonths(1));
TTsim_moy= groupsummary(TTsim, 'Var1', 'monthofyear', 'mean');
TTsim_min= groupsummary(TTsim, 'Var1', 'monthofyear', 'min');
TTsim_max= groupsummary(TTsim, 'Var1', 'monthofyear', 'max');
%TTsim_Q25= groupsummary(TTsim, 'Var1', 'monthofyear', Q(TTsim.Var2,0.25));
if i==1
    TTsim1=TTsim;
elseif i==2
    TTsim2=TTsim;
elseif i==3
    TTsim3=TTsim;
end 

figure_handle=figure('Name',['ET' num2str(i) ]);
  subplot(3,1,1)
% x=categorical({'JFM' 'AMJ' 'JAS' 'OND' })
% x = reordercats(x,{'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' });
 bar(TTobs.TIMESTAMP, TTobs.P)
 legend('1 months total precipitation')
ylabel('P (mm)')

grid on
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

subplot(3,1,2)
for ii=1:length(TTobs.TIMESTAMP)
    %y= TTsim.Var5(i);
    y(ii,:) = [TTsim.Var2(ii) TTsim.Var3(ii) TTsim.Var4(ii)];
    hold on 
end

bar(TTobs.TIMESTAMP,y,'stacked')
hold on
scatter(TTobs.TIMESTAMP,TTobsET.Var2,50,'k','filled')
legend('Transpi. sim','Soil evap. sim','Leaf evap. sim','ET obs')
ylabel('ET (mm)')
  dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on 
subplot(3,1,3)           
scatter(TTsim.Var1,TTsim.Var8)
legend('missing data LE hf')
grid on
           all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
% ETobs_sim=TTsim.Var1;
% timetable(TTsim.Var1,)

ETobs_sim(:,i)=TTsim.Var5-TTobsET.Var2;

end
ETobs_sim=array2table(ETobs_sim);

ETobs_sim = table2timetable(ETobs_sim, 'RowTimes', TTsim.Var1');
%ETobs_sim=[TTsim.Var1 ETobs_sim];
ET_sim_moy=groupsummary(ETobs_sim, 'Time', 'monthofyear', 'mean');
ET_sim_min= groupsummary(ETobs_sim, 'Time','monthofyear', 'min');
ET_sim_max= groupsummary(ETobs_sim, 'Time', 'monthofyear', 'max');
%ET_sim_max=table2array(ET_sim_max);
%ET_sim_min=table2array(ET_sim_min);
%% 
figure_handle=figure('Name',['ET' num2str(i) ])
  subplot(4,1,1)
% x=categorical({'JFM' 'AMJ' 'JAS' 'OND' })
% x = reordercats(x,{'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' 'JFM' 'AMJ' 'JAS' 'OND' });
 bar(TTobs.TIMESTAMP, TTobs.P)
 legend('1 months total precipitation')
ylabel('P (mm)')
grid on
 dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
     ax = gca;
ax.FontSize = 10; datetime
ax.FontWeight = 'bold'; 

subplot(4,1,2)
for ii=1:length(TTobs.TIMESTAMP)
    %y= TTsim.Var5(i);
    y(ii,:) = [TTsim1.Var2(ii) TTsim1.Var3(ii) TTsim1.Var4(ii)];
    hold on 
end

bar(TTobs.TIMESTAMP,y,'stacked')
hold on
scatter(TTobs.TIMESTAMP,TTobsET.Var2,50,'k','filled')
legend('Transpi. sim','Soil evap. sim','Leaf evap. sim','ET obs')
ylabel('ET (mm)')
  dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on 

subplot(4,1,3)
for ii=1:length(TTobs.TIMESTAMP)
    %y= TTsim.Var5(i);
    y(ii,:) = [TTsim2.Var2(ii) TTsim2.Var3(ii) TTsim2.Var4(ii)];
    hold on 
end

bar(TTobs.TIMESTAMP,y,'stacked')
hold on
scatter(TTobs.TIMESTAMP,TTobsET.Var2,50,'k','filled')
legend('Transpi. sim','Soil evap. sim','Leaf evap. sim','ET obs')
ylabel('ET (mm)')
  dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on 

subplot(4,1,4)
for ii=1:length(TTobs.TIMESTAMP)
    %y= TTsim.Var5(i);
    y(ii,:) = [TTsim3.Var2(ii) TTsim3.Var3(ii) TTsim3.Var4(ii)];
    hold on 
end

bar(TTobs.TIMESTAMP,y,'stacked')
hold on
scatter(TTobs.TIMESTAMP,TTobsET.Var2,50,'k','filled')
legend('Transpi. sim','Soil evap. sim','Leaf evap. sim','ET obs')
ylabel('ET (mm)')
  dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on 

           all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


% calcul de excédent ou deficit d'ET mensuel simulé par rapport aux obs
%ETobs_sim=array2table(ETobs_sim);
%ETobs_sim=[TTsim.Var1 ETobs_sim];
figure()
yy=[ETobs_sim.ETobs_sim1  ETobs_sim.ETobs_sim2 ETobs_sim.ETobs_sim3];
bar(TTobs.TIMESTAMP,yy)
grid on
ylabel('Surplus ou déficit ET simulé (mm)')
legend('sim. standard', 'sim. karst ')
dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
           grid on 

           all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


figure()
D = datetime(2021,1,1):calmonths(1):datetime(2021,12,31);
D.Format = 'MMM';
yyy=[ET_sim_moy(:,"mean_ETobs_sim1")  ET_sim_moy(:,"mean_ETobs_sim3")];
bar(string(D),table2array(yyy))
hold on 
ngroups = size(yyy, 1);
nbars = size(yyy, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    %x=x';
    %errorbar(x, yy(:,i), err(:,i), '.');
    %errorbar(x,yy(:,i),(table2array(ET_sim_max(:,i+2))+table2array(ET_sim_min(:,i+2)))/2,(table2array(ET_sim_max(:,i+2))-table2array(ET_sim_min(:,i+2)))/2,'.','Color','k');
    errorbar(x,table2array(yyy(:,i)),table2array(yyy(:,i))-table2array(ET_sim_min(:,i+2)),table2array(ET_sim_max(:,i+2))-table2array(yyy(:,i)),'.','Color','k');
end
hold off
grid on
ylabel('Surplus ou déficit ET simulé (mm)')
legend('sim. standard', 'sim. karst ','min et max')
 ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

%mettre dans boucle 
scoresLE1s=vr_scoresstr(obs_FB.LE_sf(i_obs),temp2(i_sim,14));
scoresLE2s=vr_scoresstr(obs_FB.LE_sf(i_obs),temp3(i_sim,14));
scoresLE3s=vr_scoresstr(obs_FB.LE_sf(i_obs),temp4(i_sim,14));
scoresLE1sc=vr_scoresstr(LE_corr,temp2(i_sim,14));
scoresLE1h=vr_scoresstr(obs_FB.LE_hf(i_obs),temp2(i_sim,14));
scoresLE2h=vr_scoresstr(obs_FB.LE_hf(i_obs),temp3(i_sim,14));
scoresLE3h=vr_scoresstr(obs_FB.LE_hf(i_obs),temp4(i_sim,14));
scoresLE1hc=vr_scoresstr(LE_corr,temp2(i_sim,14));
s_LE_1=[scoresLE1s.rmse scoresLE1s.mbias scoresLE1s.cor2 sum(ind_LEs);
    scoresLE1sc.rmse scoresLE1sc.mbias scoresLE1sc.cor2 sum(ind_LEs);
    scoresLE1h.rmse scoresLE1h.mbias scoresLE1h.cor2 sum(ind_LEh);
    scoresLE1hc.rmse scoresLE1hc.mbias scoresLE1hc.cor2 sum(ind_LEh)];

% nombre de valeurs manquantes pour LE obs par mois 
figure()
scatter(TTsim.Var1,TTsim.Var7)
hold on
scatter(TTsim.Var1,TTsim.Var8)
legend('LE sf', 'LE hf')
ylabel('missing data')
grid on