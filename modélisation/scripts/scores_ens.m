%% calcul de scores des ensembles
clear;
close all;
clc;
addpath('/home/brune/SISPAT/Fonctions matlab/')

%quelle année considérée?
annee=2014;
% annee=2017;
Lv=2.484*10^6;
dt=1800;
% afficher figures oui=1   non=0
fig_ET=0;

%% chargement des mesures obs_FB

obs_FBt=importdata(['/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat']);
obs_FBt=table2timetable(obs_FBt);
obs_FB=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');

G_obs=obs_FB.G_5;
Pcumuljour=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/Pcumuljour.mat');
%% Lecture du (des) fichier(s) atm
numsim=2000;
for isim = 1:numsim
    numsim=sprintf('%04d',isim);
        %numsim=3;
    % path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/fev/roote2/';
    % files_atm2=dir([path_res2 'atm_FB_14_21_20m.out']);

    path_out='/media/brune/Ultra Touch/SISPAT/out6';
   % path_out='/home/brune/SISPAT/execution/sortie_modele/FB/';
    %path_sol=strcat(path_out,sprintf('/sol/sol__%s.txt',numsim));
    path_atm=sprintf('/atm/atm_%s.out',numsim);
    path_atm=[path_out sprintf('/atm/atm_%s.out',numsim)];
% simulation de référence 
    %path_ref=;
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
% file_atm2=files_atm2.name;
%     fid = fopen([path_res2 file_atm2],'r');
%     temp2 = fread(fid,'single');
%     nbline= length(temp2)/86; % 86 =nombre de variables par ligne
%     temp = reshape(temp2,[86,nbline])';

 %file_a=path_atm.name;   
    fid = fopen(path_atm,'r');
     temp = fread(fid,'single');
    nbline= length(temp)/86; % 86 =nombre de variables par ligne
    temp = reshape(temp,[86,nbline])';
    fclose(fid);
    dates=datenum(annee,8,1)-1+temp(:,1)+(1/24)*temp(:,2)+(1/(24*60))*temp(:,3);
    dates=dates(7344:end);  %tronquer les premiers mois de simulation (1er aout au 1er janvier)    
   
    %% SÃ©lection des dates communes + processing eddyflux
    obs_FB.TIMESTAMP=datenum(obs_FB.TIMESTAMP);
        ind_t=find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end));
    % Variable Ã  12h pour albedo notamment
    ind_t_12=find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
      % Intersection pour les scatterplot
    [~,i_obs,i_sim]=intersect(round(obs_FB.TIMESTAMP*48),round(dates*48));
        % retirer les premiers mois de la simulation (erreurs d'initialisation)
    % variables cumulatives comme ETR cumulé (ETR TOT (23))
% %     % on enlève aout au 1er janvier
    temp=temp(7344:end,1:end);
   
    % processing H et LE
    %hard filter (range filter+rainfilter) % on rajoute + spike flag activé + Test SS_ITC 
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
        G0j_sim(j)=mean(temp(j:48:end,15),1);
        G5j_sim(j)=mean(temp(j:48:end,67),1);
    end
    [A0,A5]= deal ( max(G0j_sim),max(G5j_sim)) ;
    [ia0,ia5]=deal(  find(G0j_sim==A0),find(G5j_sim==A5));
     G0_obs=(G_obs(ind_t)+ia5-ia0)*A0/A5;
    %G0_obs=G_obs(i_obs+ia5-ia0)*A0/A5;
      
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
    p=[2,5,15,25,35,50,75,100,150,200]; %profondeur
    col=[45,46,47,48,49,50,51,52,53]; % colonnes de temp correspondantes
    cailloux=[0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1,0.1]; % pourcentage de sol
     c= [0.85,0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1];  
    for s=1:9
    temp(:,col(s))=temp(:,col(s))/cailloux(s);
    end

    %% 
   ind_LE=isnan(obs_FB.LE_hf(ind_t));
   LEbis=temp(:,14);
      for i=1:length(obs_FB.LE_hf(ind_t))
        if ind_LE(i) ==1
             LEbis(i)=NaN;
        end 
    end 
    
      % ETR journalière
    %N_LE nombre de valeurs qu'il y a. On peut filtrer les jours avec trop
    %de valeurs manquantes si plus de 4 h manquantes par exemple nan  pour
    %la journée pour les obs et les simu aux jours correspondant 
    [ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),'day',@nansum);
    ETR_day_obs=ETR_day_obs*dt/Lv;
    [ETR_day_sim,days_sim]=timeavgCK(dates,temp(:,14),'day',@nansum);
    ETR_day_sim=ETR_day_sim*dt/Lv;
    for i=1:length(ETR_day_obs)
        if N_LE_obs(i) <24 % si inférieur à 40 mesures par jour (20 h/24) pas du cumul journalier
            ETR_day_obs(i)=NaN;          
        end 
    end 
    %% Calcul de scores obs / simu
    scores.('RN').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.NETRAD(i_obs),temp(i_sim,7));
    scores.('LEb').(['sim_' num2str(isim)])=vr_scoresstr(LE_corr,temp(i_sim,14));
    scores.('LE').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.LE_hf(i_obs),temp(i_sim,14));
    scores.('Hb').(['sim_' num2str(isim)])=vr_scoresstr(H_corr,temp(i_sim,10));
    scores.('ET_j').(['sim_' num2str(isim)])=vr_scoresstr(ETR_day_obs,ETR_day_sim);
    scores.('H').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.H_hf(i_obs),temp(i_sim,10));
    scores.('G').(['sim_' num2str(isim)])=vr_scoresstr(G_obs(i_obs),temp(i_sim,67));
    scores.('G0').(['sim_' num2str(isim)])=vr_scoresstr(G0_obs,temp(i_sim,15));
    scores.('T5').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.TS_TD_5_1(i_obs),temp(i_sim,36));
    scores.('T15').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.TS_TD_15_1(i_obs),temp(i_sim,37));
    scores.('T25').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.TS_TD_25_1(i_obs),temp(i_sim,38));
    scores.('T35').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.TS_TD_35_1(i_obs),temp(i_sim,39));
    scores.('T50').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.TS_TD_50_1(i_obs),temp(i_sim,40));
    scores.('W5').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.SWC_TD_5_1(i_obs),temp(i_sim,45)/100);
    scores.('W15').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.SWC_TD_15_1(i_obs),temp(i_sim,46)/100);
    scores.('W25').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.SWC_TD_25_2(i_obs),temp(i_sim,47)/100);
    scores.('W35').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.SWC_TD_35_2(i_obs),temp(i_sim,48)/100);
    scores.('W50').(['sim_' num2str(isim)])=vr_scoresstr(obs_FB.SWC_TD_50_2(i_obs)-0.2,temp(i_sim,49)/100);
file_score=['scores_' num2str(isim) '.mat']
%scores_ensemble1=scores;
  %  save([path_out '/' 'scores_ensemble2.mat'])
  
%% Calcul de scores ref / simu
  %scoresSWout=vr_scoresstr((obs_FB.SW_IN(i_obs)).*(obs_FB.ALB(i_obs)),(obs_FB.SW_IN(i_obs)).*(temp(i_sim,65)));
    % ref_scores.('RN').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,7),temp(i_sim,7));
    % ref_scores.('LE').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,14),temp(i_sim,14));
    % %scoresH2.(['sim_' num2str(isim)])=vr_scoresstr(H_corr,temp(i_sim,10));
    % ref_scores.('H').([sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,10),temp(i_sim,10));
    % ref_scores.('G').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,67),temp(i_sim,67));
    % ref_scores.('G0').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,15)),temp(i_sim,15));
    % ref_scores.('T5').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,36)),temp(i_sim,36));
    % ref_scores.('T15').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,37),temp(i_sim,37));
    % ref_scores.('T25').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,38),temp(i_sim,38));
    % ref_scores.('T35').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,39),temp(i_sim,39));
    % ref_scores.('T50').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,40),temp(i_sim,40));
    % ref_scores.('W5').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,45)/100,temp(i_sim,45)/100);
    % ref_scores.('W15').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,46)/100,temp(i_sim,46)/100);
    % ref_scores.('W25').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,47)/100,temp(i_sim,47)/100);
    % ref_scores.('W35').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,48)/100,temp(i_sim,48)/100);
    % ref_scores.('W50').(['sim_' num2str(isim)])=vr_scoresstr(ref_temp(i_sim,49)/100,temp(i_sim,49)/100);

%% REPARTITION PAR MOIS DES COMPOSANTES DE L'ET
TTobs = retime(obs_FBt(ind_t,:),'regular','sum','TimeStep',calmonths(1));
[ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),LE_corr,'day',@nansum);
ETR_day_obs=ETR_day_obs*dt/Lv;
TTobsET=array2table([days_obs ETR_day_obs]);
TTobsET.Var1=datetime(TTobsET.Var1,'ConvertFrom','datenum');
TTobsET= table2timetable(TTobsET);
TTobsET=retime(TTobsET,'regular','sum','TimeStep',calmonths(1));
%TTobsET= groupsummary(TTobsET, 'Var1', 'monthofyear', 'mean');
 ind_LEs=isnan(obs_FB.LE_sf(ind_t));
 ind_LEh=isnan(obs_FB.LE_hf(ind_t));
 ind_LE=isnan(obs_FB.LE(ind_t));
   
TTsim=1;
i=1;
ETobs_sim=[];
%temp=eval(['temp',num2str(i+1)]);
    LEf=temp(:,12)*dt/Lv;% evap feuille (eau interceptée)
    LEtr=temp(:,13)*dt/Lv;%transpi
    LEs=temp(:,11)*dt/Lv;% evap sol
    LEtot=temp(:,14)*dt/Lv;% total

      for j=1:length(obs_FB.LE_hf(ind_t))
        if ind_LEh(j) ==1
             LEf(j,1)=NaN;
             LEs(j,1)=NaN;
             LEtr(j,1)=NaN;
             LEtot(j,1)=NaN;
        end 
      end

TTsim=array2table([dates LEtr(:,i) LEs(:,i) LEf(:,i) LEtot(:,i) ind_LE ind_LEs ind_LEh]);
TTsim.Var1=datetime(TTsim.Var1,'ConvertFrom','datenum');
TTsim= table2timetable(TTsim);
TTsimj=retime(TTsim,'regular','sum','TimeStep',days(1));
TTsim=retime(TTsim,'regular','sum','TimeStep',calmonths(1));
TTsim_moy= groupsummary(TTsim, 'Var1', 'monthofyear', 'mean');
TTsim_min= groupsummary(TTsim, 'Var1', 'monthofyear', 'min');
TTsim_max= groupsummary(TTsim, 'Var1', 'monthofyear', 'max');
%TTsim_Q25= groupsummary(TTsim, 'Var1', 'monthofyear', Q(TTsim.Var2,0.25));
if fig_ET==1     
figure_handle=figure('Name',['ET' num2str(i) ]);
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
% linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
end

ETobs_sim=TTsim.Var5-TTobsET.Var2; % écart ET mensuel mesuré et simulé
ETobs_sim=array2table(ETobs_sim);
ETobs_sim = table2timetable(ETobs_sim, 'RowTimes', TTsim.Var1');
ET_scores.('obs_sim_mensuel').(['sim_' num2str(isim)])=ETobs_sim;
%ETobs_sim=[TTsim.Var1 ETobs_sim];
% moyenne, max,min des mois de janvier, des mois de février... 
ET_sim_moy=groupsummary(ETobs_sim, 'Time', 'monthofyear', 'mean');
ET_sim_min= groupsummary(ETobs_sim, 'Time','monthofyear', 'min');
ET_sim_max= groupsummary(ETobs_sim, 'Time', 'monthofyear', 'max');

ET_sim_mens= horzcat(double(string(ET_sim_max.monthofyear_Time)), ET_sim_moy.mean_ETobs_sim , ET_sim_min.min_ETobs_sim, ET_sim_max.max_ETobs_sim);
ET_sim_mens=array2table(ET_sim_mens);
ET_sim_mens=renamevars(ET_sim_mens,["ET_sim_mens2","ET_sim_mens3","ET_sim_mens4"],["mean","min","max"]);
ET_scores.('mois_obs_sim_mensuel').(['sim_' num2str(isim)])=ET_sim_mens;







%% Autres données drainage, stock, 

scores.('drainage').(['sim_' num2str(isim)])=temp(end,28)-temp(1,28);


end

%%
save([path_out '/' 'scores_ensemble6_1.mat'],'scores')
%clear('scores')
save([path_out '/' 'scores_ET6_1.mat'],'ET_scores')
