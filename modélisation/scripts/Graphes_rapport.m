% Creation graphes du rapport 
clear;
close all;
clc;

 
addpath('/home/brune/SISPAT/Fonctions matlab')

%annee=2017;
annee=2017;
Lv=2.484*10^6;
dt=1800;
%% 


%% Chargement des obs

obs_FBt=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');
obs_FBt=table2timetable(obs_FBt);
obs_FB=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/tableau_obs_FB_14_21.mat');

G_obs=obs_FB.G_5;
Pcumuljour=importdata('/home/brune/SISPAT/tableaux matlab/Font-Blanche/Pcumuljour.mat');

%% Lecture du (des) fichier(s) atm
path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/K1/';
path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/K1d2/';
path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6/';
%[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_17_21_20m_28.out']),dir([path_res3 'atm_FB_17_21_20m_28.out']),dir([path_res4 'atm_FB_17_21_20m_28.out']));
[files_atm2,files_atm3,files_atm4] = deal(dir([path_res2 'atm_FB_14_21_20m.out']),dir([path_res3 'atm_FB_14_21_20m.out']),dir([path_res4 'atm_FB_14_21_20m.out']));


for i=1:length(files_atm4)
    
    file_atm2=files_atm2(i).name;
    fid = fopen([path_res2 file_atm2],'r');
    temp2 = fread(fid,'single');
    nbline= length(temp2)/86; % 86 =nombre de variables par ligne
    temp2 = reshape(temp2,[86,nbline])';
    %temp=temp4;
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

    % Ordre des variables
    % (01)JOUR      (02)HEURE     (03)MINUTE    (04)HEURE     (05)RN SOL    (06)RN VEG
    % (07)RN TOT   F (08)H SOL     (09)H VEG     (10)H TOT     (11)LE SOL    (12)LE FEUIL
    % (13)LE VEG    (14)LE TOT    (15)G         (16)ETR/ETP   (17)TF        (18)HF
    % (19)Trad      (20)TS        (21)TAF       (22)QAF       (23)ETR TOT    (24)TR
    % (25)E SOL     (26)E FEUIL   (27)RUIS SUR  (28)PERCOL    (29)PLUIE     (30)STOCK
    % (31)Rb        (32)Ras       (33)Rsto      (34)ustar     (35)Uaf       (36)TS 5
    % (37)TS 15cm   (38)TS 25cm   (39)TS 35cm  (40)TS 50cm  (41)TS 200cm  (42)TS 500cm
    % (43)TS 1000   (44)TS 2000   (45)W 5    (46)W 15cm    (47)W 25cm    (48)W 35cm
    % (49)W 50cm   (50)W 200cm   (51)W 500cm   (52)W 1000    (53)W 2000    (54)WI 0-20cm
    % (55)WI 20-30  (56)WI 30-40 (57)WI 40-60 (58)WI 60-100 (59)WI 100-200(60)WI200-500   35
    %(61)WI total  (62)Stock     (63)Q surf    (64)Albsol    (65)Albcouv   (66)Ra
    % (67)G 5cm     (68)ETP       (69)Flxliqsurf              (70)Flxliq10cm
    % (71)Flxliq50cm              (72)Flxliq100cm             (73) Flxliq150cm
    % (74)Flxliq200cm             (75)Flxliq250cm             (76)Flxliqrien
    % (77)Flxliqrien              (78)Flxvapsurf              (79)Flxvap10cm
    % (80)Flxvap50cm              (81)Flxvap100cm             (82) Flxvap150cm
    % (83)Flxvap200cm             (84)Flxvap250cm             (85)Flxvaprien
    % (86)Flxvaprien
end 
    dates=datenum(annee,8,1)-1+temp2(:,1)+(1/24)*temp2(:,2)+(1/(24*60))*temp2(:,3);
    dates3=datenum(annee,8,1)-1+temp3(:,1)+(1/24)*temp3(:,2)+(1/(24*60))*temp3(:,3);
    dates4=datenum(annee,8,1)-1+temp4(:,1)+(1/24)*temp4(:,2)+(1/(24*60))*temp4(:,3);
%tronquer les premiers mois de simulation (6MOIS°    
    dates=dates(7344:70127);
    dates3=dates3(7344:70127);
    dates4=dates4(7344:70127);
   
    %% SÃ©lection des dates communes
    obs_FB.TIMESTAMP=datenum(obs_FB.TIMESTAMP);
    
    [ind_t,ind_t3,ind_t4]=deal (  find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end)),find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end))  , find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end))   );
    % Variable Ã  12h pour albedo notamment
    ind_t_12=find(obs_FB.TIMESTAMP>=dates(1) & obs_FB.TIMESTAMP<=dates(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
    ind_t_123=find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
    ind_t_124=find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end) & abs(hour(datetime(obs_FB.TIMESTAMP,'ConvertFrom','datenum'))-12)<(5/60));
   
    % Intersection pour les scatterplot
    [~,i_obs,i_sim]=intersect(round(obs_FB.TIMESTAMP*48),round(dates4*48));
    
    % retirer les premiers mois de la simulation (erreurs d'initialisation)
    % variables cumulatives comme ETR cumulé (ETR TOT (23))
    % on enlève aout au 1er janvier
    temp2=temp2(7344:70127,1:end);
    temp3=temp3(7344:70127,1:end);
    temp4=temp4(7344:70127,1:end);

   
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


    %i_obs=i_obs1(8600:end);
    %i_sim=i_sim1(8600:end);
    %dates_b=dates(8600:end);
    % Filtrage des valeurs aberrantes
    obs_FB.SWC_TD_5_1(obs_FB.SWC_TD_5_1<0 | obs_FB.SWC_TD_5_1>0.5)=NaN;     % ????????
    obs_FB.ALB(obs_FB.ALB<0 | obs_FB.ALB>1)=NaN;
    obs_FB.H_hf(obs_FB.H<-50)=NaN;
    obs_FB.LE_hf(obs_FB.LE_hf<-50)=NaN;
    
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
     G0_obs=G_obs(ind_t+ia5-ia0)*A0/A5;
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

    
       %% CORRECTION PRISE EN COMPTE PROPORTION CAILLOUX
    %les teneur en eau mesurées sont pour un sol 100% fin et les teneurs
    %simulés sont pour un sol avec des cailloux
    % pour pouvoir comparer rescaler les teneurs en eau simulés  
    
    % Pourcentage de cailloux FB 50,60,65,70,80,85,90,90%
    % humidités à 
        col=[45,46,47,48,49,50,51,52,53]; % colonnes de temp correspondantes
    cailloux=[0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1,0.1]; % pourcentage de sol
%     c= [0.85,0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1];  
%     bd=[0.55 0.32 0.5 0.5 0.48 0.45 0.45 0.45 0.45];
%     bres=[0.08 0.1 0.12 0.12 0.12 0.12 0.12 0.12 0.12];
    for s=1:9
    %temp(:,col(s))=temp(:,col(s))/cailloux(s);
    temp2(:,col(s))=temp2(:,col(s))/cailloux(s);
    temp3(:,col(s))=temp3(:,col(s))/cailloux(s);
    temp4(:,col(s))=temp4(:,col(s))/cailloux(s);
    end 
    %% Calcul de scores
    %scoresSWout=vr_scoresstr((obs_FB.SW_IN(i_obs)).*(obs_FB.ALB(i_obs)),(obs_FB.SW_IN(i_obs)).*(temp(i_sim,65)));
    scoresRN1=vr_scoresstr(obs_FB.NETRAD(i_obs),temp2(i_sim,7));
    scoresRN2=vr_scoresstr(obs_FB.NETRAD(i_obs),temp3(i_sim,7));
    scoresRN3=vr_scoresstr(obs_FB.NETRAD(i_obs),temp4(i_sim,7));
    scoresLE1=vr_scoresstr(obs_FB.LE_hf(i_obs),temp2(i_sim,14));
    scoresLE2=vr_scoresstr(obs_FB.LE_hf(i_obs),temp3(i_sim,14));
    scoresLE3=vr_scoresstr(obs_FB.LE_hf(i_obs),temp4(i_sim,14));
    
    scoresH1=vr_scoresstr(obs_FB.H_hf(i_obs),temp2(i_sim,10));
    scoresH2=vr_scoresstr(obs_FB.H_hf(i_obs),temp3(i_sim,10));
    scoresH3=vr_scoresstr(obs_FB.H_hf(i_obs),temp4(i_sim,10));
    scoresG1=vr_scoresstr(G_obs(i_obs),temp2(i_sim,67));
    scoresG2=vr_scoresstr(G_obs(i_obs),temp3(i_sim,67));
    scoresG3=vr_scoresstr(G_obs(i_obs),temp4(i_sim,67));
    scoresG01=vr_scoresstr(G0_obs,temp2(i_sim,15));
    scoresG02=vr_scoresstr(G0_obs,temp3(i_sim,15));
    scoresG03=vr_scoresstr(G0_obs,temp4(i_sim,15));
       for i=1:3
       temp_=eval(['temp',num2str(i+1)]);
    scoresT5.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.TS_TD_5_1(i_obs),temp_(i_sim,36));
    scoresT15.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.TS_TD_15_1(i_obs),temp_(i_sim,37));
    scoresT25.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.TS_TD_25_1(i_obs),temp_(i_sim,38));
    scoresT35.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.TS_TD_35_1(i_obs),temp_(i_sim,39));
    scoresT50.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.TS_TD_50_1(i_obs),temp_(i_sim,40));
           scoresW5.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.SWC_TD_5_1(i_obs),temp_(i_sim,45)/100);
    scoresW15.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.SWC_TD_15_1(i_obs),temp_(i_sim,46)/100);
    scoresW25.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.SWC_TD_25_2(i_obs),temp_(i_sim,47)/100);
    scoresW35.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.SWC_TD_35_2(i_obs),temp_(i_sim,48)/100);
    scoresW50.(['sim_' num2str(i)])=vr_scoresstr(obs_FB.SWC_TD_50_2(i_obs)-0.2,temp_(i_sim,49)/100);
     end 

%% GRAPHES 


%% Scatter plot bilan d'energie (avec saisons)

% filtrer valeurs journalières seulement avec SWin
ind_day=find(obs_FBt.SW_IN(i_obs) > 0.5 );
temp__day=temp4(ind_day,:);
% hour_=hour(obs_FBt.TIMESTAMP(i_obs));
% ind_day=find(hour_> 8& hour <19);
obs_=obs_FB(i_obs(ind_day),:);

% filtrer par mois 
month_d=month(obs_FBt.TIMESTAMP(i_obs));
%ind_w1=find(month_d==12);
%ind_w2=find(month_d<6);
%ind_sp=find(month_d >=3 & month_d<6);
ind_s=find(month_d >=6 & month_d<9 );
ind_s_day=intersect(ind_s,ind_day);
%ind_a=find(month_d >=9 & month_d<12);
 s=1;

figure_handle=figure('Name',['Scatter energy balance'])
 subplot(2,2,1)
score=vr_scoresstr(obs_.NETRAD ,temp__day(:,07));
score_RN=score;
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
 message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
hold on
    plot(obs_.NETRAD,temp__day(:,07),'b+')
      hold on
   plot(obs_FB.NETRAD(i_obs(ind_s_day)),temp4(ind_s_day,7),'r+')%    plot(obs_FB.NETRAD(ind_w1),temp_(ind_w1,7),'b+')
        plot(0:1000,0:1000,'k--','LineWidth',2)
        plot(0:1000,pente.*[0:1000]+ordorg,'k:','LineWidth',2)
        xlim([-100 1000])
        ylim([-100 1000])
        grid on 
        xlabel('Obs. (W.m^{-2})')
        ylabel('Sim. (W.m^{-2})')
       legend('All year','Summer','','Regression')
        title(['R_{NET} ' num2str(s)])
        text(40,800,message,"BackgroundColor",[1 1 1],'FontSize',12)
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

subplot (2,2,3)%
%score=vr_scoresstr(obs_.LE_hf ,temp__day(:,14));
score=vr_scoresstr(LE_corr(ind_day) ,temp__day(:,14));
 score_LE=score;
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
   plot(LE_corr(ind_day),temp__day(:,14),'b+')
   hold on 
  plot(obs_FB.LE_hf(i_obs(ind_s_day)),temp4(ind_s_day,14),'r+')
       hold on 
        plot(-100:800,-100:800,'k--','linewidth',2)
        plot(-100:800,pente.*[-100:800]+ordorg,'k:','linewidth',2)
        xlim([-100 600])
        ylim([-100 600])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['LE ' ])
   legend('All year','Summer','','Regression')
        grid on 
        text(-80,650,message,"BackgroundColor",[1 1 1],'FontSize',12)
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

subplot (2,2,4)%
%score=vr_scoresstr(obs_.H_hf ,temp__day(:,10));   
score=vr_scoresstr(H_corr(ind_day) ,temp__day(:,10));    
score_H=score;
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
   plot(H_corr(ind_day),temp__day(:,10),'b+')
   hold on 
  plot(H_corr(ind_s_day),temp4(ind_s_day,10),'r+')
       hold on 
       plot(-300:800,-300:800,'k--','LineWidth',2)
        plot(-300:800,pente.*[-300:800]+ordorg,'b:','LineWidth',2)
        xlim([-300 800])
        ylim([-300 800])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['H ' ])
        legend('All year','Summer','','Regression')
        grid on 
        text(-280,650,message,"BackgroundColor",[1 1 1],'FontSize',12)
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

subplot (2,2,2)%
score=vr_scoresstr(G0_obs(ind_day) ,temp__day(:,15));   
 score_G=score;
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
   plot(G0_obs(ind_day),temp__day(:,15),'b+')
   hold on 
  plot(G0_obs(ind_s_day),temp4(ind_s_day,15),'r+')
       hold on 
        plot(-120:120,-120:120,'k--','LineWidth',2)
        plot(-120:120,pente.*[-120:120]+ordorg,'k:','LineWidth',2)
        xlim([-100 120])
        ylim([-100 120])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['G ' ])
     legend('All year','Summer','','Regression')
        grid on 
        text(-110,60,message,"BackgroundColor",[1 1 1],'FontSize',12)
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 



%% Bilan de masse cumulé 

figure()
temp2(:,29)=temp2(:,29)-temp2(1,29);
temp2(:,23)=temp2(:,23)-temp2(1,23);
temp2(:,28)=temp2(:,28)-temp2(1,28);
temp2(:,30)=temp2(:,30)-temp2(1,30);
temp2(:,27)=temp2(:,27)-temp2(1,27);
temp4(:,29)=temp4(:,29)-temp4(1,29);
temp4(:,23)=temp4(:,23)-temp4(1,23);
temp4(:,28)=temp4(:,28)-temp4(1,28);
temp4(:,30)=temp4(:,30)-temp4(1,30);
temp4(:,27)=temp4(:,27)-temp4(1,27);

 figure('Name',['Mass balance' ] )
  plot(dates,temp2(:,29),'k-','linewidth',2)
  hold on
  errBar=transpose([temp2(:,23),temp4(:,23)]);
  %s=shadedErrorBar1(transpose(dates),transpose(temp2(:,23)),errBar,'lineprops','-r','transparent',true,'patchSaturation',0.3);
  hold on 
  plot(dates,temp2(:,23),'LineStyle','-','Color','r', 'LineWidth',2)
  hold on 
  plot(dates,temp4(:,23),'LineStyle','--','Color','r', 'LineWidth',2)


  hold on
  errBar=transpose([temp2(:,27),temp4(:,27)]);
  %s=shadedErrorBar1(transpose(dates),transpose(temp2(:,27)),errBar,'lineprops','-b','transparent',true,'patchSaturation',0.3);
  hold on 
  plot(dates,temp2(:,27),'LineStyle','-','Color',[0.07,0.62,1.00], 'LineWidth',2)
  hold on 
  plot(dates,temp4(:,27),'LineStyle','--','Color',[0.07,0.62,1.00], 'LineWidth',2)
  %s.patch.FaceColor = [0.07,0.62,1.00];
  
  hold on
  errBar=transpose([temp2(:,30),temp4(:,30)]);
  %shadedErrorBar1(transpose(dates),transpose(temp2(:,30)),errBar,'lineprops','-g','transparent',true,'patchSaturation',0.3);
  hold on 
  plot(dates,temp2(:,30),dates,temp4(:,30),'LineStyle','-','Color','g', 'LineWidth',2)
  
   hold on
  errBar=transpose([temp2(:,28),temp4(:,28)]);
  %shadedErrorBar1(transpose(dates),transpose(temp2(:,28)),errBar,'lineprops','--b','transparent',true,'patchSaturation',0.3);
  hold on 
  plot(dates,temp2(:,28),'LineStyle','-','Color','b', 'LineWidth',2)
  hold on 
  plot(dates,temp4(:,28),'LineStyle','--','Color','b', 'LineWidth',2)
  
  
  
 
% dstock2=temp2(:,29)-temp2(:,23)-temp2(:,27)-temp2(:,28)-temp2(:,30)+temp2(1,30);
% dstock4=temp4(:,29)-temp4(:,23)-temp4(:,27)-temp4(:,28)-temp4(:,30)+temp4(1,30);
%     hold on
%   errBar=transpose([dstock2,dstock4]);
%   %shadedErrorBar1(transpose(dates),transpose(dstock2),errBar,'lineprops','--m','transparent',true,'patchSaturation',0.3);
%   hold on 
%   plot(dates,dstock2,'LineStyle','-','Color','m', 'LineWidth',2)
%   hold on 
%   plot(dates,dstock4,'LineStyle','--','Color','m', 'LineWidth',2)

%         plot(dates,t(:,29)-t(:,23)-t(:,27)-t(:,28)-t(:,30)+t(1,30),'m-','linewidth',2)
        dynamicDateTicks(); datetick('x','mmm-yy','keepticks')
        set(gca, 'YGrid', 'on', 'XGrid', 'off')
        ylabel('mm')
        title(['Bilan de masse L'])
        legend('rain','ET','','RUN-OFF','','','d STOCK','PERCOL','','','location','northwest')
        
        %% Température avec des valeurs de conductivité différentes
        
figure_handle=figure('Name','Tsol')
subplot(2,1,1)
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,36),'r-','linewidth',1)
%   
%     hold on
%     plot(dates4,temp4(:,36),'b-','linewidth',1)
    ylim([0 30])
     dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('T (°C)')
    title(['T soil 5 cm ' ])
    grid on
   
        ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

    subplot(2,1,2)
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_50_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,40),'r-','linewidth',1.5)
%        hold on
%     plot(dates4,temp4(:,40),'b-','linewidth',1.5)
    ylim([0 30])
     dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('T (°C)')
    title(['Tsoil 50 cm ' ])
    grid on
    legend('Observed','Simulated','location','northeast')
    all_ha = findobj( figure_handle, 'type', 'axes', 'tag' ,'');
linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

%% Humidités dans le sol + pluies 
dates_d=datetime(dates,'ConvertFrom','datenum');
    figure_handle=figure('Name','teneur en eau ')
    %figure(6)
    subplot(5,1,1)
    Precip= Pcumuljour;
    Pcumuljourtime=datenum(Pcumuljour.TIMESTAMP);
     ind_t_Ps=find(Pcumuljourtime==dates(1));
     ind_t_Pe=find(Pcumuljourtime==floor(dates(end)));
dates_p=Pcumuljourtime(ind_t_Ps:ind_t_Pe);
     bar(dates_p,Precip.P(ind_t_Ps:ind_t_Pe))
     hold on 
     grid on 
%dynamicDateTicks(); 
%datetick('x','mmm-yy')
     legend('Dayly rainfall')
 ylabel('P (mm')   

       dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
    ax.YDir = 'reverse';
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
   
  subplot(5,1,2)  
  
  plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_5_1(ind_t),'k-','linewidth',1.5)
     hold on
    plot(dates,temp2(:,45)/100,'b-','linewidth',1.5)
   ylim([0.079,0.666])
    datetick('x')
    grid on 
        ylabel('\theta (m^3.m^{-3})')
    hold on   
grid on 
 dynamicDateTicks(); 
datetick('x','mmm-yy','keepticks','keeplimits')
    title(['5cm ' ],'FontSize',15)
     legend('Observed ','Simulated ','location','north')
       dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
   
    
    subplot(5,1,3)

    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_25_2(ind_t),'k-','linewidth',1.5)
    
    hold on
    plot(dates,temp2(:,47)/100,'b-','linewidth',1.5)
   ylim([0.079,0.666])
    datetick('x')
    grid on 
        ylabel('\theta (m^3.m^{-3})')
    title(['25 cm' ],'FontSize',15)
      dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    
    
    subplot(5,1,4)
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_50_2(ind_t)-0.2,'k-','linewidth',1.5)
     hold on
     plot(dates,temp2(:,49)/100,'b-','linewidth',1.5)
    ylim([0.079,0.666])
    grid on 
    datetick('x')
    ylabel('\theta (m^3.m^{-3})')
    title(['50 cm' ],'FontSize',15)
      dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
  
    
   subplot(5,1,5)
    plot(dates,temp2(:,52)/100,'b-','linewidth',1.5)
    hold on
   ylim([0.079,0.666])
    datetick('x')
    ylabel('\theta (m^3.m^{-3})')
    grid on 
    title(['150 cm ' ],'FontSize',15)
    %legend('obs','sim 4','sim 3','location','eastoutside')
      dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')    
ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    
    all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );


        
        %% subplot Rnet SWout et LW out 
        
 figure_handle= figure ()
subplot(4,1,1)
plot(dates,temp2(:,07),'r-','linewidth',2)
hold on 
plot(dates3,temp3(:,07),'m-','linewidth',2)
hold on 
plot(dates4,temp4(:,07),'b-','linewidth',2)
hold on
scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.NETRAD(ind_t),5,'k','filled')
datetick('x')
ylabel('W.m^{-2}')
title(['Rn ' ])
grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','location','eastoutside')

subplot(4,1,2)
plot(dates,sim_SWout1,'r-','linewidth',2)
hold on
plot(dates,sim_SWout2,'m-','linewidth',2)
hold on
plot(dates,sim_SWout3,'b-','linewidth',2)
hold on
scatter(obs_FB.TIMESTAMP(i_obs),obs_FB.SW_OUT(i_obs),5,'k','filled')
datetick('x')
ylabel('W.m^{-2}')
title(['SW out ' ])
grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','location','eastoutside')

    
sim_LW_out1=temp2(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout1-obs_FB.LW_IN(i_obs);
sim_LW_out2=temp3(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout2-obs_FB.LW_IN(i_obs);
sim_LW_out3=temp4(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout3-obs_FB.LW_IN(i_obs);
    
subplot(4,1,3)

plot(dates,-sim_LW_out1,'r-','linewidth',2)
    hold on 
    plot(dates3,-sim_LW_out2,'m-','linewidth',2)
    hold on 
    plot(dates4,-sim_LW_out3,'b-','linewidth',2)
    hold on
    LWout=obs_FB.LW_OUT(ind_t);
    LWout2=obs_FB.NETRAD(ind_t)-obs_FB.SW_IN(ind_t)+obs_FB.SW_OUT(ind_t)-obs_FB.LW_IN(ind_t)
    scatter(obs_FB.TIMESTAMP(ind_t),-LWout2,5,'r','filled')
    hold on 
    scatter(obs_FB.TIMESTAMP(ind_t),LWout,5,'k','filled')
    datetick('x')
    ylabel('W.m^{-2}')
    title(['LW out ' ])
    grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','obs','location','eastoutside')

    subplot(4,1,4)
    plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',2)
    hold on
    plot(dates4,temp2(:,65)*100,'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,65)*100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,65)*100,'b-','linewidth',1)
    datetick('x')
    grid on
    ylabel('%')
    title(['Albedo ' ])
    legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
    

all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots

%% partition evap transpi potentiel foliaire

   figure('Name','potentiel foliaire')
    %plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',1.5)
    hold on
    plot(dates4,temp4(:,18)/101.9,'b-','linewidth',1.5)
    dynamicDateTicks();
   
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('Potentiel hydrique foliaire (Mpa)')
    title([' '])
    %legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

%ET=cumsum(((t(:,14))- t(1,14))*dt/Lv,'omitnan');
figure()
t=temp2;
  plot(dates,t(:,29)-t(1,29),'k-','linewidth',2)
        hold on
     plot(dates,t(:,25)-t(1,25),'b-','linewidth',2)
      hold on
     plot(dates,t(:,24)-t(1,24),'r-','linewidth',2)
      hold on
     plot(dates,t(:,26)-t(1,26),'g-','linewidth',2)
 hold on

plot(dates(122736),t(:,23)-t(1,23))
 legend('pluie','Evap sol','Transpi','Evap feuille')
 dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    

%% scores pour tableau 
 
s_RN_1=[scoresRN1.rmse, scoresRN1.mbias, scoresRN1.cor2];
s_LE_1=[scoresLE1.rmse, scoresLE1.mbias, scoresLE1.cor2];
s_H_1=[scoresH1.rmse scoresH1.mbias scoresH1.cor2];
s_G_1=[scoresG1.rmse scoresG1.mbias scoresG1.cor2];
s_T5_1=[scoresT5.sim_1.rmse scoresT5.sim_1.mbias scoresT5.sim_1.cor2];
s_T15_1=[scoresT15.sim_1.rmse scoresT15.sim_1.mbias scoresT15.sim_1.cor2];
s_T25_1=[scoresT25.sim_1.rmse scoresT25.sim_1.mbias scoresT25.sim_1.cor2];
s_T35_1=[scoresT35.sim_1.rmse scoresT35.sim_1.mbias scoresT35.sim_1.cor2];
s_T50_1=[scoresT50.sim_1.rmse scoresT50.sim_1.mbias scoresT50.sim_1.cor2];
s_W5_1=[scoresW5.sim_1.rmse scoresW5.sim_1.mbias scoresW5.sim_1.cor2];
s_W15_1=[scoresW15.sim_1.rmse scoresW15.sim_1.mbias scoresW15.sim_1.cor2];
s_W25_1=[scoresW25.sim_1.rmse scoresW25.sim_1.mbias scoresW25.sim_1.cor2];
s_W35_1=[scoresW35.sim_1.rmse scoresW35.sim_1.mbias scoresW35.sim_1.cor2];
s_W50_1=[scoresW50.sim_1.rmse scoresW50.sim_1.mbias scoresW50.sim_1.cor2];

s_RN_2=[scoresRN2.rmse, scoresRN2.mbias, scoresRN2.cor2];
s_LE_2=[scoresLE2.rmse, scoresLE2.mbias, scoresLE2.cor2];
s_H_2=[scoresH2.rmse scoresH2.mbias scoresH2.cor2];
s_G_2=[scoresG2.rmse scoresG2.mbias scoresG2.cor2];
s_T5_2=[scoresT5.sim_2.rmse scoresT5.sim_2.mbias scoresT5.sim_2.cor2];
s_T15_2=[scoresT15.sim_2.rmse scoresT15.sim_2.mbias scoresT15.sim_2.cor2];
s_T25_2=[scoresT25.sim_2.rmse scoresT25.sim_2.mbias scoresT25.sim_2.cor2];
s_T35_2=[scoresT35.sim_2.rmse scoresT35.sim_2.mbias scoresT35.sim_2.cor2];
s_T50_2=[scoresT50.sim_2.rmse scoresT50.sim_2.mbias scoresT50.sim_2.cor2];
s_W5_2=[scoresW5.sim_2.rmse scoresW5.sim_2.mbias scoresW5.sim_2.cor2];
s_W15_2=[scoresW15.sim_2.rmse scoresW15.sim_2.mbias scoresW15.sim_2.cor2];
s_W25_2=[scoresW25.sim_2.rmse scoresW25.sim_2.mbias scoresW25.sim_2.cor2];
s_W35_2=[scoresW35.sim_2.rmse scoresW35.sim_2.mbias scoresW35.sim_2.cor2];
s_W50_2=[scoresW50.sim_2.rmse scoresW50.sim_2.mbias scoresW50.sim_2.cor2];


s_RN_3=[scoresRN3.rmse, scoresRN3.mbias, scoresRN3.cor2];
s_LE_3=[scoresLE3.rmse, scoresLE3.mbias, scoresLE3.cor2];
s_H_3=[scoresH3.rmse scoresH3.mbias scoresH3.cor2];
s_G_3=[scoresG3.rmse scoresG3.mbias scoresG3.cor2];
s_T5_3=[scoresT5.sim_3.rmse scoresT5.sim_3.mbias scoresT5.sim_3.cor2];
s_T15_3=[scoresT15.sim_3.rmse scoresT15.sim_3.mbias scoresT15.sim_3.cor2];
s_T25_3=[scoresT25.sim_3.rmse scoresT25.sim_3.mbias scoresT25.sim_3.cor2];
s_T35_3=[scoresT35.sim_3.rmse scoresT35.sim_3.mbias scoresT35.sim_3.cor2];
s_T50_3=[scoresT50.sim_3.rmse scoresT50.sim_3.mbias scoresT50.sim_3.cor2];
s_W5_3=[scoresW5.sim_3.rmse scoresW5.sim_3.mbias scoresW5.sim_3.cor2];
s_W15_3=[scoresW15.sim_3.rmse scoresW15.sim_3.mbias scoresW15.sim_3.cor2];
s_W25_3=[scoresW25.sim_3.rmse scoresW25.sim_3.mbias scoresW25.sim_3.cor2];
s_W35_3=[scoresW35.sim_3.rmse scoresW35.sim_3.mbias scoresW35.sim_3.cor2];
s_W50_3=[scoresW50.sim_3.rmse scoresW50.sim_3.mbias scoresW50.sim_3.cor2];

tab_scores_1=[s_RN_1;s_LE_1;s_H_1;s_G_1;s_T5_1;s_T15_1;s_T25_1;s_T35_1;s_T50_1;s_W5_1;s_W15_1;s_W25_1;s_W35_1;s_W50_1];
tab_scores_2=[s_RN_2;s_LE_2;s_H_2;s_G_2;s_T5_2;s_T15_2;s_T25_2;s_T35_2;s_T50_2;s_W5_2;s_W15_2;s_W25_2;s_W35_2;s_W50_2];
tab_scores_3=[s_RN_3;s_LE_3;s_H_3;s_G_3;s_T5_3;s_T15_3;s_T25_3;s_T35_3;s_T50_3;s_W5_3;s_W15_3;s_W25_3;s_W35_3;s_W50_3];
%% profil des pourcentages de cailloux 
    p=[2,5,15,25,35,50,75,100,150,200];
    cailloux=[0.15,0.5,0.6,0.65,0.7,0.8,0.85,0.9,0.9,0.9]; % pourcentage de sol6
    figure()
    plot(cailloux*100,p)
   
%% climato temperature pluies
fid=figure()
StartDate = datetime('2013-01-01');
EndDate =datetime('2013-12-31'); 
StartDate.TimeZone='Z';
EndDate.TimeZone='Z';
 xData = linspace(StartDate,EndDate,12);
%ax = gca;
% ax.XTick = xData;
Date = StartDate:caldays(1):EndDate;
Date.Format = 'dd-MMM';
x=Date';


subplot(2,1,1)
y = transpose(Moyenneparjour.TA);
errBar=transpose([Minparjour.TA,Maxparjour.TA]);
hold on 
shadedErrorBar1(x,y,errBar,'lineprops','-b','transparent',true,'patchSaturation',0.3);
hold on 
plot(x,y,'Color','k')
hold on 
%plot(x,Minparjour.TA,x,Maxparjour.TA,'LineStyle','--','Color','1.00,0.59,0.01','Marker','o','MarkerEdgeColor','none','MarkerSize',1.8,'MarkerFaceColor','1.00,0.59,0.01');
ylabel('Temperature ( °C )','FontWeight','bold')
%plot(Moyenneparjour.TIMESTAMP,Minparjour.TA,Moyenneparjour.TIMESTAMP,Maxparjour.TA,'LineStyle','--','Color','b','Marker','o','MarkerEdgeColor','none','MarkerSize',1.8,'MarkerFaceColor','blue');
%title("Température de l'air (TA) ")
%xlabel('FontWeight','bold')
grid on
ax = gca;
legend
ax.FontWeight='bold';
%set(gca,'fontsize',12,'fontweight','bold')

subplot(2,1,2)
x = ['Jan','Feb','Mar','Apr', 'May', 'Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
ymoy = [70 58 50 52 50 28 12 11 72 86 110 90];
ymax=[205 180 105 90 100 73 50 90 138 300 220 182];
ymin=[10 15 0 15 0 1 0 0 1 0 30 24];
X = categorical({'Jan','Feb','Mar','Apr', 'May', 'Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
X = reordercats(X,{'Jan','Feb','Mar','Apr', 'May', 'Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
bar(X,ymoy)
hold on 
ylabel('Precipitations (mm)','FontWeight','bold')
scatter(X,ymin,'Marker','_')
hold on
scatter(X,ymax,'Marker','_')
grid on
legend('Mean monthly ', 'Min-Max monthly' ,'')
ax = gca;
legend
ax.FontWeight='bold';

%% Transpi et evap cumulées comparaison version karst 

 ind_LE=isnan(obs_FB.LE_hf(ind_t));
    [Ebis2,Ebis3,Ebis4]=deal(  temp2(:,25),temp3(:,25),temp4(:,25));
    [Tbis2,Tbis3,Tbis4]=deal(  temp2(:,24),temp3(:,24),temp4(:,24));
      for i=1:length(obs_FB.LE_hf(ind_t))
        if ind_LE(i) ==1
 
            
             Ebis2(i)=NaN;
             Ebis4(i)=NaN;
             Ebis3(i)=NaN;
             Tbis4(i)=NaN;
              Tbis2(i)=NaN;
              Tbis3(i)=NaN;
        end 
      end 
      % période choisie mars septembre 2019 ici 
 %améliorer et déterminer l'indice à partir d'une date 
T1=Tbis2(1700*48)-Tbis2(1502*48);
E1=Ebis2(1700*48)-Ebis2(1502*48);
T2=Tbis3(1700*48)-Tbis3(1502*48);
E2=Ebis3(1700*48)-Ebis3(1502*48);
T3=Tbis4(1700*48)-Tbis4(1502*48);
E3=Ebis4(1700*48)-Ebis4(1502*48);
ob=sum(obs_FB.LE_hf(ind_t(1502*48:1700*48))*dt/Lv,'omitnan')
   
ET2=Tbis2(1502*48:1700*48)+Ebis2(1502*48:1700*48);
ET3=Tbis3(1502*48:1700*48)+Ebis3(1502*48:1700*48);
ET4=Tbis4(1502*48:1700*48)+Ebis4(1502*48:1700*48);
etobs=obs_FB.LE_hf(ind_t(1502*48:1700*48))*dt/Lv;

e2=vr_scoresstr(ETR_day_obs(1502:1700),ETR_day_sim1(1502:1700));
e3=vr_scoresstr(ETR_day_obs(1502:1700),ETR_day_sim2(1502:1700));
e4=vr_scoresstr(ETR_day_obs(1502:1700),ETR_day_sim3(1502:1700));
% T1=temp2(1700*48,24)-temp2(1502*48,24)
% E1=temp2(1700*48,25)-temp2(1502*48,25)
% T2=temp4(1700*48,24)-temp4(1502*48,24);
% E2=temp4(1700*48,25)-temp4(1502*48,25);
y=[T1 E1 0; T2 E2 0; T3 E3 0; 0 0 ob];
X = categorical({'Standard' ,'Karst1','Karst2', 'Obs.'});
X = reordercats(X,{'Standard' ,'Karst1','Karst2', 'Obs.'});
figure ()
b=bar(X,y,'stacked','FaceColor','flat')
%b.CData(3) = 'k';
ylabel('mm')

%pour regarder la chronique sur la periode considérée 
figure()
plot (dates(1502*48:1700*48),temp2(1502*48:1700*48,24))
dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')

    
%% SUBPLOT TENEUR EN EAU A PLUSIEURS PROFONDEURS ET EVAPOTRANSPI 

 [ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),obs_FB.LE(ind_t),'day',@nansum);
    ETR_day_obs=ETR_day_obs*dt/Lv;
    [ETR_day_sim1,days_sim1]=timeavgCK(dates,temp2(:,14),'day',@nansum);
    [ETR_day_sim2,days_sim2]=timeavgCK(dates3,temp3(:,14),'day',@nansum);
    [ETR_day_sim3,days_sim3]=timeavgCK(dates4,temp4(:,14),'day',@nansum);
    [ETR_day_sim1,ETR_day_sim2,ETR_day_sim3]=deal( ETR_day_sim1*dt/Lv, ETR_day_sim2*dt/Lv, ETR_day_sim3*dt/Lv);
    for i=1:length(ETR_day_obs)
        if N_LE_obs(i) <40
            ETR_day_obs(i)=NaN;
            ETR_day_obs3(i)=NaN;
           
        end 
    end 


   figure_handle=figure('Name','teneur en eau ')
   
     subplot(4,1,1)  
 plot(days_sim1,ETR_day_sim1,'b-','linewidth',2)
    hold on 
    plot(days_sim3,ETR_day_sim2,'m-','linewidth',2)
    hold on 
    plot(days_sim2,ETR_day_sim3,'r-','linewidth',2) 
    hold on 
    plot(days_obs,ETR_day_obs,'k-','linewidth',2)
    grid on
    ylabel('Daily ET (mm)')
     dynamicDateTicks();
          ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    datetick('x','mmm-yy','keepticks')
%        % legend(['sim1: ' num2str(sum(ETR_day_sim1,'omitnan')) 'mm'],['sim2: ' num2str(sum(ETR_day_sim2,'omitnan')) 'mm'],['sim3: ' num2str(sum(ETR_day_sim3,'omitnan')) 'mm'],...
%     ['obs: ' num2str(sum(ETR_day_obs,'omitnan')) 'mm'],'location','eastoutside');
% legend('Shallow (2m)','Deep (20m)','Karst (20m+cracks)')
% %    
  subplot(4,1,2)  
     plot(dates,temp2(:,45)/100,'b-','linewidth',2)
    hold on
    plot(dates3,temp3(:,45)/100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,45)/100,'r-','linewidth',2)
    hold on
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_5_1(ind_t),'k-','linewidth',1)
    %ylim([0 0.14])
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('\theta (m^3.m^{-3})')
    title(['SMC 5 ' ])
     %legend('sim ref','sim2','sim3','obs','location','eastoutside')
     ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold';  
   
    
    subplot(4,1,3)

  hold on
    plot(dates,temp2(:,47)/100,'b-','linewidth',1.5)
    hold on
    plot(dates3,temp3(:,47)/100,'m-','linewidth',1.5)
    hold on
   plot(dates4,temp4(:,47)/100,'r-','linewidth',1.5)
    hold on 
     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_25_2(ind_t),'k-','linewidth',1.5)
    %ylim([0 0.14])
    grid on 
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('\theta (m^3.m^{-3})')
    title(['25 cm' ])
    legend('sim ref','sim2','sim3','obs','location','east')
    %legend('/theta r min /theta s min', '/theta r max /theta s min','/theta r max /theta s max', 'obs')
    ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    
    
    subplot(4,1,4)
     hold on
    plot(dates,temp2(:,49)/100,'b-','linewidth',2)
    hold on
    plot(dates3,temp3(:,49)/100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,49)/100,'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_50_2(ind_t)-0.2,'k-','linewidth',1)
    grid on 
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks','keeplimits')
    ylabel('\theta (m^3.m^{-3})')
    title(['SMC 50 '])
    %legend('sim ref','sim2','sim3','obs','location','eastoutside')
    ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

  
    
all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );


        
        %% subplot Rnet SWout et LW out 
        
 figure_handle= figure ()
subplot(4,1,1)
plot(dates,temp2(:,07),'r-','linewidth',2)
hold on 
plot(dates3,temp3(:,07),'m-','linewidth',2)
hold on 
plot(dates4,temp4(:,07),'b-','linewidth',2)
hold on
scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.NETRAD(ind_t),5,'k','filled')
datetick('x')
ylabel('W.m^{-2}')
title(['Rn ' ])
grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','location','eastoutside')

subplot(4,1,2)
plot(dates,sim_SWout1,'r-','linewidth',2)
hold on
plot(dates,sim_SWout2,'m-','linewidth',2)
hold on
plot(dates,sim_SWout3,'b-','linewidth',2)
hold on
scatter(obs_FB.TIMESTAMP(i_obs),obs_FB.SW_OUT(i_obs),5,'k','filled')
datetick('x')
ylabel('W.m^{-2}')
title(['SW out ' ])
grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','location','eastoutside')

    
sim_LW_out1=temp2(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout1-obs_FB.LW_IN(i_obs);
sim_LW_out2=temp3(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout2-obs_FB.LW_IN(i_obs);
sim_LW_out3=temp4(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout3-obs_FB.LW_IN(i_obs);
    
subplot(4,1,3)

% plot(dates,-sim_LW_out1,'r-','linewidth',2)
    hold on 
    plot(dates3,-sim_LW_out2,'m-','linewidth',2)
    hold on 
    plot(dates4,-sim_LW_out3,'b-','linewidth',2)
    hold on
    LWout=obs_FB.LW_OUT(ind_t);
    LWout2=obs_FB.NETRAD(ind_t)-obs_FB.SW_IN(ind_t)+obs_FB.SW_OUT(ind_t)-obs_FB.LW_IN(ind_t);
    scatter(obs_FB.TIMESTAMP(ind_t),-LWout2,5,'r','filled')
    hold on 
    scatter(obs_FB.TIMESTAMP(ind_t),LWout,5,'k','filled')
    datetick('x')
    ylabel('W.m^{-2}')
    title(['LW out ' ])
    grid on
legend('sim L1 ref','sim L1B','sim L1C','obs','obs','location','eastoutside')

    subplot(4,1,4)
    plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',2)
    hold on
    plot(dates4,temp2(:,65)*100,'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,65)*100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,65)*100,'b-','linewidth',1)
    datetick('x')
    grid on
    ylabel('%')
    title(['Albedo ' ])
    legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
    

all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
