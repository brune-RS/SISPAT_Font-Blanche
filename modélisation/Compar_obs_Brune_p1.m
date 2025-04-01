%% Compar_obs_.m
% Ce script lit un (ou plusieurs) fichier atm de sortie de SiSPAT et le compare aux obs de
% février- juin 2023 Brune
% pour données de FONT-BLANCHE
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
path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/K0bis/';
path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/K1d6s1/'; 
path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/tests/a2/';
%path_res2='/home/brune/SISPAT/execution/sortie_modele/FB/fev/runoff1/';
%path_res3='/home/brune/SISPAT/execution/sortie_modele/FB/fev/runoff2/';
%path_res4='/home/brune/SISPAT/execution/sortie_modele/FB/fev/runoff3/';
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
    %% Calcul de scores
    %scoresSWout=vr_scoresstr((obs_FB.SW_IN(i_obs)).*(obs_FB.ALB(i_obs)),(obs_FB.SW_IN(i_obs)).*(temp(i_sim,65)));
    scoresRN1=vr_scoresstr(obs_FB.NETRAD(i_obs),temp2(i_sim,7));
    scoresRN2=vr_scoresstr(obs_FB.NETRAD(i_obs),temp3(i_sim,7));
    scoresRN3=vr_scoresstr(obs_FB.NETRAD(i_obs),temp4(i_sim,7));
    
    scoresbLE1=vr_scoresstr(LE_corr,temp2(i_sim,14));
    scoresbLE2=vr_scoresstr(LE_corr,temp3(i_sim,14));
    scoresbLE3=vr_scoresstr(LE_corr,temp4(i_sim,14));
    scoresLE1=vr_scoresstr(obs_FB.LE_hf(i_obs),temp2(i_sim,14));
    scoresLE2=vr_scoresstr(obs_FB.LE_hf(i_obs),temp3(i_sim,14));
    scoresLE3=vr_scoresstr(obs_FB.LE_hf(i_obs),temp4(i_sim,14));
    
    scoresbH1=vr_scoresstr(H_corr,temp2(i_sim,10));
    scoresbH2=vr_scoresstr(H_corr,temp3(i_sim,10));
    scoresbH3=vr_scoresstr(H_corr,temp4(i_sim,10));
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
%% Graphiques :albedo, Rnet, SWout, LWout (chroniques et scatterplot+stats)
    % Albedo
   figure('Name','Albedo')
    plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',1.5)
    hold on
    plot(dates4,temp4(:,65)*100,'b-','linewidth',1.5)
    dynamicDateTicks();
    hold on
    plot(dates3,temp3(:,65)*100,'m-','linewidth',2)
    hold on
   plot(dates,temp2(:,65)*100,'r-','linewidth',2) 
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('Albedo (%)')
    title([' '])
    legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
  



    % Rn
    figure('Name','Rn')
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
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    
    figure('Name',['Rn scatter' ])
    for s=1:3
       subplot(1,3,s)
        score=eval(['scoresRN',num2str(s)]);
        tempb=eval(['temp',num2str(s+1)]);
%         [ia0,ia5]=deal(  find(G0j_sim==A0),find(G5j_sim==A5));
        [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
         message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
         hold on
        plot(obs_FB.NETRAD(i_obs),tempb(i_sim,07),'r+')
        plot(0:1200,0:1200,'k--')
        plot(0:1200,pente.*[0:1200]+ordorg,'b:')
        xlim([-100 1200])
        ylim([-100 1200])
        grid on 
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['Rn ' num2str(s)])
        text(50,1000,message)
    end
%     s=1;
%      figure('Name',['R_{NET}' ])
%  score=eval(['scoresRN',num2str(s)]);
%         tempb=eval(['temp',num2str(s+1)]);
% %         [ia0,ia5]=deal(  find(G0j_sim==A0),find(G5j_sim==A5));
%         [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
%          message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
%          hold on
%         plot(obs_FB.NETRAD(i_obs),tempb(i_sim,07),'b+')
%         plot(-100:1200,-100:1200,'k--','LineWidth',2)
%         plot(-100:1200,pente.*[-100:1200]+ordorg,'k:','LineWidth',2)
%         xlim([-100 1000])
%         ylim([-100 1000])
%         grid on 
%         xlabel('Obs. W.m^{-2}')
%         ylabel('Sim. W.m^{-2}')
%         title(['R_{NET} ' num2str(s)])
%         text(50,800,message)
%             ax = gca;
% ax.FontSize = 10; 
% ax.FontWeight = 'bold'; 
       %SW out
   
    figure('Name',['Sw out scatter' ])
    %obs_SWout=(obs_FB.SW_IN(i_obs)).*(obs_FB.ALB(i_obs));
    sim_SWout1=(obs_FB.SW_IN(i_obs)).*(temp2(i_sim,65));
    sim_SWout2=(obs_FB.SW_IN(i_obs)).*(temp3(i_sim,65));
    sim_SWout3=(obs_FB.SW_IN(i_obs)).*(temp4(i_sim,65));
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        sim_SWout=eval(['sim_SWout',num2str(s)]);
        scoresSWout=vr_scoresstr(obs_FB.SW_OUT(i_obs),(obs_FB.SW_IN(i_obs)).*(temp(i_sim,65)));   
        subplot(1,3,s)
        score=eval(['scoresSWout']);
         [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.SW_OUT(i_obs),sim_SWout,'r+')
        plot(0:150,0:150,'k--','LineWidth',2)
        plot(0:150,pente.*[0:150]+ordorg,'b:','LineWidth',2)
        xlim([-10 150])
        ylim([-10 150])
        grid on 
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['Sw out ' num2str(s)])
        text(0,120,message)    
        
    end
    
    figure('Name','SWout')
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
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    %IR out  (lW out)
    sim_LW_out1=temp2(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout1-obs_FB.LW_IN(i_obs);
    sim_LW_out2=temp3(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout2-obs_FB.LW_IN(i_obs);
    sim_LW_out3=temp4(:,07)-obs_FB.SW_IN(i_obs)+sim_SWout3-obs_FB.LW_IN(i_obs);
    figure('Name','LWout')
    plot(dates,-sim_LW_out1,'r-','linewidth',2)
    hold on 
    plot(dates3,-sim_LW_out2,'m-','linewidth',2)
    hold on 
    plot(dates4,-sim_LW_out3,'b-','linewidth',2)
    hold on
    LWout=obs_FB.LW_OUT(ind_t);
    scatter(obs_FB.TIMESTAMP(ind_t),LWout,5,'k','filled')
    datetick('x')
    ylabel('W.m^{-2}')
    title(['LW out ' ])
    grid on
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    LWnbis=LWout(1:35087);
    scoresLWout1=vr_scoresstr(LWout,(-sim_LW_out1));
    scoresLWout2=vr_scoresstr(LWout,(-sim_LW_out2));
    scoresLWout3=vr_scoresstr(LWout,(-sim_LW_out3));
    figure('Name',['LWout scatter' ])
    for s=1:3
       subplot(1,3,s)
        score=eval(['scoresLWout',num2str(s)]);
        sim_LW_out=eval(['sim_LW_out',num2str(s)]);
        [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
         hold on
        plot(LWout,-sim_LW_out,'r+')
        plot(250:600,250:600,'k--','LineWidth',2)
        plot(250:600,pente.*[250:600]+ordorg,'b:','LineWidth',2)
        xlim([250 600])
        ylim([250 600])
        grid on 
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['LWout ' num2str(s)])
        text(260,550,message)
    end
    
%close all
%% SUBPLOT BILIAN RADIATIF SURFACE ( Rnet, SWout, LWout, albedo)

figure_handle= figure ()
subplot(4,1,1)
plot(dates,temp2(:,07),'r-','linewidth',2)
hold on 
plot(dates3,temp3(:,07),'m-','linewidth',2)
hold on 
plot(dates4,temp4(:,07),'b-','linewidth',2)
hold on
scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.NETRAD(ind_t),5,'k','filled')
dynamicDateTicks();
datetick('x','mmm-yy','keepticks')
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
dynamicDateTicks();
datetick('x','mmm-yy','keepticks')
ylabel('W.m^{-2}')
title(['SW out ' ])
grid on
%legend('sim L1 ref','sim L1B','sim L1C','obs','location','eastoutside')

    
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
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('W.m^{-2}')
    title(['LW out ' ])
    grid on
%legend('sim L1 ref','sim L1B','sim L1C','obs','obs','location','eastoutside')

    subplot(4,1,4)
    plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',2)
    hold on
    plot(dates4,temp2(:,65)*100,'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,65)*100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,65)*100,'b-','linewidth',1)
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('%')
    title(['Albedo ' ])
    %legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
    

all_ha = findobj( figure_handle, 'type', 'axes', 'tag' ,'');
linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
     
%% Graphes Bilan d'énergie (Gsurface, G5cm, H, LE)   
    % G 5cm
    figure('Name','G')
    plot(dates,temp2(:,67),'r-','linewidth',2)
    hold on 
    plot(dates3,temp3(:,67),'m-','linewidth',2)
    hold on 
    plot(dates4,temp4(:,67),'b-','linewidth',2)
    hold on 
    scatter(obs_FB.TIMESTAMP(ind_t),G_obs(ind_t),5,'k','filled')
    dynamicDateTicks('','','dd/mm');
    datetick('x','mmm-yy','keepticks')
    ylabel('G 5cm (W.m^{-2})','fontweight','bold')
    %title(['G 5cm ' ])
    grid on 
    %legend('Simulation','Obs','location','eastoutside')
      legend('sim ref','sim2','sim3','obs','location','eastoutside')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

    figure('Name',['G scatter' ])
    for s=1:3
        score=eval(['scoresG',num2str(s)]);
        subplot(1,3,s)
        temp=eval(['temp',num2str(s+1)]);
        [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(G_obs(i_obs),temp(i_sim,67),'r+')
        hold on
        plot(-100:300,-100:300,'k--','linewidth',2)
        plot(-100:300,pente.*[-100:300]+ordorg,'b:','linewidth',2)
        xlim([-40 40])
        ylim([-40 40])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['G ' num2str(s)])
        text(-35,30,message)
    end    
   
    
    % G surf
    figure('Name','G surf')
    plot(dates,temp2(:,15),'r-','linewidth',2)
    hold on 
    plot(dates3,temp3(:,15),'m-','linewidth',2)
    hold on 
    plot(dates4,temp4(:,15),'b-','linewidth',2)
    hold on
    scatter(obs_FB.TIMESTAMP(ind_t),G0_obs,5,'k','filled')
    datetick('x')
    ylabel('W.m^{-2}')
    title(['G surf ' ])
    grid on
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    figure('Name',['G surf scatter' ])
    for s=1:3
        subplot(1,3,s)
        score=eval(['scoresG0',num2str(s)]);
        temp_=eval(['temp',num2str(s+1)]);
        [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(G_obs(i_obs),temp_(i_sim,67),'r+')
        hold on
        plot(-50:50,-50:50,'k--','LineWidth',2)
        plot(-50:50,pente.*[-50:50]+ordorg,'b:','LineWidth',2)
        xlim([-50 50])
        ylim([-50 50])
        grid on 
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['G ' num2str(s)])
        text(-40,35,message)
    end    

    
    % H
    figure('Name','H')
    plot(dates,temp2(:,10),'r-','linewidth',2)
    hold on 
    %plot(dates3,temp3(:,10),'m-','linewidth',2)
     plot(dates3,H_corr,'m-','linewidth',2)
    hold on 
    plot(dates4,temp4(:,10),'b-','linewidth',2)
    %scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.H(i_obs),5,'r','filled')
    hold on
    scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.H_hf(i_obs),5,'k','filled')
    hold on
    
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('W.m^{-2}')
    title(['H ' ])
    grid on
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    figure('Name',['H scatter' ])
    for s=1:3
       subplot(1,3,s)
        score=eval(['scoresbH',num2str(s)]);
        temp_=eval(['temp',num2str(s+1)]);
         [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
         hold on
        plot(obs_FB.H_hf(ind_t),temp_(i_sim,10),'r+')
        plot(-300:800,-300:800,'k--','LineWidth',2)
        plot(-300:800,pente.*[-300:800]+ordorg,'b:','LineWidth',2)
        xlim([-300 800])
        ylim([-300 800])
        grid on 
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['H ' num2str(s)])
        text(-50,650,message)
    end
%     subplot(1,3,2)
%     plot(H_corr,temp2(i_sim,10),'r+')
%     
    % LE
    figure('Name','LE')
    plot(dates,temp2(:,14),'b-','linewidth',2)
    hold on 
        scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),5,'m','filled')
    hold on
plot(dates4,temp4(:,14),'r-','linewidth',2)
    %plot(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),'k-','linewidth',1)%,'Marker','o','MarkerSize',5)
    scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),5,'k','filled')
    dynamicDateTicks([],[],'dd/mm');
    datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('W.m^{-2}')
    title(['LE ' ])
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
    
    figure('Name',['LE scatter' ])
    for s=1:3
       subplot(1,3,s)
        score=eval(['scoresbLE',num2str(s)]);   
        temp_=eval(['temp',num2str(s+1)]);
         [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
        message=sprintf('pente: %.2f \n ordorg: %.2f \n r2: %.2f \n biais: %.2f \n rmse: %.2f',pente,ordorg,r2,biais,rmse);
         hold on
        plot(LE_corr,temp_(i_sim,14),'r+')
        plot(-100:800,-100:800,'k--','linewidth',2)
        plot(-100:800,pente.*[-100:800]+ordorg,'b:','linewidth',2)
        xlim([-200 800])
        ylim([-200 800])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['LE sim' num2str(s) ])
        grid on
        text(-80,690,message)
    end



%% EVAP/TRANSPI

%partitionnement de l'évapotranspiration
temp_=temp3;
LEf=temp_(:,12);% evap feuille (eau interceptée)
LEtr=temp_(:,13);%transpi
LEs=temp_(:,11);% evap sol
LEtot=temp_(:,14);% total

  figure_handle=figure('Name',['Evap/Transpi' ])
  subplot(3,1,1)
        plot(dates,temp_(:,14),'r-','linewidth',2)
        hold on
        plot(dates,temp_(:,11),'b-','linewidth',2)
        hold on
        plot(dates,temp_(:,13),'color',[0 0.5 0],'linewidth',2)
        hold on 
         plot(dates,temp_(:,12),'g','linewidth',2)
        hold on
        scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),5,'k','filled')
        datetick('x','mmm-yy','keepticks')
        ylabel('W.m^{-2}')
        ylim([-100 750])
        grid on 
        title(['Partition E/T' num2str(s+1) ])
        legend('LE total','LE Evap.','LE Transp.','LE feuille','Observed LE','location','eastoutside')
         dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
 % pluies (journalières pour vision à l'échelle saionnière et semi horaire pour une vision à l'échelle de quelques jours)          
    subplot(3,1,2)
    Precip= Pcumuljour;
    Pcumuljourtime=datenum(Pcumuljour.TIMESTAMP);
    [ind_t_P,ind_t3,ind_t4]=deal (  find(Pcumuljourtime>=dates(1) & Pcumuljourtime<=dates(end)),find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end))  , find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end))   );
     bar(Pcumuljourtime(ind_t_P),Precip.P(ind_t_P))
     hold on 
     bar(obs_FB.TIMESTAMP(ind_t),obs_FB.P(ind_t))% 30min
     grid on 
 dynamicDateTicks(); 
datetick('x','mmm-yy','keepticks','keeplimits')
     legend('Dayly rainfall','30 min rainfall','location','eastoutside')
 ylabel('P (mm')   
       dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

 % pourcentage d'évapotranspiration rempli par l'évaporation du sol et de
 % l'eau interceptée par les feuilles 
       	subplot(3,1,3)

   % enlever les valeurs négatives, regarder juste les % d'évaporation et non pas de condensation     	
   posLEf=LEf(LEtot >LEf & LEtot >0.5 & LEf>=0);
   posLEs=LEs(LEtot >LEf & LEtot >0.5 & LEf>=0);
   posYY=LEtot(LEtot >LEf &LEtot >0.5 & LEf>=0);
posX=dates(LEtot >LEf & LEtot >0.5 & LEf>=0);
       % plot(dates,(temp_(:,11)./temp_(:,14))*100,'b-','linewidth',2)
        hold on
        %plot(dates,(temp_(:,13)./temp_(:,14))*100,'color',[0 0.5 0],'linewidth',2)
        hold on 
       scatter(posX,(posLEf./posYY)*100,19,'g')
        %plot(dates,(temp_(:,12)./temp_(:,14))*100,'g','linewidth',2)
        hold on
        scatter(posX,(posLEs./posYY)*100,19,'b')
        %scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),5,'k','filled')
        datetick('x','mmm-yy','keepticks')
        ylabel('% of LE total')
        %ylim([-100 750])
        grid on 
        legend('Evap feuille','evap sol')
        %legend('LE Evap.','LE Transp.','LE feuille','Observed LE','location','eastoutside')
         dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')
            ylim([-10 150])

    all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
  
    % 
    %   figure_handle=figure('Name',['Evap: sol/interception' ])
    % for s=1:3
    %    	subplot(3,1,s)
    %     temp_=eval(['temp',num2str(s+1)]);
    %      LEf=temp_(:,12);% evap feuille (eau interceptée)
    %     %
    %     LEtot=temp_(:,14);
    % 
    %     plot(dates,LEf,'color',[0 0.5 0],'linewidth',2)
    %     hold on
    %     plot(dates,LEtot,'linewidth',2)
    %     datetick('x','mmm-yy','keepticks')
    %     ylabel('W.m^{-2}')
    %     ylim([-100 750])
    %     grid on 
    %     title(['Partition E/T' num2str(s+1) ])
    %     legend('LEsol','LEtot','interp','location','eastoutside')
    %      dynamicDateTicks()
    %        datetick('x','mmm-yy','keepticks')
    % end
    % all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    % linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots



    
    % pause
    
   % close all
    
    %% Evapotranspiration 
    
    % ETR
    % on a tronqué le debut pour ne pas prendre en compte l'initialisation
    % mais ETR tot fait la somme des etr depuis le début de la periode de
    % simulation, il faut partir de 0 !!  et retirer de la somme les
    % moments ou il n'y a pas de mesures pour pouvoir comparer avec les obs
   ind_LE=isnan(obs_FB.LE_hf(ind_t));
   LEbis2=temp2(:,14);
    [LEbis2,LEbis3,LEbis4]=deal(  temp2(:,14),temp3(:,14),temp4(:,14));
      for i=1:length(obs_FB.LE_hf(ind_t))
        if ind_LE(i) ==1
             LEbis2(i)=NaN;
             LEbis3(i)=NaN;
             LEbis4(i)=NaN;
        end 
    end 
    
    figure('Name','ETR')
    plot(dates,cumsum(LEbis2*dt/Lv,'omitnan'),'b-','linewidth',2)
     hold on 
     plot(dates,cumsum(LEbis3*dt/Lv,'omitnan'),'m-','linewidth',2)
     hold on
     plot(dates,cumsum(LEbis4*dt/Lv,'omitnan'),'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),cumsum(obs_FB.LE_hf(ind_t)*dt/Lv,'omitnan'),'k-','linewidth',2)
    hold on
    %plot(dates,cumsum(temp2(:,14)*dt/Lv,'omitnan'),'r--','linewidth',2)
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('mm')
    title(['ETR '])
    legend('20m ref','karst 20m','karst étalonné','obs','sim non corr''location','eastoutside')
    %legend('Simulated ET','Observed ET','Simulated ET (non corrigée)')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
    
    % ETR journalière
    %N_LE nombre de valeurs qu'il y a. On peut filtrer les jours avec trop
    %de valeurs manquantes si plus de 4 h manquantes par exemple nan  pour
    %la journée pour les obs et les simu aux jours correspondant 
    [ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),'day',@nansum);
    % [ETR_day_obs,days_obs,N_LE_obs]=timeavgCK(obs_FB.TIMESTAMP(ind_t),LE_corr,'day',@nansum);
    
    ETR_day_obs=ETR_day_obs*dt/Lv;
    [ETR_day_sim1,days_sim1]=timeavgCK(dates,temp2(:,14),'day',@nansum);
    [ETR_day_sim2,days_sim2]=timeavgCK(dates3,temp3(:,14),'day',@nansum);
    [ETR_day_sim3,days_sim3]=timeavgCK(dates4,temp4(:,14),'day',@nansum);
    [ETR_day_sim1,ETR_day_sim2,ETR_day_sim3]=deal( ETR_day_sim1*dt/Lv, ETR_day_sim2*dt/Lv, ETR_day_sim3*dt/Lv);
    for i=1:length(ETR_day_obs)
        if N_LE_obs(i) <24 % si inférieur à 40 mesures par jour (20 h/24) pas du cumul journalier
            ETR_day_obs(i)=NaN;
            ETR_day_obs3(i)=NaN;
           
        end 
    end 
    
    figure('Name',['ETR dayly' ])
    plot(days_sim1,ETR_day_sim1,'b-','linewidth',2)
    hold on 
    plot(days_sim3,ETR_day_sim2,'m-','linewidth',2)
    hold on 
    plot(days_sim2,ETR_day_sim3,'r-','linewidth',2) 
    hold on 
    plot(days_obs,ETR_day_obs,'k-','linewidth',2)
    grid on
    ylabel('mm')
        legend(['sim1: ' num2str(sum(ETR_day_sim1,'omitnan')) 'mm'],['sim2: ' num2str(sum(ETR_day_sim2,'omitnan')) 'mm'],['sim3: ' num2str(sum(ETR_day_sim3,'omitnan')) 'mm'],...
    ['obs: ' num2str(sum(ETR_day_obs,'omitnan')) 'mm'],'location','eastoutside');
legend('20 m ref','20 m karst','Karst étalonné','obs')
   dynamicDateTicks([],[],'dd/mm')
   datetick('x','mmm-yy','keepticks')
    %datetick('x')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


    figure('Name',['ETR dayly scatter' ])
    for s=1:3
        ETR_day_sim=eval(['ETR_day_sim',num2str(s)]);
        scoresETR_day=vr_scoresstr(ETR_day_obs,ETR_day_sim);
        [pente,ordorg,r2,rmse,biais]=deal(scoresETR_day.pente,scoresETR_day.ordorg, scoresETR_day.cor2,scoresETR_day.rmse,scoresETR_day.mbias  );
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        subplot(3,1,s)
        plot(ETR_day_obs,ETR_day_sim,'r+')
        ylabel('mm')
        hold on
        plot(-2:7,-2:7,'k--','linewidth',2)
        hold on
        plot(-2:10,pente.*[-2:10]+ordorg,'b:','linewidth',2)
        xlim([-2 6])
        ylim([-2 6])
        xlabel('Obs. mm')
        ylabel('Sim. mm')
        title(['ETR daily scatter ' num2str(s) ])
        text(-1.5,2,message)
    end

   figure('Name','potentiel foliaire')
    %plot(obs_FB.TIMESTAMP(ind_t_12),obs_FB.ALB(ind_t_12)*100,'k-','linewidth',1.5)
    hold on
    plot(dates,temp4(:,18)/101.9,'b-','linewidth',1.5)
    dynamicDateTicks();
   
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('Potentiel hydrique foliaire (Mpa)')

    title([' '])
    %legend('obs','sim L1 ref','sim L1B','sim L1C','location','eastoutside')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 
%% pluies, potentiel hydrique, ET, humidité sol 
figure_handle=figure()

% pluies (journalières pour vision à l'échelle saionnière et semi horaire pour une vision à l'échelle de quelques jours)          
    subplot(4,1,1)
    Precip= Pcumuljour;
    Pcumuljourtime=datenum(Pcumuljour.TIMESTAMP);
    [ind_t_P,ind_t3,ind_t4]=deal (  find(Pcumuljourtime>=dates(1) & Pcumuljourtime<=dates(end)),find(obs_FB.TIMESTAMP>=dates3(1) & obs_FB.TIMESTAMP<=dates3(end))  , find(obs_FB.TIMESTAMP>=dates4(1) & obs_FB.TIMESTAMP<=dates4(end))   );
     bar(Pcumuljourtime(ind_t_P),Precip.P(ind_t_P))
     hold on 
     bar(obs_FB.TIMESTAMP(ind_t),obs_FB.P(ind_t))% 30min
     grid on 
dynamicDateTicks([],[],'dd/mm')
datetick('x','mmm-yy','keepticks','keeplimits')
     legend('Dayly rainfall','30 min rainfall')
 ylabel('P (mm')   
       dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


 %potentiel hydrique 
    subplot(4,1,2) 
    plot(dates,temp4(:,18)/101.9,'b-','linewidth',1.5)
    dynamicDateTicks([],[],'dd/mm')
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('Potentiel hydrique foliaire (Mpa)')
    title([' '])
    legend('Hydric potential')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 

% ETR 

LEf=temp_(:,12);% evap feuille (eau interceptée)
LEtr=temp_(:,13);%transpi
LEs=temp_(:,11);% evap sol
LEtot=temp_(:,14);% total
  subplot(4,1,3)
        plot(dates,LEtot,'r-','linewidth',2)
        hold on
        plot(dates,LEs,'b-','linewidth',2)
        hold on
        plot(dates,LEtr,'color',[0 0.5 0],'linewidth',2)
        hold on 
         plot(dates,LEf,'g','linewidth',2)
        hold on
        scatter(obs_FB.TIMESTAMP(ind_t),obs_FB.LE_hf(ind_t),5,'k','filled')
        datetick('x','mmm-yy','keepticks')
        ylabel('W.m^{-2}')
        ylim([-100 750])
        grid on 
        title(['Partition ETR' num2str(s+1) ])
        legend('LE total','LE Evap.','LE Transp.','LE feuille','Observed LE')
         dynamicDateTicks([],[],'dd/mm')
           datetick('x','mmm-yy','keepticks')


   subplot(4,1,4)        

      plot(dates,temp2(:,45)/100,'b-','linewidth',2)
    
    hold on
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_5_1(ind_t),'k-','linewidth',1)
hold on
    plot(dates,temp2(:,49)/100,'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_50_2(ind_t)-0.2,'k--','linewidth',1)
    hold on 
    grid on 
    dynamicDateTicks([],[],'dd/mm')
    datetick('x','mmm-yy','keepticks','keeplimits')
    ylabel('/theta (m^3.m^{-3})')
    legend('5 cm sim','5 cm obs','50 cm sim','50 cm obs')


 all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
    linkaxes( all_ha, 'x' ); % pour faire décaler l'axe des x en même temps sur tous les subplots
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


%% REPARTITION PAR PERIODE DES COMPOSANTES DE L'ET
TTobs = retime(obs_FBt(ind_t,:),'regular','sum','TimeStep',calmonths(1));

TTobsET=array2table([days_obs ETR_day_obs]);
TTobsET.Var1=datetime(TTobsET.Var1,'ConvertFrom','datenum');
TTobsET= table2timetable(TTobsET);
TTobsET=retime(TTobsET,'regular','sum','TimeStep',calmonths(1));

 ind_LEs=isnan(obs_FB.LE_sf(ind_t));
   ind_LEh=isnan(obs_FB.LE_hf(ind_t));
   ind_LE=isnan(obs_FB.LE(ind_t));
   
TTsim1=1;
TTsim2=1;
TTsim3=1;

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
if i==1
    TTsim1=TTsim;
elseif i==2
    TTsim2=TTsim;
elseif i==3
    TTsim3=TTsim;
end 


figure_handle=figure('Name',['ET' num2str(i) ])
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

ETobs_sim(:,i)=TTsim.Var5-TTobsET.Var2;

end


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
%ETobs_sim=[(TTsim.Var1, ETobs_sim]
figure()
yy=[ETobs_sim(:,1) ETobs_sim(:,2) ETobs_sim(:,3)]
bar(TTsim.Var1,yy)
grid on
ylabel('Surplus ou déficit ET simulé (mm)')
legend('sim1: 20m ref','sim2: 20m karst' ,'sim3: karst étalonné')







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


 %%   Bilan d'eau + fermeture  en cumulé
figure('Name',['Bilan'])
    for s=2:4
        t=eval(['temp' num2str(s) ]);
         subplot(1,3,s-1)
        plot(dates,t(:,29)-t(1,29),'k-','linewidth',2) % précipitations
        hold on
        plot(dates,t(:,23)-t(1,23),'r-','linewidth',2)% ET
        plot(dates,t(:,27)-t(1,27),'b-','linewidth',2) % ruissellement
        plot(dates,t(:,28)-t(1,28),'b--','linewidth',2)% Drainage
        plot(dates,t(:,30)-t(1,30),'g-','linewidth',2)% stock
        plot(dates,t(:,29)-t(:,23)-t(:,27)-t(:,28)-t(:,30)+t(1,30),'m-','linewidth',2)% bilan
        dynamicDateTicks(); datetick('x','mmm-yy','keepticks')
        set(gca, 'YGrid', 'on', 'XGrid', 'off')
        ylabel('mm')
        title(['Bilan de masse sim' num2str(s-1)])
        legend('Rain','ET','Rn-off','Drainage','/Delta Stock','Balance')
    end  
%%  Variation stock   + Bilan obs  (partie pas au point)

    % Variation du stock d'eau dans le sol 0-250cm
    % Calcul d'un stock d'eau interpolé à partir des observations
    prof=0:200;
    prof_obs=[5 10 50 100 150 200];
    hum_interp=NaN(length(ind_t),length(prof));
    for n=1:length(ind_t)
        Y=[obs_FB.SWC_TD_5_1(ind_t(n)) obs_FB.SWC_TD_50_1(ind_t(n))];
          Y(Y<=0)=NaN;
          ind_Y=find(~isnan(Y));
          if length(ind_Y)==1
              Y(ind_Y+1)=Y(ind_Y);
              hum_interp(n,:)=naninterp1(prof_obs,Y,prof,'linear');
          elseif length(ind_Y)>1
              hum_interp(n,:)=naninterp1(prof_obs,Y,prof,'linear');
          end
    end

    stock_obs=sum(hum_interp,2)*10;
    stock_obs=stock_obs-stock_obs(1);
    stock_obs(stock_obs<=-10)=NaN;
    
    scoresstock.(['sim_' num2str(i)])=vr_scoresstr(stock_obs,temp3(:,30)-temp3(1,30));
    pente=scoresstock.(['sim_' num2str(i)]).pente;
    ordorg=scoresstock.(['sim_' num2str(i)]).ordorg;
    r2=scoresstock.(['sim_' num2str(i)]).cor2;
    rmse=scoresstock.(['sim_' num2str(i)]).rmse;
    biais=scoresstock.(['sim_' num2str(i)]).mbias;
    message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
    
    figure('Name','dStock')
    plot(obs_FB.TIMESTAMP(ind_t),stock_obs,'k-','linewidth',2)
    hold on
    plot(dates,temp2(:,30)-temp2(1,30),'r-','linewidth',2)
    ylabel('mm')
    title(['Var. Stock ' ])
    legend('obs','sim','location','eastoutside')
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    text(dates(end-48*100),max(stock_obs-20),message)
    
    % Bilan Obs.
    % Fermeture du bilan
    ETR_obs=cumsum(obs_FB.LE(ind_t)*dt/Lv,'omitnan');
    res_bilan=temp2(:,29)-ETR_obs-stock_obs;
    
    figure('Name','Bilan Obs.')
    plot(dates,temp2(:,29),'k-','linewidth',2)
    hold on
    plot(obs_FB.TIMESTAMP(ind_t),ETR_obs,'r-','linewidth',2)
    plot(obs_FB.TIMESTAMP(ind_t),stock_obs,'g-','linewidth',2)
    plot(obs_FB.TIMESTAMP(ind_t),res_bilan,'b-','linewidth',2)
    ylabel('mm')
    title(['Bilan Obs. ' ])
    legend('pluie','ETR','dSTOCK','Res. bilan','location','eastoutside')
    dynamicDateTicks(); datetick('x','mmm-yy','keepticks')


    
%% 

    %% Profils de température du sol (chroniques et scatter plot+stats)
    
    figure('Name','Tsol 5')
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,36),'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,36),'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,36),'b-','linewidth',2)
    ylim([0 35])
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('°C')
    title(['TSOL 5 ' ])
    grid on
    legend('obs TD1','sim ref','sim2','sim3','location','eastoutside')
    
   figure('Name','Tsol 5 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresT5.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.TS_TD_5_1(i_obs),temp(i_sim,36));
        subplot(1,3,s)
        score=eval(['scoresT5.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.TS_TD_5_1(i_obs),temp(:,36),'r+')
        plot(-2:27,-2:27,'k--','LineWidth',2)
        plot(-2:27,pente.*[-2:27]+ordorg,'b:','LineWidth',2)
        xlim([-2 27])
        ylim([-2 27])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 5cm ' num2str(s)])
        text(0,22,message)    
    end
%     
    figure('Name','Tsol 15')
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_15_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_15_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,37),'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,37),'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,37),'b-','linewidth',2)
    ylim([0 35])
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('°C')
    title(['TSOL 10 ' ])
    legend('obs TD1','sim ref','sim2','sim3','location','eastoutside')
    
    
       figure('Name','Tsol 15 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresT15.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.TS_TD_15_1(i_obs),temp(i_sim,37));
        subplot(1,3,s)
        score=eval(['scoresT15.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.TS_TD_15_1(i_obs),temp(:,37),'r+')
        plot(2:25,2:25,'k--','LineWidth',2)
        plot(2:25,pente.*[2:25]+ordorg,'b:','LineWidth',2)
        xlim([2 25])
        ylim([2 25])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 15cm ' num2str(s)])
        text(2.2,21,message)    
    end
    
    figure('Name','Tsol 25')
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_25_1(ind_t),'k-','linewidth',2)
    hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_25_2(ind_t),'k--','linewidth',2)
%     hold on
    plot(dates,temp2(:,38),'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,38),'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,38),'b-','linewidth',2)
    ylim([0 35])
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('°C')
    title(['TSOL 25 ' ])
    legend('obs TD1','sim ref','sim2','sim3','location','eastoutside')
    
    
          figure('Name','Tsol 25 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresT25.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.TS_TD_25_1(i_obs),temp(i_sim,38));
        subplot(1,3,s)
        score=eval(['scoresT25.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.TS_TD_25_1(i_obs),temp(:,38),'r+')
        plot(2:25,2:25,'k--','LineWidth',2)
        plot(2:25,pente.*[2:25]+ordorg,'b:','LineWidth',2)
        xlim([2 25])
        ylim([2 25])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 25cm ' num2str(s)])
        text(2.2,21,message)    
    end
    
    figure('Name','Tsol 35')
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_35_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_35_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,39),'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,39),'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,39),'b-','linewidth',2)
    ylim([0 35])
    datetick('x','mmm-yy')
    grid on 
    ylabel('°C')
    title(['TSOL 35 ' ])
    legend('obs TD1','sim ref','sim2','sim3','location','eastoutside')
%     
          figure('Name','Tsol 35 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresT35.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.TS_TD_35_1(i_obs),temp(i_sim,39));
        subplot(1,3,s)
        score=eval(['scoresT35.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.TS_TD_35_1(i_obs),temp(:,39),'r+')
        plot(2:25,2:25,'k--','LineWidth',2)
        plot(2:25,pente.*[2:25]+ordorg,'b:','LineWidth',2)
        xlim([2 25])
        ylim([2 25])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 35cm ' num2str(s)])
        text(2.2,21,message)    
    end

    figure('Name','Tsol 50')
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_50_1(ind_t),'k-','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_50_2(ind_t),'k--','linewidth',2)
   hold on
    plot(dates,temp2(:,40),'r-','linewidth',2)
    hold on
    plot(dates3,temp3(:,40),'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,40),'b-','linewidth',2)
    ylim([0 35])
    dynamicDateTicks(); datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('°C')
    title(['TSOL 50 '])
    legend('obs TD1','sim ref','sim2','sim3','location','eastoutside')
    
    
     figure('Name','Tsol 50')
         plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_1(ind_t),'k--','linewidth',2)
%     hold on
%     plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_5_2(ind_t),'k--','linewidth',2)
    hold on
    plot(dates,temp2(:,36),'r--','linewidth',2)
     hold on
    plot(dates4,temp4(:,36),'b--','linewidth',2)
    hold on
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.TS_TD_50_1(ind_t),'k-','linewidth',2)
   hold on
    plot(dates,temp2(:,40),'r-','linewidth',2)
    hold on
    plot(dates4,temp4(:,40),'b-','linewidth',2)
    hold on 
  
    ylim([0 35])
    dynamicDateTicks(); datetick('x','mmm-yy','keepticks')
    grid on
    ylabel('°C')
    title(['TSOL 50 '])
    legend('obs 5cm','ref 5cm','smaller Ktherm 5cm','obs 50cm','ref 50cm','smaller Ktherm 50cm','location','east')
       
%     
            figure('Name','Tsol 50 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresT50.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.TS_TD_50_1(i_obs),temp(i_sim,40));
        subplot(1,3,s)
        score=eval(['scoresT50.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        plot(obs_FB.TS_TD_50_1(i_obs),temp(:,40),'r+')
        plot(2:25,2:25,'k--','LineWidth',2)
        plot(2:25,pente.*[2:25]+ordorg,'b:','LineWidth',2)
        xlim([2 25])
        ylim([2 25])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 50cm ' num2str(s)])
        text(2.2,21,message)    
    end

%close all
%% Profils de teneur en eau dans le sol 

    figure('Name','SMC 5')
  
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
    ylabel('/theta (m^3.m^{-3})')
    title(['SMC 5 ' ])
     legend('sim ref','sim2','sim3','obs','location','eastoutside')
   legend('/theta r min /theta s min', '/theta r max /theta s min','/theta r max /theta s max')
   %ax=gca;
%ax.FontSize = 10; 
%ax.FontWeight = 'bold';  
%   
    figure('Name','SMC 15')
    hold on
    plot(dates,temp2(:,46)/100,'b-','linewidth',2)
    hold on
    plot(dates3,temp3(:,46)/100,'m-','linewidth',2)
    hold on
   plot(dates4,temp4(:,46)/100,'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_15_1(ind_t),'k-','linewidth',1)
    %ylim([0 0.14])
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    grid on 
    ylabel('/theta (m^3.m^{-3})')
    title(['15 cm' ])
     legend('sim ref','sim2','sim3','obs','location','eastoutside')
%     legend('/beta min','/beta max','obs')
%    %ax=gca;
% ax.FontSize = 10; 
% ax.FontWeight = 'bold'; 

 figure('Name','SMC 25')
%     hold on
%     plot(dates,temp(:,47)/100,'r-','linewidth',2)
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
    ylabel('/theta (m^3.m^{-3})')
    title(['25 cm' ])
    legend('sim ref','sim2','sim3','obs','location','east')
    legend('hg min','hg max','obs')
    %legend('/theta r min /theta s min', '/theta r max /theta s min','/theta r max /theta s max', 'obs')
%     ax=gca;
% ax.FontSize = 10; 
% ax.FontWeight = 'bold'; 

figure('Name','SMC 35')
%     hold on
%     plot(dates,temp(:,48)/100,'r-','linewidth',2)
    hold on
    plot(dates,temp2(:,48)/100,'b-','linewidth',2)
    hold on
    plot(dates3,temp3(:,48)/100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,48)/100,'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_35_2(ind_t),'k-','linewidth',1)
    %ylim([0 0.14])
    dynamicDateTicks();
    datetick('x','mmm-yy','keepticks')
    ylabel('/theta (m^3.m^{-3})')
    grid on 
    title(['SMC 35 ' ])
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
%     ax=gca;
% ax.FontSize = 10; 
% ax.FontWeight = 'bold'; 
%     
    
    figure('Name','SMC 50')
    %     hold on
%     plot(dates,temp(:,49)/100,'r-','linewidth',2)
    hold on
    plot(dates,temp2(:,49)/100,'b-','linewidth',2)
    hold on
    plot(dates3,temp3(:,49)/100,'m-','linewidth',2)
    hold on
    plot(dates4,temp4(:,49)/100,'r-','linewidth',2)
    hold on 
    plot(obs_FB.TIMESTAMP(ind_t),obs_FB.SWC_TD_50_2(ind_t)-0.2,'k-','linewidth',1)
    hold on 
%     yyaxis right
%     bar(obs_FB.TIMESTAMP(ind_t),obs_FB.P(ind_t))
% set(h1, 'Ydir', 'reverse')
% set(h1, 'YAxisLocation', 'Right')
%     %ylim([0 0.14])
    grid on 
    dynamicDateTicks(); 
    datetick('x','mmm-yy','keepticks','keeplimits')
    ylabel('/theta (m^3.m^{-3})')
    title(['SMC 50 '])
    legend('sim ref','sim2','sim3','obs','location','eastoutside')
%     ax=gca;
% ax.FontSize = 10; 
% ax.FontWeight = 'bold'; 

  figure('Name','Tsol 50 scatter')
    for s=1:3
        temp=eval(['temp',num2str(s+1)]);
        scoresW50.(['sim_' num2str(s)])=vr_scoresstr(obs_FB.SWC_TD_50_1(i_obs),temp(i_sim,40));
        subplot(1,3,s)
        score=eval(['scoresW50.sim_' num2str(s)]);
        pente=score.pente;
        ordorg=score.ordorg;
        r2=score.cor2;
        rmse=score.rmse;
        biais=score.mbias;
        message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
        hold on
        %plot(obs_FB.SWC_TD_50_1(i_obs),temp(:,49)/100,'r+')
         plot(0:0.7,0:0.7,'k--','LineWidth',2)
%         plot(1:1,pente.*[0:1]+ordorg,'b:','LineWidth',2)
% %         xlim([2 25])
% %         ylim([2 25])
        grid on 
        xlabel('Obs. °C')
        ylabel('Sim. °C')
        title(['Temp 50cm ' num2str(s)])
        text(2.2,2,message)    
    end


%% Scatter plot bilan d'energie (avec saisons)

% filtrer valeurs journalières seulement avec SWin
ind_day=find(obs_FBt.SW_IN(i_obs) > 0.5 );
temp__day=temp2(ind_day,:);
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

figure_handle=figure('Name',['Scatter energy balance karst version'])
 subplot(2,2,1)
score=vr_scoresstr(obs_.NETRAD ,temp__day(:,07));
score_RN=score;
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
 message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
hold on
    plot(obs_.NETRAD,temp__day(:,07),'b+')
      hold on
   %plot(obs_FB.NETRAD(i_obs(ind_s_day)),temp4(ind_s_day,7),'r+')%    plot(obs_FB.NETRAD(ind_w1),temp_(ind_w1,7),'b+')
        plot(0:1000,0:1000,'k--','LineWidth',2)
        plot(0:1000,pente.*[0:1000]+ordorg,'k:','LineWidth',2)
        xlim([-100 1000])
        ylim([-100 1000])
        grid on 
        xlabel('Obs. (W.m^{-2})')
        ylabel('Sim. (W.m^{-2})')
       %legend('All year','Summer','','Regression')
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
  %plot(obs_FB.LE_hf(i_obs(ind_s_day)),temp4(ind_s_day,14),'r+')
       hold on 
        plot(-100:800,-100:800,'k--','linewidth',2)
        plot(-100:800,pente.*[-100:800]+ordorg,'k:','linewidth',2)
        xlim([-100 600])
        ylim([-100 600])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['LE ' ])
   %legend('All year','Summer','','Regression')
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
  %plot(H_corr(ind_s_day),temp4(ind_s_day,10),'r+')
       hold on 
       plot(-300:800,-300:800,'k--','LineWidth',2)
        plot(-300:800,pente.*[-300:800]+ordorg,'b:','LineWidth',2)
        xlim([-300 800])
        ylim([-300 800])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['H ' ])
        %legend('All year','Summer','','Regression')
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
  %plot(G0_obs(ind_s_day),temp4(ind_s_day,15),'r+')
       hold on 
        plot(-120:120,-120:120,'k--','LineWidth',2)
        plot(-120:120,pente.*[-120:120]+ordorg,'k:','LineWidth',2)
        xlim([-100 120])
        ylim([-100 120])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['G ' ])
     %legend('All year','Summer','','Regression')
        grid on 
        text(-110,60,message,"BackgroundColor",[1 1 1],'FontSize',12)
    ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


%%  scatter saisons
% %indexation par mois  pour selectionner les saisons 
month_d=month(obs_FBt.TIMESTAMP(i_obs));
ind_w1=find(month_d==12);
ind_w2=find(month_d<6);
%ind_sp=find(month_d >=3 & month_d<6);
ind_s=find(month_d >=6 & month_d<9);
ind_a=find(month_d >=9 & month_d<12);

%à faire pour les journées et nuits 

figure_handle=figure('Name',['Scatter energy balance'])

 subplot(2,2,1)
 s=1
score=eval(['scoresRN',num2str(s)]);
  tempb=eval(['temp',num2str(s+1)]);       
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
 message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
hold on
    plot(obs_FB.NETRAD(i_obs),tempb(i_sim,07),'b+')
%    plot(obs_FB.NETRAD(ind_w1),temp_(ind_w1,7),'b+')
%    hold on 
%    plot(obs_FB.NETRAD(ind_w2),temp_(ind_w2,7),'b+')
%     hold on 
% %    plot(obs_FB.LE_sf(ind_sp),temp_(ind_sp,14),'b+')
% %     hold on   
%      plot(obs_FB.NETRAD(ind_a),temp_(ind_a,7),'b+')
      hold on
   plot(obs_FB.NETRAD(i_obs(ind_s)),temp_(ind_s,7),'r+')
        plot(0:1200,0:1200,'k--')
        plot(0:1200,pente.*[0:1200]+ordorg,'b:')
        xlim([-100 1200])
        ylim([-100 1200])
        grid on 
        xlabel('Obs. (W.m^{-2})')
        ylabel('Sim. (W.m^{-2})')
         legend('All year','','','Summer')
        title(['R_{NET} ' num2str(s)])
        text(50,1000,message)

subplot (2,2,3)%
s=1
 score=eval(['scoresbLE',num2str(s)]);   
 temp_=eval(['temp',num2str(s+1)]);
 [pente,ordorg,r2,rmse,biais]=deal(score.pente,score.ordorg, score.cor2,score.rmse,score.mbias  );
message=sprintf('pente: %.2f\nordorg: %.2f\nr2: %.2f\nbiais: %.2f\nrmse: %.2f',pente,ordorg,r2,biais,rmse);
   plot(obs_FB.LE_hf(ind_w1),temp_(ind_w1,14),'b+')
   hold on 
   plot(obs_FB.LE_hf(ind_w2),temp_(ind_w2,14),'b+')
    hold on 
%    plot(obs_FB.LE_hf(ind_sp),temp_(ind_sp,14),'b+')
%     hold on   
     plot(obs_FB.LE_hf(ind_a),temp_(ind_a,14),'b+')
     hold on
   plot(LE_corr,temp_(ind_s,14),'r+')
       hold on 
        plot(-100:800,-100:800,'k--','linewidth',2)
        plot(-100:800,pente.*[-100:800]+ordorg,'b:','linewidth',2)
        xlim([-100 800])
        ylim([-100 800])
        xlabel('Obs. W.m^{-2}')
        ylabel('Sim. W.m^{-2}')
        title(['LE ' ])
        legend('All year','','','Summer')
        grid on 
        text(-80,690,message)

   
