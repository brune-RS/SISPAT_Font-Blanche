%% Programme de mise en place du profil initial d'humidite du sol

%  Auteur : Cecile VELLUET - 17 Mars 2010
% modifié: février 2023 Brune
% pour données de FONT-BLANCHE

clear
clc
close all
disp('Conditions initiales : Profil initial de temperature et d''humidite du sol - Format SiSPAT')



%% Parametres

%choisir dates de initialisation (1er août 2017 à FB / 2014)

annee=2014;
mois_ci=8;
jour_ci=1;
heure_ci=0;
minute_ci=0;
date_ci=datenum(annee,mois_ci,jour_ci,heure_ci,minute_ci,0);

annee_str=num2str(annee);
mois_str=num2str(mois_ci,'%02d');
jour_str=num2str(jour_ci,'%02d');

% observations humidité dans le sol (m3/m3) à Font-Blanche
path_lv='C:\SISPAT\tableaux matlab\Font-Blanche\';
file_obs_l2=[path_lv 'FB_SoilMoistureL2A_TD1_TD2.mat'];
file_obs_h=[path_lv 'FB_SoilMoistureL1_TD1_TD2.mat'];
file_obs_T=[path_lv 'FB_SoilTemp.mat'];
FB_SoilMoistureL1_TD1_TD2=importdata(file_obs_h);
load(file_obs_T)
load(file_obs_l2)

% Profondeur des observations (cm) à FB 
prof_obs = [5,15,25,35,50]; 

% Choix d'initialisation
% iWr=1, les simulations sont intilialisees aux observations
% iWr=0, W(sim,i)=W(obs,i)-W(obs,min)+Wr(sim)
iWr=1;
         
% Maillage de FB   profondeur des noeuds 
path_EN='C:\SISPAT\execution\entrees\';
%file_Maillage=[path_EN 'prof_maillageFB.csv']; %maillage à 20m
file_Maillage=[path_EN 'profmaillageFB2m.mat']; %maillage à 2m 
load(file_Maillage);
%profmaillageP1=prof_maillageFB;  %déja en cm
profmaillageP1=profmaillageFB2m

% Maillage repris de celui de Puechabon
% path_EN='C:\SISPAT\execution\entrees\';
% file_Maillage=[path_EN 'prof_maillageP.csv'];
% profmaillageP1=importdata(file_Maillage);
% profmaillageP1=100*profmaillageP1;  % m -------> cm
               
%% Lecture des donnees de temperature et humidite
FB_SoilMoistureL1_TD1_TD2 = standardizeMissing(FB_SoilMoistureL1_TD1_TD2, -9999);
moist=FB_SoilMoistureL1_TD1_TD2;
name=moist.Properties.VariableNames;
% for i=2:13
%     d=name{i};
%     moist{:,i}=moist{:,i}*1;  % Prise en compte des 80% d'éléments grossiers 
% end 

moistb=removevars(moist,{'TIMESTAMP'});
%moistb=array2table(table2array(moistb)*0.2);
FBnSoilTemperatureL1Copie = standardizeMissing(FBnSoilTemperatureL1Copie, -9999);
temp=FBnSoilTemperatureL1Copie;

%% Récupérer profil à la date choisie 
% Recherche du pas de temps le plus proche   

% 2017i
%ci=find((abs(datenum(moist.TIMESTAMP)-date_ci))==min(abs(datenum(moist.TIMESTAMP)-date_ci)));
%j_ci=find((abs(datenum(temp.TIMESTAMP)-date_ci))==min(abs(datenum(temp.TIMESTAMP)-date_ci)));
%i_ci=find(moist.date_time == date_ci);
%j_ci=find(temp.TIMESTAMP == date_ci);
% TD52=moist.SWC_TD_5_2(i_ci);
% TD51=moist.SWC_TD_5_1(i_ci);
% TD151=moist.SWC_TD_15_2(i_ci);
% TD152=moist.SWC_TD_15_2(i_ci);
% moist_FB=[(TD51+TD52)/2;TD151;moist.SWC_TD_25_2(i_ci);moist.SWC_TD_35_2(i_ci);moist.SWC_TD_50_2(i_ci)];
% TSOL_FB=[temp.TS_TD_5_2(j_ci);temp.TS_TD_15_2(j_ci);temp.TS_TD_25_2(j_ci);
%      temp.TS_TD_35_2(j_ci);temp.TS_TD_50_2(j_ci)];
 

%2014
TD52=moist.SWC_TD_5_3(i_ci);
TD51=moist.SWC_TD_5_4(i_ci);
TD151=moist.SWC_TD_15_1(i_ci);
TD152=moist.SWC_TD_15_1(i_ci);
moist_FB=[(TD51+TD52)/2;TD151;moist.SWC_TD_25_1(i_ci);moist.SWC_TD_35_1(i_ci);moist.SWC_TD_50_1(i_ci)];

TSOL_FB=[temp.TS_TD_5_1(j_ci);temp.TS_TD_15_1(j_ci);temp.TS_TD_25_1(j_ci);
   temp.TS_TD_35_1(j_ci);temp.TS_TD_50_1(j_ci)];

%% Stockage des humidités minimales pour déterminer l'humidité résiduelle à chaque prof ou il y a des obs

Wr_FB=[min(moist.SWC_TD_5_1(moist.SWC_TD_5_1>0));min(moist.SWC_TD_15_1(moist.SWC_TD_15_1>0));
    min(moist.SWC_TD_25_2(moist.SWC_TD_25_2>0));min(moist.SWC_TD_35_2(moist.SWC_TD_35_2>0));
    min(moist.SWC_TD_50_2(moist.SWC_TD_50_2>0))];

%% Ajout de la profondeur 2000 cm pour éviter une extrapolation peu réaliste:


% On fait une première inteprolation linéaire jusqu'à 250cm pour supprimmer
% les NaN
moist_FB=interp1(prof_obs,moist_FB,prof_obs,'linear');
TSOL_FB=interp1(prof_obs,TSOL_FB,prof_obs,'linear');

% on applique les valeurs à 250cm à la profondeur à 2000cm
moist_FB=[moist_FB';moist_FB(end)];
Wr_FB=[Wr_FB;Wr_FB(end)];
TSOL_FB=[TSOL_FB';13.3]; % température moyenne annuelle à 50 cm en condition limite du fond

% Ajout de la profondeur 2000cm /200cm a prof_obs fichier maillage Puech
prof_obs2=[prof_obs 200];

%% Interpolation et extraction des humidités résiduelles 

Wr_FB_i=interp1(prof_obs2,Wr_FB,profmaillageP1,'linear');

% Min des horizons du modèle (0-1,1-20,20-70,70-120 et 120-400cm)  (pas utilisé )
% index == index des noeuds definissant chaque horizon ??  à changer ??
%Wr_FB_H=[min(Wr_FB_i(1:11)) min(Wr_FB_i(12:47)) min(Wr_FB_i(48:91)) min(Wr_FB_i(92:133)) min(Wr_FB_i(134:end))];

%Wr_FB_H=[min(Wr_FB_i(1:11)) min(Wr_FB_i(12:47)) min(Wr_FB_i(48:91)) min(Wr_FB_i(92:133)) min(Wr_FB_i(134:end))];

%% Check Graphique

% Humidite Nord
id=figure
plot(moist_FB,prof_obs2,'k-','linewidth',2)
xlabel('m^3.m^{-3}')
ylabel('Prof. (cm)')

% Temperature Nord
id=figure
plot(TSOL_FB,prof_obs2,'k-','linewidth',2)
xlabel('°C')
ylabel('Prof. (cm)')


%% CREATION DU FICHIER CI 

disp('Nord')
% Preparation des fichiers ci_N.dat (Conditions Initiales)
% - Choix des conditions initiales        (Ligne 1)

disp(' 1. Choix des conditions initiales')
disp ('Information : Se reporter au manuel de l''utilisateur pour definir les conditions initiales')
prompt  = {'Condition a  la limite superieure : itypsup','Condition a  la limite inferieure : itypinf','Choix de condition initiale : ichoixci'}; % Question a  l'utilisateur
title   = 'Conditions initiales - Nord';
lines   = 1;
def     = {'2','3','1'};
answer  = inputdlg(prompt,title,lines,def);
ci      = str2num(char(answer))';

% Condition e la limite superieure
disp('    + Condition a  la limite superieure :')

switch ci(1)
    
    case {1}
        disp ('     h et T imposes au noeud 1')
        prompt  = {'h impose au noeud 1 (m)','T impose au noeud 1 (eC)'};
        title   = 'Conditions e la limite superieure';
        lines   = 1;
        def     = {'Potentiel à  partir de l''humidite de surface initiale + Equation de VG','Temperature de surface initiale'};
        answer  = inputdlg(prompt,title,lines,def);
        itypsup = str2num(char(answer))';
        
    case {2}
        disp ('     flux de masse et de chaleur imposes au noeud 1')
        
    case {3}
        disp ('     flux de masse et T imposes au noeud 1')
        prompt  = {'T imposee au noeud 1 (eC)'};
        title   = 'Conditions a  la limite superieure';
        lines   = 1;
        def     = {'Temperature de surface initiale'};
        answer  = inputdlg(prompt,title,lines,def);
        itypsup = [str2num(char(answer))'];
        
end

% Condition a  la limite inferieure
disp('    + Condition a  la limite inferieure :')

switch ci(2)
    
    case {1}
        disp ('      h et T imposes au noeud n')
        disp('definir h(n)')
        itypinf = [h(n) TSOL_FB(end)];
        
    case {2}
        disp ('      flux de masse et T imposes au noeud n')
        disp('definir fluxn')
        itypinf = [fluxn TSOL_FB(end)];
        
    case {3}
        disp ('      flux de masse gravitaire et T imposee au noeud n')
        itypinf = TSOL_FB(end);
        
    case {4}
        disp ('      T et h imposes au fond mais variables dans le temps')
        
end

% - Preparation des temperatures
disp(' 2. Profil des temperatures initiales')

% Interpolation des donnees Tsoil sur toute la hauteur de sol
% Pour la surface, on extrapole de manière linéaire

TSOL_moist_FBci=interp1(prof_obs2,TSOL_FB,profmaillageP1,'linear');
% Pour 0 à 5cm, on prends la température à 5cm
TSOL_moist_FBci(isnan(TSOL_moist_FBci))=TSOL_FB(1);

% - Preparation des humidites
disp(' 3. Profil des humidites initiales')

switch ci(3)
    
    case {0}
        disp ('    + ichoixci : potentiels')
        disp ('Remarque : Preparer le fichier des potentiels observes')
        
    case {1}
        disp ('    + ichoixci : humidites')
        % Interpolation des donnees SMC sur toute la hauteur de sol
        % Pour la surface, on prend la valeur à 10cm
        
        moist_FB_ci=interp1(prof_obs2,moist_FB,profmaillageP1,'linear');
        % Pour 0 à 5cm, on prends l'humidité à 5cm
        moist_FB_ci(isnan(moist_FB_ci))=moist_FB(1);
        % Pour éviter des rééquilibrages en début de simulation au sein de l'horizon 5, on
        % l'initialise l'humidité à l'humidité résiduelle
        moist_FB_ci(118:end)=Wr_FB(end);
        
        
end


%% Pourcentage de cailloux FB 15,35,50,60,70,80,85,90,90% correction 
%15,50,60,65,70,80,85,90,90%

moist_FB_ci(1:15)=moist_FB_ci(1:15)*0.15;
moist_FB_ci(16:44)=moist_FB_ci(16:44)*0.5;
moist_FB_ci(45:64)=moist_FB_ci(45:64)*0.4;
moist_FB_ci(65:79)=moist_FB_ci(65:79)*0.35;
moist_FB_ci(80:94)=moist_FB_ci(80:94)*0.3;
moist_FB_ci(95:117)=moist_FB_ci(95:117)*0.2;
moist_FB_ci(118:145)=moist_FB_ci(118:145)*0.15;
moist_FB_ci(146:end)=moist_FB_ci(146:end)*0.1;
%% 

% Nombre de noeuds :
nb_noeuds=length(profmaillageP1)

% Enregistrement des fichiers ci_N.dat

disp(' 4. Enregistrement des fichiers de conditions initiales')

disp('    + Enregistrement des fichiers de conditions initiales')

% Enregistrement du fichier ci.dat
file_jachere_ci = ['ci_' annee_str mois_str jour_str 'FB_2m.dat'];
fid = fopen(file_jachere_ci,'w+');
switch ci(1)
    case {1}
        fprintf(fid,'%1.0f %1.0f %1.0f \n',ci);
        fprintf(fid,'%4.2f %4.2f \n',itypsup);
    case {2}
        fprintf(fid,'%1.0f %1.0f %1.0f \n',ci);
    case {3}
        fprintf(fid,'%1.0f %1.0f %1.0f \n',ci);
        fprintf(fid,'%4.2f \n',itypsup);
end

switch ci(2)
    case {1}
        fprintf(fid,'%4.2f %4.2f \n',itypinf);
    case {2}
        fprintf(fid,'%4.2f %4.2f \n',itypinf);
    case {3}
        fprintf(fid,'%4.2f \n',itypinf);
end

for i=1:length(TSOL_moist_FBci)
    fprintf(fid,'%4.2f \n',TSOL_moist_FBci(i));
end

for i=1:length(moist_FB_ci)
    fprintf(fid,'%4.4f \n',moist_FB_ci(i));
end

fclose(fid);