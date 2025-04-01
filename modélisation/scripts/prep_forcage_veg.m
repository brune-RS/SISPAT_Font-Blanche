%% prep_forcage_veg.m

clear all
close all
clc

addpath('C:\Users\brune\Desktop\stage\Stage 2023\SISPAT\')

%% Parametres

annee=2019;
annee_str=num2str(annee);
datedeb=datenum(annee,1,1);
dateend=datenum(annee,12,31);
dates_j=datedeb:dateend;

path_veg=['C:\Users\Etchanchu Jordi\Documents\Niger\Obs\fesa\Data_vegetation\Campagne_' annee_str];
path_res='C:\Users\Etchanchu Jordi\Documents\Niger\Obs\lvl3_SiSPAT';
load('C:\Users\brune\Desktop\stage\Stage 2023\SISPAT\tableaux matlab\laiN0105201701052018.mat');
fic_LAI=laiN0105201701052018
fic_hauteur=[path_veg '\Hauteur_' num2str(annee) '.xlsx'];

%% Lecture des fichiers LAI et hauteur

% LAI

Data_I=fic_LAI;
%Data_J=readtable(fic_LAI,'Sheet','Parcelle J');

dates_LAI=Data_I.Date;
dates_LAI=datenum(datestr(dates_LAI,'dd/mm/yyyy'),'dd/mm/yyyy');
dmin_LAI=min(dates_LAI);

LAI_I=cellfun(@str2num,Data_I.LAI);
LAI_J=cellfun(@str2num,Data_J.LAI);

% Hauteurs

Data_I=readtable(fic_hauteur,'Sheet','Parcelle_I','ReadVariableNames',false,'TextType','char','DatetimeType','text');
Data_I=table2cell(Data_I);
Data_J=readtable(fic_hauteur,'Sheet','Parcelle_J','ReadVariableNames',false,'TextType','char','DatetimeType','text');
Data_J=table2cell(Data_J);

dates_I=datenum(Data_I(1,2:end),'dd/mm/yyyy');
dates_J=datenum(Data_J(1,2:end),'dd/mm/yyyy');
dmin_I=min(dates_I);
dmin_J=min(dates_J);

h_I=cellfun(@str2num,Data_I(end-1,2:end));
h_J=cellfun(@str2num,Data_J(end-1,2:end));

h_I=h_I/100;
h_J=h_J/100;

% Ecriture des obs en .mat 
% Après 2013, Nord=J=Mil   , Sud=I=Jach
% if annee>2013
%     save(['LAI_N_' annee_str '.mat'],'dates_LAI','LAI_J');
%     save(['LAI_S_' annee_str '.mat'],'dates_LAI','LAI_I');
%     save(['hauteur_N_' annee_str '.mat'],'dates_J','h_J');
%     save(['hauteur_S_' annee_str '.mat'],'dates_I','h_I');
% end
%% Ajouts de points fixes à 0 : au jour 1, 15j avant le début des mesures et au jour 365 pour contraindre l'interpolation hors période de mesure

dates_LAI=[datedeb;dmin_LAI-15;dates_LAI;dateend];
dates_I=[datedeb;dmin_I-15;dates_I;dateend];
dates_J=[datedeb;dmin_J-15;dates_J;dateend];

LAI_I=[0;0;LAI_I;0];
LAI_J=[0;0;LAI_J;0];
h_I=[0.05 0.05 h_I 0.05];
h_J=[0.05 0.05 h_J 0.05];

% Ecriture  en .mat 
% Après 2013, Nord=J=Mil   , Sud=I=Jach
if annee>2013
    save(['LAI_N_' annee_str '_lvl2.mat'],'dates_LAI','LAI_J');
    save(['LAI_S_' annee_str '_lvl2.mat'],'dates_LAI','LAI_I');
    save(['hauteur_N_' annee_str '_lvl2.mat'],'dates_J','h_J');
    save(['hauteur_S_' annee_str '_lvl2.mat'],'dates_I','h_I');
end

% %% Interpolation journalière
% 
% LAI_I_j=naninterp1(dates_LAI,LAI_I,dates_j,'PCHIP');
% LAI_J_j=naninterp1(dates_LAI,LAI_J,dates_j,'PCHIP');
% h_I_j=naninterp1(dates_I,h_I,dates_j,'PCHIP');
% h_J_j=naninterp1(dates_J,h_J,dates_j,'PCHIP');
% 
% %% Check Graphique
% 
% % LAI
% figure
% plot(dates_j,LAI_I_j,'b-')
% hold on
% plot(dates_LAI,LAI_I,'b+')
% plot(dates_j,LAI_J_j,'r-')
% plot(dates_LAI,LAI_J,'r+')
% datetick('x')
% ylabel('LAI')
% title(annee_str)
% legend('I (interp.)','I (obs.)','J (interp.)','J (obs.)','location','eastoutside')
% 
% % Hauteur
% figure
% plot(dates_j,h_I_j,'b-')
% hold on
% plot(dates_I,h_I,'b+')
% plot(dates_j,h_J_j,'r-')
% plot(dates_J,h_J,'r+')
% datetick('x')
% ylabel('Hauteur (m)')
% title(annee_str)
% legend('I (interp.)','I (obs.)','J (interp.)','J (obs.)','location','eastoutside')
% 
% %% Ecriture des fichiers
% 
 % Après 2013, Nord=J=Mil   , Sud=I=Jach
 
 fic1=[path_res '\lai_N_' annee_str '.dat'];
 fic2=[path_res '\lai_S_' annee_str '.dat'];
 
 % Calcul des DoY
 DoY=dates_j-datedeb+1;
 
 if annee>2013
     fid1=fopen(fic1,'w+');
     fid2=fopen(fic2,'w+');
     
     for i=1:length(DoY)
         fprintf(fid1,'%i %7.2f %7.2f %7.2f\n',DoY(i),LAI_J_j(i),0,h_J_j(i));
         fprintf(fid2,'%i %7.2f %7.2f %7.2f\n',DoY(i),LAI_I_j(i),0,h_I_j(i));
     end
 end
 fclose(fid1);
 fclose(fid2);
