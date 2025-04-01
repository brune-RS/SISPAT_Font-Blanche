%% Programme de mise en forme des donnees : Format SiSPAT (pas de tps 30mn) 

%   - stations Mil et Jachere - 2005-2011
%  Auteur : Cecile VELLUET - 07 Fevrier 2010 - Modifié en Janvier 2013 (AMMA Niger)
% Modification: Brune RAYNAUD--SCHELL février 2022
% pour données de FONT-BLANCHE

clear
clc
close all
disp('Mise en forme des fichiers de données Font-Blanche - Format SiSPAT')


%% Paramètres
% -----------
%Fichiers de forçage meteo Font-Blanche

path_lvl2='C:\SISPAT\tableaux matlab\Font-Blanche\';
file_meteo=[path_lvl2 'Meteo2014-2021(timetable).mat' ];
file_radiation=[path_lvl2 'FBradiation20142021.mat'];


%% Lecture des donnees meteo (level 2, format matlab)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
meteo_FBtt=importdata(file_meteo); %timetable
time_index=meteo_FBtt.TIMESTAMP;
meteo_FB=table2array(removevars((timetable2table(meteo_FBtt)),{'TIMESTAMP'}));

I=isnan(meteo_FB); 
J=find(~I); 
meteo_FB(I)=interp1(J,meteo_FB(J),find(I));

% (index)time jj/MM/aaaa hh:mm   (1)Rg  (2)PPFD  (3)Patm  (4)Ta  (5)RH
% (6)VPD  (7)Precip  (8)Wind Speed  (9)WS max (10)Wdir (11)ETP
rad_FBtt=importdata(file_radiation);
rad_FB=table2array(removevars(rad_FBtt,{'TIMESTAMP'}));

I=isnan(rad_FB); 
J=find(~I); 
rad_FB(I)=interp1(J,rad_FB(J),find(I));


% % ()DoY  (2)Heure  (3)Min  (4)Rg  (5)Rg(diffus)  (6)Ra  (7)Ta  (8)RH
% (9)Ua  (10)Precip  (11)Patm  (12)ZaU  (13)ZaT (14)CO2
%% Preparation des fichiers SiSPAT
% --------------------------------------------------------------------------------------------

%Jour julien                           					(Col 1, [-])
% ------------------------------------------------------------------------
time_index.Format='dd-MM-yyyy';
t=datenum(time_index);
SiSPAT_climate_FB(:,1)  = fix(t-t(1)+1);             %Unitless

%Heure	                           						(Col 2, [-])
% ------------------------------------------------------------------------
time_index=meteo_FBtt.TIMESTAMP;
time_index.Format='hh';
SiSPAT_climate_FB(:,2)  = str2num(  datestr(time_index,'hh'));               %Unitless

%Minute                          						(Col 3, [-])
% ------------------------------------------------------------------------
time_index=meteo_FBtt.TIMESTAMP;
SiSPAT_climate_FB(:,3)  =minute(time_index);               %Unitless

%% RAYONNEMENTS

%Rayonnement solaire incident (Rg)          			(Col 4, [W.m-2])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,4) = meteo_FB(:,1);               %W.m-2
SiSPAT_climate_FB((SiSPAT_climate_FB(:,2))<0,4)=0.0;

%Rayonnement solaire incident diffus (Rg-diffus)        (Col 5, [W.m-2])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,5)  = 0.0;                                 %W.m-2
%SiSPAT_climate_S(:,5) = 0.0;        	                      %W.m-2

%Rayonnement atmospherique incident (Ra)    		    (Col 6, [W.m-2])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,6)  = rad_FB(:,3);               %W.m-2
%SiSPAT_climate_S(:,6) = mto_sud_2(:,6);      	  %W.m-2

%Rayonnement atmospherique incident entre 8 et 14 (Ra)  (Col 7, [W.m-2])
% ------------------------------------------------------------------------
%IDSO (1981) et Olioso (AFM, 1995)
f8_14_N  = -0.6732 + 0.6240e-2*SiSPAT_climate_FB(:,8)  - 0.9140e-5*(SiSPAT_climate_FB(:,8).^2); 

e8_14_N  = 0.15 + 5.03e-6 .* (ew/100)  .* exp(2450./SiSPAT_climate_FB(:,8));

SiSPAT_climate_FB(:,7)   = 5.67e-8*f8_14_N.*e8_14_N.*(SiSPAT_climate_FB(:,8).^4);   %W.m-2

%% METEO

%Temperature de l'air (Ta)                  		    (Col 8, [K])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,8) = meteo_FB(:,4)+273.15;         %Â°C --> K
%SiSPAT_climate_S(:,8)= mto_sud_2(:,7)+273.15; 	  %Â°C --> K


%Humidité spécifique de l'air (HS)          		    (Col 9, [kg.kg-1])
% ------------------------------------------------------------------------
ewsat  = 610.7*(1+sqrt(2)*sin(pi*((SiSPAT_climate_FB(:,8)-273.15)/540))).^8.827;

ew    = (meteo_FB(:,5).*ewsat)./100; % colonne 5 : RH ?? à vérifier

Pa	   = meteo_FB(:,3)*100;
r=(0.622*ew);
s=(Pa-0.378*ew); % Patm			   %hPa  --> Pa vérifier 
SiSPAT_climate_FB(:,9)=r./s;


%Vitesse du vent (Ua)                       		    (Col 10, [m.s-1])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,10)  = meteo_FB(:,8);      	      %m.s-1 Ã  2.88m
%SiSPAT_climate_S(:,10) = mto_sud_2(:,8);  	      %m.s-1 Ã  2.88m


%Precipitations (PrÃ©cip)                         		(Col 11, [mm/30mn])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,11) = meteo_FB(:,7);              %mm/30mn
%SiSPAT_climate_S(:,11)= mto_sud_2(:,7);          %mm/30mn


%Pression atmosphérique (Patm)                        	(Col 12, [hPa])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,12) = meteo_FB(:,3);              %hPa
%SiSPAT_climate_S(:,12)= mto_sud_2(:,3);          %hPa


%% MESURES 

%Hauteur de mesure du vent (ZaU)                        (Col 13, [m])
% ------------------------------------------------------------------------
SiSPAT_climate_FB(:,13) = 11.5;                                %m
%SiSPAT_climate_S(:,13)= 2.88;                                %m

%Hauteur de mesure de la température (ZaT)              (Col 14, [m])
%------------------------------------------------------------------------
SiSPAT_climate_FB(:,14) = 10;                                %m
%SiSPAT_climate_S(:,14)= 2.88;                                %m

%% Sauvegarde des données
% ---------------------------------------------------------------------------------------------------------------------------------------------

%Format matlab
%save(SiSPAT_climate_FB);
% SiSPAT_climate_FB :  (1)DoY  (2)Heure  (3)Min  (4)Rg  (5)Rg(diffus)  (6)Ra  (7)Ra(8-14)  (8)Ta  (9)HS  (10)Ua  (11)Prï¿½cip  (12)Patm  (13)ZaU  (14)ZaT
% SiSPAT_climate_S : (1)DoY  (2)Heure  (3)Min  (4)Rg  (5)Rg(diffus)  (6)Ra  (7)Ra(8-14)  (8)Ta  (9)HS  (10)Ua  (11)Prï¿½cip  (12)Patm  (13)ZaU  (14)ZaT

%fid1 = fopen(file_SISPAT_nord,'w+');
%for i=1:length(SiSPAT_climate_FB(:,1))
    %fprintf(fid1,'%i %i %i %9.2f %9.2f %9.2f %9.2f %9.2f %10.5f %9.3f %11.2f %9.2f %9.2f %9.2f\n',SiSPAT_climate_FB(i,:));
%end
%fclose(fid1);

%% Découpage pour période de calibration
% indices des fichiers 
annee_s=2014;
mois_s=8;
jour_s=1;
heure_s=0;
minute_s=0;
date_s=datenum(annee_s,mois_s,jour_s,heure_s,minute_s,0);

annee_e=2021;
mois_e=12;
jour_e=31;
heure_e=23;
minute_e=30;
date_e=datenum(annee_e,mois_e,jour_e,heure_e,minute_e,0);

i_start=find((abs(datenum(meteo_FBtt.TIMESTAMP)-date_s))==min(abs(datenum(meteo_FBtt.TIMESTAMP)-date_s)));
j_start=find((abs(datenum(rad_FBtt.TIMESTAMP)-date_s))==min(abs(datenum(rad_FBtt.TIMESTAMP)-date_s)));
i_end=find((abs(datenum(meteo_FBtt.TIMESTAMP)-date_e))==min(abs(datenum(meteo_FBtt.TIMESTAMP)-date_e)));
j_end=find((abs(datenum(rad_FBtt.TIMESTAMP)-date_e))==min(abs(datenum(rad_FBtt.TIMESTAMP)-date_e)));


SiSPAT_climate_FB_14_21=SiSPAT_climate_FB(i_start:i_end,:);
SiSPAT_climate_FB_14_21(:,1)=SiSPAT_climate_FB_14_21(:,1)-SiSPAT_climate_FB_14_21(1,1)+1;

fid = fopen('SiSPAT_clim_FB_14_21.dat', 'w');

for i=1:length(SiSPAT_climate_FB_14_21(:,1))
    fprintf(fid,'%i %i %i %9.2f %9.2f %9.2f %9.2f %9.2f %10.5f %9.3f %11.2f %9.2f %9.2f %9.2f\n',SiSPAT_climate_FB_14_21(i,:));
end
fclose(fid);
