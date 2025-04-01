
%% Routine de création de l'ensemble des fichiers de paramètres pour MCI
clear;
close all;
clc;
addpath('/media/brune/Ultra Touch/SISPAT/Fonctions matlab/')
addpath('SISPAT/Fonctions matlab/')

% root0='/home/brune/SISPAT/execution/';
% root1= 'param/ensemble/ensemble_6/';
root0='/media/brune/Ultra Touch/SISPAT/';
root1='param_ens6/';
% root_iter 
root_ref='Aout_modele.dat';
% paramref %?
% nameparc
% annee_ini
% inbsimu  %nombre de simulations 
% nparam
% matparam % tableau avec une ligne pour chaque simu et les colonnes pour les différentes variables 

%% %% matrice paramètres à la main  _________________________________

% 1) prof fractures min  prof fractures max  2) fraction pluie reprise fraction ruissellement   3) porosité   4)Wsat  5)hg 
% % 6)q   7)Ksat   8) beta  9)ctherm  10)FDR max 11)%pmr 12)res stomatique 
% 13) Rtot

% 1) prof max fractures   2) fraction pluie reprise    3)FDR max   4)%pmr
% racines   5) resistance stomatique
% les variables que nous modifions 
p1=[1.5 10; 1.5 15 ;1.5 20];% prof max fissures
%p1num=[1 2 3];
p2=[0.10 0.95;0.20 0.95;0.30 0.95];% pourcentage reprise pluie
%p2num=[1 2 3];
p10=[300 500 700];
p11=[0.10 0.25 0.40];
p12=[1.2 1.5 2];
matparam=allcomb(p1(:,2),p2(:,1),p10,p11,p12); % toutes les combinaisons de paramètres 3 puissance 5 possibilités 
p=[1 2 10 11 12];
% creer des colonnes vides pour les paramètres à ne pas changer 
% pour garder des indices constants
%variables non changées seaddpath('/home/brune/SISPAT/Fonctions matlab/')ulement desactivés à la main dans func
for j=1:12
    if ismember(j,p)
       jd=1;
    else
        matparamav=matparam(:,1:j-1) ;
        matparamap=matparam(:,j:end);
        matze=zeros(size(matparam(:,1))) ;
        matparam=[matparamav matze matparamap];% 12 colonnes 
    end
end
inbsimu=length(matparam);
% ??        

%%   matrices paramètres hypercube 
variables=[" prof. fractures " "fraction reprise "   "porosité"   "Wsat"  "hg" "q"   "Ksat"   "beta"  "ctherm"  "FDR max" "%pmr" 'res.st ' 'res.tot.'];

% les variables que nous modifions 
%p2num=[1 2 3];
p1=[0 3]; % prof fracture
p2=[0.5 1];  % reprise ruiss
p10=[100 2000]; % FDR 
p13=[9.5 12]; % Rtot
p12=[250 500];% R stomates min
%matparam=allcomb(p1(:,2),p2(:,1),p10,p11,p12); % toutes les combinaisons de paramètres 3 puissance 5 possibilités 
p=[1 2 10 12 13];
% fonction hypercube 
[X_param,X_normalized]=lhsdesign_modified(2000,[p1(1) p2(1) p10(1) p12(1) p13(1)],[p1(2) p2(2) p10(2) p12(2) p13(2)]);
matparam=X_param;     
for j=1:13
    if ismember(j,p)
       jd=1;
    else
        matparamav=matparam(:,1:j-1) ;
        matparamap=matparam(:,j:end);
        matze=zeros(size(matparam(:,1))) ;
        matparam=[matparamav matze matparamap];% 12 colonnes 
    end
end
inbsimu=length(matparam);
       

% espace de paramètres illustré en 3D 
% figure()
%  scatter3(matparam(:,10),matparam(:,12) ,matparam(:,13))
%  xlabel(variables(10))
%  ylabel(variables(12))
%  zlabel(variables(13))
%  grid on
%  ax = gca;
%  ax.FontSize = 10; 
%  ax.FontWeight = 'bold'; 
               
               

% Write the results to a CSV file or .mat    save
%matparam=table2array(mat_param.matparam);
save(strcat(root0,root1,'/matparam6.mat'),'matparam')
tab_param=array2table(matparam);
tab_param.Properties.VariableNames(1:13)= {' prof.fractures ' ,'frac.pluie ',  'porosité ' ,  'Wsat '  ,'hg ', 'q ' ,  'Ksat '  , 'beta ',  'ctherm ' , 'FDRmax ' ,'%pmr ', 'res.st ', 'res.tot.'};
writetable(tab_param,strcat(root0,root1,'/matparam6.csv')); 
% fid = fopen('matparam2.csv', 'w');
% fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', header{:});
% fclose(fid);
% dlmwrite(strcat(root0,root1,'/matparam2.csv'), matparam);
%%
%function gener_paramfiles_v2(root0,root1,root_iter,root_ref,paramref,nameparc,annee_ini,inbsimu,nparam,matparam)
root00='/mnt/c/Users/Brune/SISPAT/';
root1='param_ens6/';
gener_paramfiles_v3(root0,root00,root1,root_ref,inbsimu,matparam)% tous les p necessaires ou non

%% 
function gener_paramfiles_v3(root0,root1,root_ref,inbsimu,matparam)%
%% Preparation et stockage des fichiers de parametre de SiSPAT

%Localisation du fichier type de parametres
% temp = strcat(paramref(1:6),nameparc(1),paramref(8),paramref(13:end));
file_ref = strcat(root0,root1,root_ref);
% file_ref
params_ref=Read_Sispat_param(file_ref);

%Preparation des fichiers de parametres
for isim = 1:2
    %Localisation du fichier de parametres en cours
    temp = strcat('param','_',sprintf('%04d',isim),'.dat');
    file_new =  strcat(root0,root1,temp);
    
    params=params_ref;
    
    numsim=sprintf('%05d',isim);

    % Fichiers de sorties
    path_out='/media/brune/Ultra Touch/SISPAT/out';
    params{13}=strcat(path_out,sprintf('/sol/sol_%s.txt',numsim));
    params{15}=strcat(path_out,sprintf('/atm/atm_%s.out',numsim));
    params{17}=strcat(path_out,sprintf('/fin/fin_%s.out',numsim));
    params{19}=strcat(path_out,sprintf('/flx/flx_%s.out',numsim));
    params{end}=strcat(path_out,sprintf('/term/term_%s.txt',numsim));
    
   % Version karst, prise en compte des fissures
   %profondeur des fissures 
    prof=strsplit(params{60},' ');
    prof=cellfun(@str2double,prof);
    
    prof(2)=matparam(isim,1);
    params{60}=num2str(prof(1));

    %pourcentage de reprise pluie et ruissellement 
    reprise=strsplit(params{62},' ');
    reprise=cellfun(@str2double,reprise);
    
    reprise(1:3)=matparam(isim,2);
    params{62}=num2str(reprise(1));
   
    % % Porosité
    % por=strsplit(params{70},' ');
    % por=cellfun(@str2double,por);
    % 
    % por(1:3)=matparam(isim,3);
    % params{70}=num2str(por);
    % 
    % % Wsat
    % params{76}=num2str(0.9.*por);
    % 
    % % hg (Van-Genuchten)
    % hg=strsplit(params{72},' ');
    % hg=cellfun(@str2double,hg);
    % 
    % hg(1:3)=matparam(isim,5};
    % params{72}=num2str(hg);
    % 
    % % q (Van-Genuchten)
    % q=strsplit(params{74},' ');
    % q=cellfun(@str2double,q);
    % 
    % q(1:3)=matparam(isim,6);
    % params{74}=num2str(q);
    % 
    % % Ksat (Brooks-Corey)
    % Ks=strsplit(params{90},' ');
    % Ks=cellfun(@str2double,Ks);
    % 
    % Ks(1:3)=exp(matparam(isim,7).*log(10));
    % params{90}=num2str(Ks);   
    % 
    % % beta
    % beta=strsplit(params{96},' ');
    % beta=cellfun(@str2double,beta);
    % 
    % beta(1:3)=matparam(isim,8);
    % params{96}=num2str(beta); 
    % 
    % % ctherm
    % ctherm=strsplit(params{118},' ');
    % ctherm=cellfun(@str2double,ctherm);
    % 
    % ctherm(1:3)=matparam(isim,9};
    % params{118}=num2str(ctherm); 
    
    % FDRmax
    for i=169:251
        PR=strsplit(params{i},' ');
        PR=cellfun(@str2double,PR);

        PR(end)=matparam(1,10)*PR(end)/25000; % on conserve la proportionalité de la dynamique racinaire
        params{i}=num2str(PR);
    end
    
    % Rtot
    Rtot=exp(matparam(2,11).*log(10));
    params{152}=num2str(Rtot);

    % R stomates min
    Rst_min=exp(matparam(isim,12));
    params{148}=num2str(Rst_min);

    
    fid=fopen(file_new,'w');
    for i=1:length(params)
        fprintf(fid,'%s\n',params{i});
    end
    fprintf(fid,'\n');
    fclose(fid);
end
return
end
