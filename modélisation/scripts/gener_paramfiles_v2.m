%% Routine de création de l'ensemble des fichiers de paramètres pour MCIP
%

function gener_paramfiles_v2(root0,root1,root_iter,root_ref,paramref,nameparc,annee_ini,inbsimu,nparam,matparam)
%% Preparation et stockage des fichiers de parametre de SiSPAT

%Localisation du fichier type de parametres
temp = strcat(paramref(1:6),nameparc(1),paramref(8),num2str(annee_ini),paramref(13:end));
file_ref = strcat(root0,root1,root_ref,temp);
params_ref=Read_Sispat_param(file_ref);

%Preparation des fichiers de parametres
for isim = 1:inbsimu
    %Localisation du fichier de parametres en cours
    temp = strcat(paramref(1:6),nameparc(1),paramref(8),num2str(annee_ini),'_',sprintf('%05d',isim),'.dat');
    file_new =  strcat(root0,root1,'Source/File_param_sispat/',temp);
    
    params=params_ref;
    
    numsim=sprintf('%05d',isim);
    
    % Fichiers de sorties
    params{13}=sprintf('DIR_RES/sol/sol_%s_%s_%s.txt',nameparc(1),num2str(annee_ini),numsim);
    params{15}=sprintf('DIR_RES/atm/atm_%s_%s_%s.out',nameparc(1),num2str(annee_ini),numsim);
    params{17}=sprintf('DIR_RES/fin/fin_%s_%s_%s.out',nameparc(1),num2str(annee_ini),numsim);
    params{19}=sprintf('DIR_RES/flx/flx_%s_%s_%s.out',nameparc(1),num2str(annee_ini),numsim);
    params{end}=sprintf('DIR_RES/term/term_%s_%s_%s.txt',nameparc(1),num2str(annee_ini),numsim);
    
    % Porosité
    por=strsplit(params{57},' ');
    por=cellfun(@str2double,por);
    
    por(1:3)=matparam(isim,1:3);
    params{57}=num2str(por);
    
    % Wsat
    params{63}=num2str(0.9.*por);
    
    % hg (Van-Genuchten)
    hg=strsplit(params{59},' ');
    hg=cellfun(@str2double,hg);
    
    hg(1:3)=matparam(isim,4:6);
    params{59}=num2str(hg);
    
    % q (Van-Genuchten)
    q=strsplit(params{61},' ');
    q=cellfun(@str2double,q);
    
    q(1:3)=matparam(isim,7:9);
    params{61}=num2str(q);
    
    % Ksat (Brooks-Corey)
    Ks=strsplit(params{77},' ');
    Ks=cellfun(@str2double,Ks);
    
    Ks(1:3)=exp(matparam(isim,10:12).*log(10));
    params{77}=num2str(Ks);   
    
    % beta
    beta=strsplit(params{83},' ');
    beta=cellfun(@str2double,beta);
    
    beta(1:3)=matparam(isim,13:15);
    params{83}=num2str(beta); 
    
    % ctherm
    ctherm=strsplit(params{107},' ');
    ctherm=cellfun(@str2double,ctherm);
    
    ctherm(1:3)=matparam(isim,16:18);
    params{107}=num2str(ctherm); 
    
    % FDRmax
    for i=149:156
        PR=strsplit(params{i},' ');
        PR=cellfun(@str2double,PR);
        PR(end)=matparam(isim,19)*PR(end)/25000; % on conserve la proportionalité de la dynamique racinaire
        params{i}=num2str(PR);
    end
    
    % Rtot
    Rtot=exp(matparam(isim,20).*log(10));
    params{139}=num2str(Rtot);
    
    fid=fopen(file_new,'w');
    for i=1:length(params)
        fprintf(fid,'%s\n',params{i});
    end
    fprintf(fid,'\n');
    fclose(fid);
end
return
end