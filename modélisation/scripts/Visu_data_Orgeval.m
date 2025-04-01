%% Init
clear;
clc;
close all;
disp('Script de visualisation des données Orgavel - Auteur Jerome 12/12/22');


%% Initialisation
SMCp_FLX_1mn_P1 = []; SMCp_FLX_1mn_P2 = []; SMCp_FLX_1mn_P3 = []; SMCp_FLX_1mn_P4 = []; idate_SMC_FLX_1mn = [];
SMC_Boissyaxel_P1 = []; SMC_Boissyaxel_P2 = []; SMC_Boissyaxel_P3 = []; SMC_Boissyaxel_P4 = []; idate_SMC_Boissyaxel = [];


%% Relation de convertions Periode en humidité

%Relations classiques (electrical conductivity<0.5 dS/m, bulk density<1.55 g/cm3 and %Clay<30%)								
RC_l = [-0.467 0.0283];	          %linéaire		
RC_q = [-0.0663 -0.0063 0.0007];  %quadratique

%Relations modifiées 1 : Sandy Clay loam ; bulk density =1.6g/cm3 ; elec cond = 0.4dS/m								
%RM1_l = [-0.62 0.0329];	          %linéaire		
%RM1_q = [0.095 -0.0211 0.001];    %quadratique
RM1_l = [-0.467 0.0283];            %linéaire		
RM1_q = [-0.0663 -0.0063 0.0007];   %quadratique
Bias1_q = [-0.038 -0.135 -0.21 -0.235];
Bias1_l = [0.0 -0.11 -0.19 -0.23];
								
%Relations modifiées 2 : Sandy Clay loam ; bulk density =1.6g/cm3 ; elec cond = 0.75dS/m								
RM2_l = [-0.447 0.0254];          %linéaire		
RM2_q = [-0.018 -0.007 0.0006];   %quadratique
Bias2_q = [0.04 -0.04 -0.11 -0.12];
Bias2_l = [0.04 -0.02 -0.09 -0.095];

%Relations modifiées 3 (papier fourni par Axel)
%https://doi.org/10.1016/j.geoderma.2009.01.002
RM3_q = [-0.207 0.097 0.288];


%%  Configuration
itrait_SMC_FLX_1mn = 1;      % 0 = Lecture et sauvegarde des données SMC - station FLX - 1mn
itrait_SMC_FLX_30mn= 0;      % 0 = Lecture et sauvegarde des données SMC - station FLX - 30mn - BDOH
itrait_SMC_ANX_30mn= 1;      % 0 = Lecture et sauvegarde des données SMC - station ANX - 1mn  - BDOH
itrait_SMC_BOISSY  = 1;      % 0 = Lecture et sauvegarde des données SMC - station BOISSY - 12h  - BDOH
itrait_SMC_SUIZY   = 1;      % 0 = Lecture et sauvegarde des données SMC - station SUISSY - 12h  - BDOH
itrait_T_FLX_30mn  = 1;      % 0 = Lecture et sauvegarde des données Tsol - station FLX - 30mn - BDOH
itrait_T_ANX_30mn  = 1;      % 0 = Lecture et sauvegarde des données Tsol - station ANX - 1mn  - BDOH
itrait_G_FLX_30mn  = 1;      % 0 = Lecture et sauvegarde des données flux G - station FLX - 30mn - BDOH
itrait_G_ANX_30mn  = 1;      % 0 = Lecture et sauvegarde des données flux G - station ANX - 1mn  - BDOH
itrait_LE_FLX_30mn = 1;      % 0 = Lecture et sauvegarde des données flux H et LE - station FLX - 30mn - Mise en forme Axel
itrait_LE2_FLX_30mn= 1;      % 0 = Lecture et sauvegarde des données flux H et LE - station FLX - 30mn - Eddypro Jerome
itrait_RN_FLX_30mn = 1;      % 0 = Lecture et sauvegarde des données flux G - station FLX - 30mn - Mise en forme Axel
itrait_Ray_FLX_30mn= 1;      % 0 = Lecture et sauvegarde des données flux G - station FLX - 30mn - BDOH
itrait_Ray_ANX_30mn= 1;      % 0 = Lecture et sauvegarde des données flux G - station ANX - 30mn - BDOH 
itrait_Ray_BOI_1h  = 1;      % 0 = Lecture et sauvegarde des données flux G - station BOISSY - 30mn - BDOH 
itrait_Pluie_1h    = 1;      % 0 = Lecture et sauvegarde des données flux G - station ANX - 30mn - BDOH 

%% SMC FLUX à 1 mn dérivé de la centrale - Lieu data HSM 
if (itrait_SMC_FLX_1mn == 0)
    disp('  - Traitement de SMC FLUX à 1 mn dérivé de la centrale');
    
    %Entetes fichiers
    rootfile = 'DataAxel/';
    filein = {'Fichier_DATA_SMSC_2016.csv'  ; 'Fichier_DATA_SMSC_2017.csv' ; 'Fichier_DATA_SMSC_2018.csv'};
    DataTable_Separateur=';';
    DataTable_Lentete=1;
    DataTable_Nbcol=13;
    DataTable_champ={'Temps','Sol_H_1_Avg','Sol_H_2_Avg','Sol_H_3_Avg','Sol_H_4_Avg','Sol_H_1_Quadratic_2016', ...
        'Sol_H_2_Quadratic_2016','Sol_H_3_Quadratic_2016','Sol_H_4_Quadratic_2016','Sol_H_1_Linear_2016',...
        'Sol_H_2_Linear_2016','Sol_H_3_Linear_2016','Sol_H_4_Linear_2016'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        c=char(Databrut.textdata(:,1));
        ijour = str2double(string(c(2:end,1:2)));
        imois = str2double(string(c(2:end,4:5)));
        ian   = str2double(string(c(2:end,7:10)));
        iheur = str2double(string(c(2:end,12:13)));
        imin  = str2double(string(c(2:end,15:16)));
        idate_SMC_FLX_1mn(size(idate_SMC_FLX_1mn,1)+1:size(idate_SMC_FLX_1mn,1)+size(c,1)-1,1)=datenum(ian,imois,ijour,iheur,imin,0);
        
        %Data
        for ivar=2:DataTable_Nbcol-1
            DataMat.(DataTable_champ{ivar})=Databrut.data(:,ivar-1);
        end
        SMCp_FLX_1mn_P1(size(SMCp_FLX_1mn_P1,2)+1:size(SMCp_FLX_1mn_P1,2)+size(DataMat.Sol_H_1_Avg,1))= DataMat.Sol_H_1_Avg;
        SMCp_FLX_1mn_P2(size(SMCp_FLX_1mn_P2,2)+1:size(SMCp_FLX_1mn_P2,2)+size(DataMat.Sol_H_2_Avg,1))= DataMat.Sol_H_2_Avg;
        SMCp_FLX_1mn_P3(size(SMCp_FLX_1mn_P3,2)+1:size(SMCp_FLX_1mn_P3,2)+size(DataMat.Sol_H_3_Avg,1))= DataMat.Sol_H_3_Avg;
        SMCp_FLX_1mn_P4(size(SMCp_FLX_1mn_P4,2)+1:size(SMCp_FLX_1mn_P4,2)+size(DataMat.Sol_H_4_Avg,1))= DataMat.Sol_H_4_Avg;
    end
    
    %Sauvegarde
    save('SMCp_FLX_1mn.mat','SMCp_FLX_1mn_P1','SMCp_FLX_1mn_P2','SMCp_FLX_1mn_P3','SMCp_FLX_1mn_P4','idate_SMC_FLX_1mn');
else
    %Chargement
    disp('  - Chargement direct de SMC FLX à 1 mn dérivée de la centrale');
    load('SMCp_FLX_1mn.mat');
end

%% SMC TAILLIS/FLX à 30 mn - Source BDOH 
if (itrait_SMC_FLX_30mn == 0)
    disp('  - Traitement de SMC FLUX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_FLX/';
    filein = {'TAILLIS_HUMS.txt';'TAILLIS_HUMS-2.txt';'TAILLIS_HUMS-3.txt';'TAILLIS_HUMS-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_SMC_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SMC_FLX_30mn_P1 = Databrut.data(:,1);
            case(2)
                SMC_FLX_30mn_P2 = Databrut.data(:,1);
            case(3)
                SMC_FLX_30mn_P3 = Databrut.data(:,1);
            case(4)
                SMC_FLX_30mn_P4 = Databrut.data(:,1);
        end
    end

    %Detection des NaN
    dump=SMC_FLX_30mn_P1<0;
    SMC_FLX_30mn_P1(dump)=NaN;
    dump=SMC_FLX_30mn_P2<0;
    SMC_FLX_30mn_P2(dump)=NaN;
    dump=SMC_FLX_30mn_P3<0;
    SMC_FLX_30mn_P3(dump)=NaN;
    dump=SMC_FLX_30mn_P4<0;
    SMC_FLX_30mn_P4(dump)=NaN;

    %Racines des SMC pour revenir aux périodes
    for i=1:length(idate_SMC_FLX_30mn)
        p = [RC_q(3) RC_q(2) RC_q(1)-SMC_FLX_30mn_P1(i)/100];
        if isnan(p(3))
            SMCp_FLX_30mn_P1(i) = NaN;
        else
            r = roots(p);
            SMCp_FLX_30mn_P1(i) = max(r);
        end
        p = [RC_q(3) RC_q(2) RC_q(1)-SMC_FLX_30mn_P2(i)/100];
        if isnan(p(3))
            SMCp_FLX_30mn_P2(i) = NaN;
        else
            r = roots(p);
            SMCp_FLX_30mn_P2(i) = max(r);
        end
        p = [RC_q(3) RC_q(2) RC_q(1)-SMC_FLX_30mn_P3(i)/100];
        if isnan(p(3))
            SMCp_FLX_30mn_P3(i) = NaN;
        else
            r = roots(p);
            SMCp_FLX_30mn_P3(i) = max(r);
        end
        p = [RC_q(3) RC_q(2) RC_q(1)-SMC_FLX_30mn_P4(i)/100];
        if isnan(p(3))
            SMCp_FLX_30mn_P4(i) = NaN;
        else
            r = roots(p);
            SMCp_FLX_30mn_P4(i) = max(r);
        end
    end
    
    %Réalignement des sondes (induit par la remise en place de la sonde le 13/11/2015)
    %2016
    delta_p1 = SMCp_FLX_30mn_P1(15394)-SMCp_FLX_30mn_P1(15190);
    SMCp_FLX_30mn_P1(1:10785) = SMCp_FLX_30mn_P1(1:10785)+delta_p1;
    %2017
    delta_p2 = 1.42;
    SMCp_FLX_30mn_P2(31904:end) = SMCp_FLX_30mn_P2(31904:end)+delta_p2;
    delta_p1 = 1.42;
    SMCp_FLX_30mn_P1(31904:end) = SMCp_FLX_30mn_P1(31904:end)+delta_p2;

    %Sauvegarde
    save('SMCp_FLX_30mn.mat','SMCp_FLX_30mn_P1','SMCp_FLX_30mn_P2','SMCp_FLX_30mn_P3','SMCp_FLX_30mn_P4','idate_SMC_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de SMC FLX à 30 mn dérivée de la BDOH');
    load('SMCp_FLX_30mn.mat');
end


%% SMC NOUE/ANX à 30 mn - Source BDOH
if (itrait_SMC_ANX_30mn == 0)
    disp('  - Traitement de SMC ANX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_ANX/';
    filein = {'NOUE_HUMS.txt';'NOUE_HUMS-2.txt';'NOUE_HUMS-3.txt';'NOUE_HUMS-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_SMC_ANX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SMC_ANX_30mn_P1 = Databrut.data(:,1);
            case(2)
                SMC_ANX_30mn_P2 = Databrut.data(:,1);
            case(3)
                SMC_ANX_30mn_P3 = Databrut.data(:,1);
            case(4)
                SMC_ANX_30mn_P4 = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('SMC_ANX_30mn.mat','SMC_ANX_30mn_P1','SMC_ANX_30mn_P2','SMC_ANX_30mn_P3','SMC_ANX_30mn_P4','idate_SMC_ANX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de SMC ANX à 30 mn dérivée de la BDOH');
    load('SMC_ANX_30mn.mat');
end

%% SMC Boissy aux profondeurs analogues que FLx - Source Axel (depuis BDOH) 
if (itrait_SMC_BOISSY == 0)
    disp('  - Traitement de SMC BOISSY à 12h dérivée des données Axel');

    %Entete du fichier
    rootfile = 'DataAxel/';
    filein = 'DATA_SOL_HMS_BOISSY.csv';
    DataTable_Separateur=';';
    DataTable_Lentete=1;
    DataTable_Nbcol=5;
    DataTable_champ={'DateHeure','H_sol1','H_sol2','H_sol3','H_sol4'};
    
    %Lecture du fichier
    Databrut=importdata(strcat(rootfile,filein),DataTable_Separateur,DataTable_Lentete);
    
    %Date
    c=char(Databrut.textdata(:,1));
    ijour = str2double(string(c(2:end,1:2)));
    imois = str2double(string(c(2:end,4:5)));
    ian   = str2double(string(c(2:end,7:10)));
    iheur = str2double(string(c(2:end,12:13)));
    imin  = str2double(string(c(2:end,15:16)));
    idate_SMC_Boissyaxel = datenum(ian,imois,ijour,iheur,imin,0);
    
    %Stockage data
    for ivar=2:DataTable_Nbcol
        DataMat.(DataTable_champ{ivar})=Databrut.data(:,ivar-1);
    end
    SMC_Boissyaxel_P1= DataMat.H_sol1;
    SMC_Boissyaxel_P2= DataMat.H_sol2;
    SMC_Boissyaxel_P3= DataMat.H_sol3;
    SMC_Boissyaxel_P4= DataMat.H_sol4;
    
    %Sauvegarde
    save('SMC_BOISSYaxel_12h.mat','SMC_Boissyaxel_P1','SMC_Boissyaxel_P2','SMC_Boissyaxel_P3','SMC_Boissyaxel_P4','idate_SMC_Boissyaxel');
    
else
    %Chargement
    disp('  - Chargement direct de SMC ANX à 30 mn dérivée de la BDOH');
    load('SMC_BOISSYaxel_12h.mat');
end

%% SMC BOISSY aux profondeurs analogues de FLx - Source BDOH 
if (itrait_SMC_BOISSY == 0)
    disp('  - Traitement de SMC BOISSY à 12h direct depuis la BDOH');

    %Entete du fichier
    rootfile = 'DataBDOH_Boissy/';
    filein = {'BOISSY-TDR_HUMS.txt';'BOISSY-TDR_HUMS-2.txt';'BOISSY-TDR_HUMS-5.txt';'BOISSY-TDR_HUMS-6.txt';'BOISSY-TDR_HUMS-9.txt';'BOISSY-TDR_HUMS-12.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_SMC_Boissy=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SMC_Boissy_P10 = Databrut.data(:,1);
            case(2)
                SMC_Boissy_P11 = Databrut.data(:,1);
            case(3)
                SMC_Boissy_P20 = Databrut.data(:,1);
            case(4)
                SMC_Boissy_P21 = Databrut.data(:,1);
            case(5)
                SMC_Boissy_P3 = Databrut.data(:,1);
            case(6)
                SMC_Boissy_P4 = Databrut.data(:,1);
        end
    end
    
    %Moyenne
    SMC_Boissy_P1 = (SMC_Boissy_P10+SMC_Boissy_P11)/2;
    SMC_Boissy_P2 = (SMC_Boissy_P20+SMC_Boissy_P21)/2;
    
    %Sauvegarde
    save('SMC_BOISSY_12h.mat','SMC_Boissy_P1','SMC_Boissy_P2','SMC_Boissy_P3','SMC_Boissy_P4','idate_SMC_Boissy');
else
    %Chargement
    disp('  - Chargement direct de SMC BOISSY à 12h dérivée de la BDOH');
    load('SMC_BOISSY_12h.mat');
end

%% SMC SUIZY aux profondeurs analogues de FLx - Source BDOH 
if (itrait_SMC_SUIZY == 0)
    disp('  - Traitement de SMC SUIZY à 12h direct de la BDOH');

    %Entete du fichier
    rootfile = 'DataBDOH_Suizy/';
    filein = {'SUIZY-TDR_HUMS-3.txt';'SUIZY-TDR_HUMS-4.txt';'SUIZY-TDR_HUMS-5.txt';'SUIZY-TDR_HUMS-6.txt';'SUIZY-TDR_HUMS-11.txt';'SUIZY-TDR_HUMS-12.txt';'SUIZY-TDR_HUMS-13.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_SMC_SUIZY_12h=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SMC_SUIZY_12h_P10 = Databrut.data(:,1);
            case(2)
                SMC_SUIZY_12h_P11 = Databrut.data(:,1);
            case(3)
                SMC_SUIZY_12h_P20 = Databrut.data(:,1);
            case(4)
                SMC_SUIZY_12h_P21 = Databrut.data(:,1);
            case(5)
                SMC_SUIZY_12h_P3 = Databrut.data(:,1);
            case(6)
                SMC_SUIZY_12h_P40 = Databrut.data(:,1);
            case(7)
                SMC_SUIZY_12h_P50 = Databrut.data(:,1);
        end
    end
    
    %Moyenne
    SMC_SUIZY_12h_P1 = (SMC_SUIZY_12h_P10+SMC_SUIZY_12h_P11)/2;
    SMC_SUIZY_12h_P2 = (SMC_SUIZY_12h_P20+SMC_SUIZY_12h_P21)/2;
    SMC_SUIZY_12h_P4 = (SMC_SUIZY_12h_P40+SMC_SUIZY_12h_P50)/2;
    
    %Sauvegarde
    save('SMC_SUIZY_12h.mat','SMC_SUIZY_12h_P1','SMC_SUIZY_12h_P2','SMC_SUIZY_12h_P3','SMC_SUIZY_12h_P4','idate_SMC_SUIZY_12h');
else
    %Chargement
    disp('  - Chargement direct de SMC SUIZY à 12h dérivée de la BDOH');
    load('SMC_SUIZY_12h.mat');
end


%% Temperature sol TAILLIS/FLX à 30 mn - Source BDOH
if (itrait_T_FLX_30mn == 0)
    disp('  - Traitement de Tsol FLUX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_FLX/';
    filein = {'TAILLIS_TEMPS.txt';'TAILLIS_TEMPS-2.txt';'TAILLIS_TEMPS-3.txt';'TAILLIS_TEMPS-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_Tsol_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                Tsol_FLX_30mn_P1 = Databrut.data(:,1);
            case(2)
                Tsol_FLX_30mn_P2 = Databrut.data(:,1);
            case(3)
                Tsol_FLX_30mn_P3 = Databrut.data(:,1);
            case(4)
                Tsol_FLX_30mn_P4 = Databrut.data(:,1);
        end
    end
    
    %Filtrage 
    Tsol_FLX_30mn_P1(Tsol_FLX_30mn_P1<=-999.0)=NaN;
    Tsol_FLX_30mn_P2(Tsol_FLX_30mn_P2<=-999.0)=NaN;
    Tsol_FLX_30mn_P3(Tsol_FLX_30mn_P3<=-999.0)=NaN;
    Tsol_FLX_30mn_P4(Tsol_FLX_30mn_P4<=-999.0)=NaN;
    
    %Sauvegarde
    save('Tsol_FLX_30mn.mat','Tsol_FLX_30mn_P1','Tsol_FLX_30mn_P2','Tsol_FLX_30mn_P3','Tsol_FLX_30mn_P4','idate_Tsol_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de Tsol FLX à 30 mn dérivée de la BDOH');
    load('Tsol_FLX_30mn.mat');
end

%% Temperature sol NOUE/ANX à 30 mn - Source BDOH
if (itrait_T_ANX_30mn == 0)
    disp('  - Traitement de Tsol ANX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_ANX/';
    filein = {'NOUE_TEMPS.txt';'NOUE_TEMPS-2.txt';'NOUE_TEMPS-3.txt';'NOUE_TEMPS-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_Tsol_ANX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                Tsol_ANX_30mn_P1 = Databrut.data(:,1);
            case(2)
                Tsol_ANX_30mn_P2 = Databrut.data(:,1);
            case(3)
                Tsol_ANX_30mn_P3 = Databrut.data(:,1);
            case(4)
                Tsol_ANX_30mn_P4 = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('Tsol_ANX_30mn.mat','Tsol_ANX_30mn_P1','Tsol_ANX_30mn_P2','Tsol_ANX_30mn_P3','Tsol_ANX_30mn_P4','idate_Tsol_ANX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de Tsol ANX à 30 mn dérivée de la BDOH');
    load('Tsol_ANX_30mn.mat');
end

%% Flux G TAILLIS/FLX à 30 mn - Source BDOH
if (itrait_G_FLX_30mn == 0)
    disp('  - Traitement de G FLUX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_FLX/';
    filein = {'TAILLIS_FCHAL.txt';'TAILLIS_FCHAL-2.txt';'TAILLIS_FCHAL-3.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_G_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                G_FLX_30mn_P1 = Databrut.data(:,1);
            case(2)
                G_FLX_30mn_P2 = Databrut.data(:,1);
            case(3)
                G_FLX_30mn_P3 = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('G_FLX_30mn.mat','G_FLX_30mn_P1','G_FLX_30mn_P2','G_FLX_30mn_P3','idate_G_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de G FLX à 30 mn dérivée de la BDOH');
    load('G_FLX_30mn.mat');
end

%% flux G NOUE/ANX à 30 mn - Source BDOH
if (itrait_G_ANX_30mn == 0)
    disp('  - Traitement de G ANX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_ANX/';
    filein = {'NOUE_FCHAL.txt';'NOUE_FCHAL-2.txt';'NOUE_FCHAL-3.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_G_ANX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                G_ANX_30mn_P1 = Databrut.data(:,1);
            case(2)
                G_ANX_30mn_P2 = Databrut.data(:,1);
            case(3)
                G_ANX_30mn_P3 = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('G_ANX_30mn.mat','G_ANX_30mn_P1','G_ANX_30mn_P2','G_ANX_30mn_P3','idate_G_ANX_30mn');
    
else
    %Chargement
        disp('  - Chargement direct de G ANX à 30 mn dérivée de la BDOH');
    load('G_ANX_30mn.mat');
end

%% flux H et LE FLX à 30 mn - Source Eddypro /traitements Jérôme
if (itrait_LE2_FLX_30mn == 0)
    disp('  - Traitement de LE et H à 30 mn dérivée des traitements Eddypro Jerome');

    % Initialisation
    H_Je_30mn_P1 =[];H_Je_30mn_P2 =[];H_Je_30mn_P3 =[];H_Je_30mn_P4 =[];H_Je_30mn_P5 =[];
    LE_Je_30mn_P1 =[];LE_Je_30mn_P2 =[];LE_Je_30mn_P3 =[];LE_Je_30mn_P4 =[];LE_Je_30mn_P5 =[];
    CO2_Je_30mn_P1 =[];CO2_Je_30mn_P2 =[];CO2_Je_30mn_P3 =[];CO2_Je_30mn_P4 =[];CO2_Je_30mn_P5 =[];
    idate_LE_Je_30mn_P1=[];idate_LE_Je_30mn_P2=[];idate_LE_Je_30mn_P3=[];idate_LE_Je_30mn_P4=[];idate_LE_Je_30mn_P5=[];
    
    %Entetes fichiers
    rootfile = '../../../Data/Orgeval/RAW_files/';
    filein = {'2015/FLX/Licor/raw_UTC1/Eddypro/eddypro_Jerome_full_output_2023-01-13T072323_adv.csv';...
              '2015/FLX/Licor/raw/Eddypro/eddypro_Jerome_full_output_2023-01-12T073802_adv.csv';...
              '2016/FLX/Licor/Eddypro/eddypro_Jerome_full_output_2023-01-10T151736_adv.csv';...
              '2017/FLX/Licor/S1/Eddypro/eddypro_Jerome_full_output_2023-01-11T071213_adv.csv';...
              '2017/FLX/Licor/S2/Eddypro/eddypro_Jerome_full_output_2023-01-11T161621_adv.csv';...
              '2018/FLX/Licor/Eddypro/eddypro_Jerome_full_output_2023-01-12T220651_adv.csv'};
          DataTable_Separateur=',';
    DataTable_Lentete=3;
    DataTable_Nbcol=470;
    DataTable_champ={'filename','date','time','DOY','daytime','file_records','used_records',...
                     'Tau','qc_Tau','H','qc_H','LE','qc_LE','co2_flux','qc_co2_flux','h2o_flux','qc_h2o_flux',...
                     'H_strg','LE_strg','co2_strg','h2o_strg','co2_v-adv','h2o_v-adv',...
                     'co2_molar_density','co2_mole_fraction','co2_mixing_ratio','co2_time_lag','co2_def_timelag','h2o_molar_density',...
                     'h2o_mole_fraction','h2o_mixing_ratio','h2o_time_lag','h2o_def_timelag','sonic_temperature','air_temperature',... %30-35
                     'air_pressure','air_density','air_heat_capacity','air_molar_volume','ET,water_vapor_density',...
                     'e','es','specific_humidity','RH','VPD','Tdew','u_unrot','v_unrot','w_unrot','u_rot','v_rot','w_rot','wind_speed',...
                     'max_wind_speed','wind_dir','yaw','pitch','roll','u*','TKE','L','(z-d)/L','bowen_ratio','T*','model','x_peak','x_offset',...
                     'x_10%','x_30%','x_50%','x_70%','x_90%',...
                     'un_Tau','Tau_scf','un_H','H_scf','un_LE','LE_scf','un_co2_flux','co2_scf','un_h2o_flux',... %73-
                     'h2o_scf','spikes_hf','amplitude_resolution_hf','drop_out_hf','absolute_limits_hf','skewness_kurtosis_hf','skewness_kurtosis_sf',...
                     'discontinuities_hf','discontinuities_sf','timelag_hf','timelag_sf','attack_angle_hf','non_steady_wind_hf',...
                     'u_spikes','v_spikes','w_spikes','ts_spikes','co2_spikes','h2o_spikes','chopper_LI-7500','detector_LI-7500','pll_LI-7500',...
                     'sync_LI-7500','mean_value_RSSI_LI-7500','u_var','v_var','w_var','ts_var','co2_var','h2o_var','w/ts_cov','w/co2_cov','w/h2o_cov'};
                
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        c1=char(Databrut.textdata(1+DataTable_Lentete:end,2));  %Jour
        c2=char(Databrut.textdata(1+DataTable_Lentete:end,3));  %Heure
        ijour = str2double(string(c1(:,9:10)));
        imois = str2double(string(c1(:,6:7)));
        ian   = str2double(string(c1(:,1:4)));
        iheur = str2double(string(c2(:,1:2)));
        imin  = str2double(string(c2(:,4:5)));
        idate_LE_Je_30mn=datenum(ian,imois,ijour,iheur,imin,0)-datenum(0,0,0,1,0,0);  %Decalage 1h
        
        %Colonnnes selectionnées
        colh   = 7; %7 ou 75
        colle  = 9; %9 ou 77
        colco2 = 11; %11 ou 79
        
        %Stockage data
        switch nbfile
            case(1)
                idate_LE_Je_30mn_P1 = idate_LE_Je_30mn;
                H_Je_30mn_P1  = Databrut.data(:,colh); 
                LE_Je_30mn_P1 = Databrut.data(:,colle); %9 ou 77
                CO2_Je_30mn_P1= Databrut.data(:,colco2); %11 ou 79
            case(2)
                idate_LE_Je_30mn_P2 = idate_LE_Je_30mn;
                H_Je_30mn_P2  = Databrut.data(:,colh);
                LE_Je_30mn_P2 = Databrut.data(:,colle);
                CO2_Je_30mn_P2= Databrut.data(:,colco2);
            case(3)
                idate_LE_Je_30mn_P3 = idate_LE_Je_30mn;
                H_Je_30mn_P3  = Databrut.data(:,colh);
                LE_Je_30mn_P3 = Databrut.data(:,colle);
                CO2_Je_30mn_P3= Databrut.data(:,colco2);
            case(4)
                idate_LE_Je_30mn_P4 = idate_LE_Je_30mn;
                H_Je_30mn_P4  = Databrut.data(:,colh);
                LE_Je_30mn_P4 = Databrut.data(:,colle);
                CO2_Je_30mn_P4= Databrut.data(:,colco2);
            case(5)
                idate_LE_Je_30mn_P5 = idate_LE_Je_30mn;
                H_Je_30mn_P5  = Databrut.data(:,colh);
                LE_Je_30mn_P5 = Databrut.data(:,colle);
                CO2_Je_30mn_P5= Databrut.data(:,colco2);
            case(6)
                idate_LE_Je_30mn_P6 = idate_LE_Je_30mn;
                H_Je_30mn_P6  = Databrut.data(:,colh);
                LE_Je_30mn_P6 = Databrut.data(:,colle);
                CO2_Je_30mn_P6= Databrut.data(:,colco2);
        end
    end
    
    %Concatenation annuelle
    idate_LE_Je_tmp=[idate_LE_Je_30mn_P1;idate_LE_Je_30mn_P2;idate_LE_Je_30mn_P3;idate_LE_Je_30mn_P4;idate_LE_Je_30mn_P5;idate_LE_Je_30mn_P6];
    H_tmp   = [H_Je_30mn_P1;H_Je_30mn_P2;H_Je_30mn_P3;H_Je_30mn_P4;H_Je_30mn_P5;H_Je_30mn_P6];
    LE_tmp  = [LE_Je_30mn_P1;LE_Je_30mn_P2;LE_Je_30mn_P3;LE_Je_30mn_P4;LE_Je_30mn_P5;LE_Je_30mn_P6];
    CO2_tmp = [CO2_Je_30mn_P1;CO2_Je_30mn_P2;CO2_Je_30mn_P3;CO2_Je_30mn_P4;CO2_Je_30mn_P5;CO2_Je_30mn_P6];
    
    %Comblage des lacunes
    idate_LE_Je_30mn = datenum(2015,1,1,0,0,0):datenum(0,0,0,0,30,0):datenum(2018,31,12,23,30,0);
    H_Je_30mn   = NaN*zeros(length(idate_LE_Je_30mn),1);
    LE_Je_30mn  = NaN*zeros(length(idate_LE_Je_30mn),1);
    CO2_Je_30mn = NaN*zeros(length(idate_LE_Je_30mn),1);
    for idate=1:length(idate_LE_Je_30mn)
        dump = find(idate_LE_Je_tmp == idate_LE_Je_30mn(idate));
        if (~isempty(dump))
            H_Je_30mn(idate,1)   = H_tmp(dump,1);
            LE_Je_30mn(idate,1)  = LE_tmp(dump,1);
            CO2_Je_30mn(idate,1) = CO2_tmp(dump,1);
        end
    end
    
    %Sauvegarde
    save('LE_Je_30mn.mat','H_Je_30mn','LE_Je_30mn','CO2_Je_30mn','idate_LE_Je_30mn');
    
else
    %Chargement
        disp('  - Chargement direct de LE et H à 30 mn dérivée des traitements eddypro Jerome');
    load('LE_Je_30mn.mat');
end

% Flux H et LE
if (itrait_LE_FLX_30mn == 0)
    disp('  - Traitement de H et LE FLUX à 30 mn');

    %Entetes fichiers
    rootfile = 'DataAxel/';
    filein = 'Flux_G_H_LE_Obs.csv';
    DataTable_Separateur=';';
    DataTable_Lentete=1;
    DataTable_Nbcol=6;
    DataTable_champ={'DateHeure','G1','G2','G3','LE','H'};
    
    %Lecture
    Databrut=importdata(strcat(rootfile,filein),DataTable_Separateur,DataTable_Lentete);
        
    %Date
    c=char(Databrut.textdata(:,1));
    ijour = str2double(string(c(2:end,1:2)));
    imois = str2double(string(c(2:end,4:5)));
    ian   = str2double(string(c(2:end,7:10)));
    iheur = str2double(string(c(2:end,12:13)));
    imin  = str2double(string(c(2:end,15:16)));
    idate_LE_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        
    %Stockage data
    for ivar=2:DataTable_Nbcol
        DataMat.(DataTable_champ{ivar})=Databrut.data(:,ivar-1);
    end
    H_FLX_30mn = DataMat.H;
    LE_FLX_30mn= DataMat.LE;

    %Sauvegarde
    save('LE_FLX_30mn.mat','H_FLX_30mn','LE_FLX_30mn','idate_LE_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de H et LE FLX à 30 mn dérivée de la BDOH');
    load('LE_FLX_30mn.mat');
end

% Rayonnements
if (itrait_RN_FLX_30mn == 0)
    disp('  - Traitement de RN - FLUX à 30 mn - AXEL');

    %Entetes fichiers
    rootfile = 'DataAxel/';
    filein = 'Rayonnement_Obs.csv';
    DataTable_Separateur=';';
    DataTable_Lentete=1;
    DataTable_Nbcol=7;
    DataTable_champ={'DateHeure','Sw_in','Sw_out','Lw_in_FLX','Lw_in_ANX','Lw_out','Sw_out_ANX'};
    
    %Lecture
    Databrut=importdata(strcat(rootfile,filein),DataTable_Separateur,DataTable_Lentete);
        
    %Date
    c=char(Databrut.textdata(:,1));
    ijour = str2double(string(c(2:end,1:2)));
    imois = str2double(string(c(2:end,4:5)));
    ian   = str2double(string(c(2:end,7:10)));
    iheur = str2double(string(c(2:end,12:13)));
    imin  = str2double(string(c(2:end,15:16)));
    idate_RN_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        
    %Stockage data
    for ivar=2:DataTable_Nbcol
        DataMat.(DataTable_champ{ivar})=Databrut.data(:,ivar-1);
    end
    SWin_FLX_30mn = DataMat.Sw_in;
    SWout_FLX_30mn= DataMat.Sw_out;
    toto=find(SWin_FLX_30mn<0);
    SWin_FLX_30mn(toto) = 0.0;
    SWout_FLX_30mn(toto) = 0.0;
    LWin_FLX_30mn = DataMat.Lw_in_FLX;
    toto2=find(LWin_FLX_30mn==-9999.0);
    LWin_FLX_30mn(toto2) = DataMat.Lw_in_ANX(toto2);
    LWout_FLX_30mn= DataMat.Lw_out;
    RN_FLX_30mn = SWin_FLX_30mn-SWout_FLX_30mn+LWin_FLX_30mn-LWout_FLX_30mn;
    Alb_FLX_30mn = 0.*SWout_FLX_30mn;
    toto=find(SWin_FLX_30mn>0);
    Alb_FLX_30mn(toto) = SWout_FLX_30mn(toto)./SWin_FLX_30mn(toto);
    
    %Sauvegarde
    save('RN_FLX_30mn.mat','RN_FLX_30mn','Alb_FLX_30mn','SWin_FLX_30mn','SWout_FLX_30mn','LWin_FLX_30mn','LWout_FLX_30mn','idate_RN_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de RN FLX à 30 mn');
    load('RN_FLX_30mn.mat');
end

%% Flux Rayonnement TAILLIS/FLX à 30 mn - Source BDOH
if (itrait_Ray_FLX_30mn == 0)
    disp('  - Traitement des rayonnements FLUX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_FLX/';
    filein = {'TAILLIS_PRAY.txt';'TAILLIS_PRAY-2.txt';'TAILLIS_PRAY-3.txt';'TAILLIS_PRAY-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_Ray_FLX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SWin_FLX_30mn = Databrut.data(:,1);
            case(2)
                SWout_FLX_30mn = Databrut.data(:,1);
            case(3)
                LWin_FLX_30mn = Databrut.data(:,1);
            case(4)
                LWout_FLX_30mn = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('Ray_FLX_30mn.mat','SWin_FLX_30mn','SWout_FLX_30mn','LWin_FLX_30mn','LWout_FLX_30mn','idate_Ray_FLX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de Rayonnement FLX à 30 mn dérivée de la BDOH');
    load('Ray_FLX_30mn.mat');
end

%% Flux Rayonnement TAILLIS/FLX à 30 mn - Source BDOH
if (itrait_Ray_ANX_30mn == 0)
    disp('  - Traitement des rayonnements ANX à 30 mn dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_ANX/';
    filein = {'NOUE_PRAY.txt';'NOUE_PRAY-2.txt';'NOUE_PRAY-3.txt';'NOUE_PRAY-4.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_Ray_ANX_30mn=datenum(ian,imois,ijour,iheur,imin,0);
        end
        
        %Stockage data
        switch nbfile
            case(1)
                SWin_ANX_30mn = Databrut.data(:,1);
            case(2)
                SWout_ANX_30mn = Databrut.data(:,1);
            case(3)
                LWin_ANX_30mn = Databrut.data(:,1);
            case(4)
                LWout_ANX_30mn = Databrut.data(:,1);
        end
    end
    
    %Sauvegarde
    save('Ray_ANX_30mn.mat','SWin_ANX_30mn','SWout_ANX_30mn','LWin_ANX_30mn','LWout_ANX_30mn','idate_Ray_ANX_30mn');
    
else
    %Chargement
    disp('  - Chargement direct de Rayonnement FLX à 30 mn dérivée de la BDOH');
    load('Ray_ANX_30mn.mat');
end

%% Flux Rayonnement TAILLIS/FLX à 30 mn - Source BDOH
if (itrait_Ray_BOI_1h == 0)
    disp('  - Traitement des rayonnements BOISSY à 1h dérivée de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_Boissy/';
    filein = {'BOISSY-METEO_PRAY.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};
    
    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        if (nbfile == 1)
            c=char(Databrut.textdata(:,1));
            ijour = str2double(string(c(4:end,1:2)));
            imois = str2double(string(c(4:end,4:5)));
            ian   = str2double(string(c(4:end,7:10)));
            iheur = str2double(string(c(4:end,12:13)));
            imin  = str2double(string(c(4:end,15:16)));
            idate_Ray_BOISSY_1h=datenum(ian,imois,ijour,iheur,imin,0)+datenum(0,0,0,0,30,0);
        end
        
        %Stockage data
        SWin_BOISSY_1h = Databrut.data(:,1);

    end
    
    %Mise à NaN des données bidons
    dump=SWin_BOISSY_1h==-9999.0;
    SWin_BOISSY_1h(dump)=NaN;
    
    %Conversion en W/m²
    SWin_BOISSY_1h = (10000/3600)*SWin_BOISSY_1h;
    
    %Sauvegarde
    save('Ray_BOISSY_1h.mat','idate_Ray_BOISSY_1h','SWin_BOISSY_1h');
    
else
    %Chargement
    disp('  - Chargement direct de Rayonnement BOISSY à 1h n dérivée de la BDOH');
    load('Ray_BOISSY_1h.mat');
end

%Précipitations
if (itrait_Pluie_1h == 0)
    disp('  - Traitement des précipitations dérivées de la BDOH');

    %Entetes fichiers
    rootfile = 'DataBDOH_pluie/';
    filein = {'BOISSY-METEO_PRCP.txt';'BOISSY-P28_PRCP.txt';'LOGE-P07_PRCP.txt';'MELARCHEZ-P35_PRCP.txt'};
    DataTable_Separateur=';';
    DataTable_Lentete=3;
    DataTable_Nbcol=1;
    DataTable_champ={'DateHeure','Valeur','Qualite','Min','Max'};

    %Lecture
    for nbfile = 1:length(filein)
        Databrut=importdata(strcat(rootfile,filein{nbfile}),DataTable_Separateur,DataTable_Lentete);
        
        %Date
        c=char(Databrut.textdata(:,1));
        ijour = str2double(string(c(4:end,1:2)));
        imois = str2double(string(c(4:end,4:5)));
        ian   = str2double(string(c(4:end,7:10)));
        iheur = str2double(string(c(4:end,12:13)));
        imin  = str2double(string(c(4:end,15:16)));
        idate_pluie=datenum(ian,imois,ijour,iheur,imin,0);
        
        %Stockage data
        switch nbfile
            case(1)
                idate_pr_Boissy_1h = idate_pluie;
                Pr_Boissy_1h = Databrut.data(:,1);
            case(2)
                idate_pr_BoissyP28_1h = idate_pluie;
                Pr_BoissyP28_1h = Databrut.data(:,1);
            case(3)
                idate_pr_Loge_1h = idate_pluie;
                Pr_Loge_1h = Databrut.data(:,1);
            case(4)
                idate_pr_Melarchez_1h = idate_pluie;
                Pr_Melarchez_1h = Databrut.data(:,1);
        end
    end
    
    dump=(Pr_Melarchez_1h<0) | (Pr_Melarchez_1h>500);
    Pr_Melarchez_1h(dump)=NaN;
    dump=(Pr_Loge_1h<0) | (Pr_Loge_1h>500);
    Pr_Loge_1h(dump)=NaN;
    dump=(Pr_BoissyP28_1h<0) | (Pr_BoissyP28_1h>500);
    Pr_BoissyP28_1h(dump)=NaN;
    dump=(Pr_Boissy_1h<0) | (Pr_Boissy_1h>500);
    Pr_Boissy_1h(dump)=NaN;
    
    %Sauvegarde
    save('Precip_BV_1h.mat','idate_pr_Melarchez_1h','Pr_Melarchez_1h','idate_pr_Loge_1h','Pr_Loge_1h','idate_pr_BoissyP28_1h','Pr_BoissyP28_1h','idate_pr_Boissy_1h','Pr_Boissy_1h');

else
    %Chargement
    disp('  - Chargement des précipitations dérivées de la BDOH');
    load('Precip_BV_1h.mat');
end

%% Calculs intermédiaires

%Conversion des périodes en humidités par toutes les relations Campbell
SMCq1_FLX_1mn_P1 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_1mn_P1 + RC_q(3).*(SMCp_FLX_1mn_P1.*SMCp_FLX_1mn_P1));
SMCq1_FLX_1mn_P2 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_1mn_P2 + RC_q(3).*(SMCp_FLX_1mn_P2.*SMCp_FLX_1mn_P2));
SMCq1_FLX_1mn_P3 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_1mn_P3 + RC_q(3).*(SMCp_FLX_1mn_P3.*SMCp_FLX_1mn_P3));
SMCq1_FLX_1mn_P4 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_1mn_P4 + RC_q(3).*(SMCp_FLX_1mn_P4.*SMCp_FLX_1mn_P4));
SMCl1_FLX_1mn_P1 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_1mn_P1);
SMCl1_FLX_1mn_P2 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_1mn_P2);
SMCl1_FLX_1mn_P3 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_1mn_P3);
SMCl1_FLX_1mn_P4 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_1mn_P4);

SMCq1_FLX_30mn_P1 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_30mn_P1 + RC_q(3).*(SMCp_FLX_30mn_P1.*SMCp_FLX_30mn_P1));
SMCq1_FLX_30mn_P2 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_30mn_P2 + RC_q(3).*(SMCp_FLX_30mn_P2.*SMCp_FLX_30mn_P2));
SMCq1_FLX_30mn_P3 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_30mn_P3 + RC_q(3).*(SMCp_FLX_30mn_P3.*SMCp_FLX_30mn_P3));
SMCq1_FLX_30mn_P4 = 100*(RC_q(1)+ RC_q(2).*SMCp_FLX_30mn_P4 + RC_q(3).*(SMCp_FLX_30mn_P4.*SMCp_FLX_30mn_P4));
SMCl1_FLX_30mn_P1 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_30mn_P1);
SMCl1_FLX_30mn_P2 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_30mn_P2);
SMCl1_FLX_30mn_P3 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_30mn_P3);
SMCl1_FLX_30mn_P4 = 100*(RC_l(1)+ RC_l(2).*SMCp_FLX_30mn_P4);

SMCq2_FLX_1mn_P1 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_1mn_P1 + RM1_q(3).*(SMCp_FLX_1mn_P1.*SMCp_FLX_1mn_P1)+Bias1_q(1));
SMCq2_FLX_1mn_P2 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_1mn_P2 + RM1_q(3).*(SMCp_FLX_1mn_P2.*SMCp_FLX_1mn_P2)+Bias1_q(2));
SMCq2_FLX_1mn_P3 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_1mn_P3 + RM1_q(3).*(SMCp_FLX_1mn_P3.*SMCp_FLX_1mn_P3)+Bias1_q(3));
SMCq2_FLX_1mn_P4 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_1mn_P4 + RM1_q(3).*(SMCp_FLX_1mn_P4.*SMCp_FLX_1mn_P4)+Bias1_q(4));
SMCl2_FLX_1mn_P1 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_1mn_P1+Bias1_l(1));
SMCl2_FLX_1mn_P2 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_1mn_P2+Bias1_l(2));
SMCl2_FLX_1mn_P3 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_1mn_P3+Bias1_l(3));
SMCl2_FLX_1mn_P4 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_1mn_P4+Bias1_l(4));

SMCq2_FLX_30mn_P1 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_30mn_P1 + RM1_q(3).*(SMCp_FLX_30mn_P1.*SMCp_FLX_30mn_P1)+Bias1_q(1));
SMCq2_FLX_30mn_P2 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_30mn_P2 + RM1_q(3).*(SMCp_FLX_30mn_P2.*SMCp_FLX_30mn_P2)+Bias1_q(2));
SMCq2_FLX_30mn_P3 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_30mn_P3 + RM1_q(3).*(SMCp_FLX_30mn_P3.*SMCp_FLX_30mn_P3)+Bias1_q(3));
SMCq2_FLX_30mn_P4 = 100*(RM1_q(1)+ RM1_q(2).*SMCp_FLX_30mn_P4 + RM1_q(3).*(SMCp_FLX_30mn_P4.*SMCp_FLX_30mn_P4)+Bias1_q(4));
SMCl2_FLX_30mn_P1 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_30mn_P1+Bias1_l(1));
SMCl2_FLX_30mn_P2 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_30mn_P2+Bias1_l(2));
SMCl2_FLX_30mn_P3 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_30mn_P3+Bias1_l(3));
SMCl2_FLX_30mn_P4 = 100*(RM1_l(1)+ RM1_l(2).*SMCp_FLX_30mn_P4+Bias1_l(4));

SMCq3_FLX_1mn_P1 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_1mn_P1 + RM2_q(3).*(SMCp_FLX_1mn_P1.*SMCp_FLX_1mn_P1)+Bias2_q(1));
SMCq3_FLX_1mn_P2 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_1mn_P2 + RM2_q(3).*(SMCp_FLX_1mn_P2.*SMCp_FLX_1mn_P2)+Bias2_q(2));
SMCq3_FLX_1mn_P3 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_1mn_P3 + RM2_q(3).*(SMCp_FLX_1mn_P3.*SMCp_FLX_1mn_P3)+Bias2_q(3));
SMCq3_FLX_1mn_P4 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_1mn_P4 + RM2_q(3).*(SMCp_FLX_1mn_P4.*SMCp_FLX_1mn_P4)+Bias2_q(4));
SMCl3_FLX_1mn_P1 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_1mn_P1+Bias2_l(1));
SMCl3_FLX_1mn_P2 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_1mn_P2+Bias2_l(2));
SMCl3_FLX_1mn_P3 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_1mn_P3+Bias2_l(3));
SMCl3_FLX_1mn_P4 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_1mn_P4+Bias2_l(4));

SMCq3_FLX_30mn_P1 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_30mn_P1 + RM2_q(3).*(SMCp_FLX_30mn_P1.*SMCp_FLX_30mn_P1)+Bias2_q(1));
SMCq3_FLX_30mn_P2 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_30mn_P2 + RM2_q(3).*(SMCp_FLX_30mn_P2.*SMCp_FLX_30mn_P2)+Bias2_q(2));
SMCq3_FLX_30mn_P3 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_30mn_P3 + RM2_q(3).*(SMCp_FLX_30mn_P3.*SMCp_FLX_30mn_P3)+Bias2_q(3));
SMCq3_FLX_30mn_P4 = 100*(RM2_q(1)+ RM2_q(2).*SMCp_FLX_30mn_P4 + RM2_q(3).*(SMCp_FLX_30mn_P4.*SMCp_FLX_30mn_P4)+Bias2_q(4));
SMCl3_FLX_30mn_P1 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_30mn_P1+Bias2_l(1));
SMCl3_FLX_30mn_P2 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_30mn_P2+Bias2_l(2));
SMCl3_FLX_30mn_P3 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_30mn_P3+Bias2_l(3));
SMCl3_FLX_30mn_P4 = 100*(RM2_l(1)+ RM2_l(2).*SMCp_FLX_30mn_P4+Bias2_l(4));

% SMCq4_FLX_1mn_P1 = 100*(RM3_q(1)+ RM3_q(2).*(SMCp_FLX_1mn_P1/1000) + RM3_q(3).*((SMCp_FLX_1mn_P1/1000).*(SMCp_FLX_1mn_P1/1000)));
% SMCq4_FLX_1mn_P2 = 100*(RM3_q(1)+ RM3_q(2).*(SMCp_FLX_1mn_P2/1000) + RM3_q(3).*((SMCp_FLX_1mn_P2/1000).*(SMCp_FLX_1mn_P2/1000)));
% SMCq4_FLX_1mn_P3 = 100*(RM3_q(1)+ RM3_q(2).*(SMCp_FLX_1mn_P3/1000) + RM3_q(3).*((SMCp_FLX_1mn_P3/1000).*(SMCp_FLX_1mn_P3/1000)));
% SMCq4_FLX_1mn_P4 = 100*(RM3_q(1)+ RM3_q(2).*(SMCp_FLX_1mn_P4/1000) + RM3_q(3).*((SMCp_FLX_1mn_P4/1000).*(SMCp_FLX_1mn_P4/1000)));
% 

%% Tracés
disp('  - Phase de tracé');

%Periodes à 1 minute
% figure(1)
% hold on
% plot(idate_SMC_FLX_1mn,SMCp_FLX_1mn_P1,'g.');
% plot(idate_SMC_FLX_1mn,SMCp_FLX_1mn_P2,'c.');
% plot(idate_SMC_FLX_1mn,SMCp_FLX_1mn_P3,'b.');
% plot(idate_SMC_FLX_1mn,SMCp_FLX_1mn_P4,'r.');
% datetick('x','QQ-YY')

% %SMC BOISSY
% figure(2)
% subplot(1,2,1)
% hold on
% plot(idate_SMC_Boissyaxel,SMC_Boissyaxel_P1,'g.');
% plot(idate_SMC_Boissyaxel,SMC_Boissyaxel_P2,'c.');
% plot(idate_SMC_Boissyaxel,SMC_Boissyaxel_P3,'b.');
% plot(idate_SMC_Boissyaxel,SMC_Boissyaxel_P4,'r.');
% datetick('x','QQ-YY')
% set(gca,'ylim',[15,65],'xlim',[idate_SMC_FLX_1mn(1),idate_SMC_FLX_1mn(end)]);
% subplot(1,2,2)
% hold on
% plot(idate_SMC_Boissy,SMC_Boissy_P1,'g.');
% plot(idate_SMC_Boissy,SMC_Boissy_P2,'c.');
% plot(idate_SMC_Boissy,SMC_Boissy_P3,'b.');
% plot(idate_SMC_Boissy,SMC_Boissy_P4,'r.');
% datetick('x','QQ-YY')
% set(gca,'ylim',[15,65],'xlim',[idate_SMC_FLX_1mn(1),idate_SMC_FLX_1mn(end)]);
% 
% %SMC Suizy
% figure(3)
% hold on
% plot(idate_SMC_SUIZY_12h,SMC_SUIZY_12h_P1,'g.');
% plot(idate_SMC_SUIZY_12h,SMC_SUIZY_12h_P2,'c.');
% plot(idate_SMC_SUIZY_12h,SMC_SUIZY_12h_P3,'b.');
% plot(idate_SMC_SUIZY_12h,SMC_SUIZY_12h_P4,'r.');
% datetick('x','QQ-YY')
% set(gca,'ylim',[15,65],'xlim',[idate_SMC_FLX_1mn(1),idate_SMC_FLX_1mn(end)]);

%SMC ANX
figure(4)
hold on
hp=plot(idate_SMC_ANX_30mn,SMC_ANX_30mn_P1,'go');
set(hp,'markersize',1,'markerfacecolor','g');
hp=plot(idate_SMC_ANX_30mn,SMC_ANX_30mn_P2,'co');
set(hp,'markersize',1,'markerfacecolor','c');
hp=plot(idate_SMC_ANX_30mn,SMC_ANX_30mn_P4,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_ANX_30mn,SMC_ANX_30mn_P3,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
ti=title('SMC - Station ANX');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
datetick('x','QQ-YY')
set(gca,'ylim',[15,65],'xlim',[idate_SMC_ANX_30mn(1)-1,idate_SMC_ANX_30mn(end)]);

%Periodes à 30 minutes
figure(10)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCp_FLX_30mn_P1,'go');
set(hp,'markersize',1,'markerfacecolor','g');
hp=plot(idate_SMC_FLX_30mn,SMCp_FLX_30mn_P2,'co');
set(hp,'markersize',1,'markerfacecolor','c');
hp=plot(idate_SMC_FLX_30mn,SMCp_FLX_30mn_P3,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_FLX_30mn,SMCp_FLX_30mn_P4,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
ti=title('Periode station FLUX ');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('Periode (µs)');
le = legend('FLX 7cm','FLX 23cm','FLX 50cm','FLX 80cm');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[10,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

% SMC à 30 minute
figure(5)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P1,'go');
set(hp,'markersize',1,'markerfacecolor','g');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P2,'co');
set(hp,'markersize',1,'markerfacecolor','c');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P3,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P4,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
ti=title('SMC station FLUX ');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
le = legend('FLX 7cm','FLX 23cm','FLX 50cm','FLX 80cm');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[10,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

% SMC 10cm multisource
figure(6)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCq2_FLX_30mn_P1,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P1,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_Boissy,SMC_Boissy_P1,'ko');
set(hp,'markersize',1,'markerfacecolor','k');
ti=title('SMC - 10cm');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
set(y,'fontsize',12,'fontweight','bold','fontname','arial');
le = legend('RQ_C_a_m_b_e_l_l v1-debias','RQ_C_a_m_b_e_l_l v3-debias','Station BOISSY','location','southwest');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[0,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

% SMC 25cm multisource
figure(7)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCq2_FLX_30mn_P2,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P2,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_Boissy,SMC_Boissy_P2,'ko');
set(hp,'markersize',1,'markerfacecolor','k');
ti=title('SMC - 25cm');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
set(y,'fontsize',12,'fontweight','bold','fontname','arial');
le = legend('RQ_C_a_m_b_e_l_l v1-debias','RQ_C_a_m_b_e_l_l v3-debias','Station BOISSY','location','southwest');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[0,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

% SMC 50cm multisource
figure(8)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCq2_FLX_30mn_P3,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P3,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_Boissy,SMC_Boissy_P3,'ko');
set(hp,'markersize',1,'markerfacecolor','k');
ti=title('SMC - 50cm');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
set(y,'fontsize',12,'fontweight','bold','fontname','arial');
le = legend('RQ_C_a_m_b_e_l_l v1-debias','RQ_C_a_m_b_e_l_l v3-debias','Station BOISSY','location','southwest');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[10,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

% SMC 80cm multisource
figure(9)
hold on
hp=plot(idate_SMC_FLX_30mn,SMCq2_FLX_30mn_P4,'ro');
set(hp,'markersize',1,'markerfacecolor','r');
hp=plot(idate_SMC_FLX_30mn,SMCq3_FLX_30mn_P4,'bo');
set(hp,'markersize',1,'markerfacecolor','b');
hp=plot(idate_SMC_Boissy,SMC_Boissy_P4,'ko');
set(hp,'markersize',1,'markerfacecolor','k');
ti=title('SMC - 80cm');
set(ti,'fontsize',12,'fontweight','bold','fontname','arial');
x=xlabel('Time');
set(x,'fontsize',10,'fontweight','bold','fontname','arial');
y=ylabel('SMC (%)');
set(y,'fontsize',12,'fontweight','bold','fontname','arial');
le = legend('RQ_C_a_m_b_e_l_l v1-debias','RQ_C_a_m_b_e_l_l v3-debias','Station BOISSY','location','southwest');
set(le,'fontsize',10,'fontweight','bold','fontname','arial');
datetick('x','QQ-YY');
set(gca,'ylim',[10,50],'xlim',[idate_SMC_FLX_30mn(2),idate_SMC_FLX_30mn(end)]);

%Save
save('SMC_Orgeval.mat','idate_SMC_FLX_1mn','SMCq1_FLX_1mn_P1','SMCq2_FLX_1mn_P2','SMCq3_FLX_1mn_P3','SMCl3_FLX_1mn_P3','SMCq3_FLX_1mn_P4','SMCl3_FLX_1mn_P4');
save('SMC_Orgeval_30mn.mat','idate_SMC_FLX_30mn','SMCq2_FLX_30mn_P1','SMCq3_FLX_30mn_P1','SMCq2_FLX_30mn_P2','SMCq3_FLX_30mn_P2','SMCq2_FLX_30mn_P3','SMCq3_FLX_30mn_P3','SMCq2_FLX_30mn_P4','SMCq3_FLX_30mn_P4');


% %% Preparation des fichiers T et W fixées au fond
% %Ecriture 2016_2018
% disp('  - Ecriture du fichier de forcage Humidité à 80cm sur la période 2016-2018');
% file_out='Simulations/Input/SMC80_2016_2018.txt';
% fid1 = fopen(file_out,'w+');
% v = datevec(idate_SMC_FLX_30mn);
% dump = find((v(:,1) >= 2016) & (v(:,4) == 12) & (v(:,5) == 0));
% ijour = (1:1:length(dump))';
% toto = isfinite(SMCq3_FLX_30mn_P4(dump));
% vh_16 = [ijour(toto) SMCq3_FLX_30mn_P4(1,dump(toto))'/100];
% a=1;
% fprintf(fid1,'%i\n',a);
% fprintf(fid1,'%i\n',size(vh_16,1));
% for i=1:size(vh_16,1)
%     fprintf(fid1,'%i %6.3f\n',vh_16(i,:));
% end
% fclose(fid1);
% disp('  - Ecriture du fichier de forcage Temperature à 80cm sur la période 2016-2018');
% file_out='Simulations/Input/T80_2016_2018.txt';
% fid1 = fopen(file_out,'w+');
% v = datevec(idate_Tsol_FLX_30mn);
% dump = find((v(:,1) >= 2016) & (v(:,4) == 12) & (v(:,5) == 0));
% ijour = (1:1:length(dump))';
% toto = isfinite(Tsol_FLX_30mn_P4(dump));
% vt_16 = [ijour(toto) Tsol_FLX_30mn_P4(dump(toto))];
% fprintf(fid1,'%i\n',size(vt_16,1));
% for i=1:size(vt_16,1)
%     fprintf(fid1,'%i %6.3f\n',vt_16(i,:));
% end
% fclose(fid1);
% 
% %Ecriture 2015_2018
% disp('  - Ecriture du fichier de forcage Humidité à 80cm sur la période 2015-2018');
% file_out='Simulations/Input/SMC80_2015_2018.txt';
% fid1 = fopen(file_out,'w+');
% v = datevec(idate_SMC_FLX_30mn);
% dump = find((v(:,1) >= 2015) & (v(:,4) == 12) & (v(:,5) == 0));
% ijour = (1:1:length(dump))';
% toto = isfinite(SMCq3_FLX_30mn_P4(dump));
% vh_15 = [ijour(toto) SMCq3_FLX_30mn_P4(1,dump(toto))'/100];
% a=1;
% fprintf(fid1,'%i\n',a);
% fprintf(fid1,'%i\n',size(vh_15,1));
% for i=1:size(vh_15,1)
%     fprintf(fid1,'%i %6.3f\n',vh_15(i,:));
% end
% fclose(fid1);
% 
% disp('  - Ecriture du fichier de forcage Temperature à 80cm sur la période 2015-2018');
% file_out='Simulations/Input/T80_2015_2018.txt';
% fid1 = fopen(file_out,'w+');
% v = datevec(idate_Tsol_FLX_30mn);
% dump = find((v(:,1) >= 2015) & (v(:,4) == 12) & (v(:,5) == 0));
% ijour = (1:1:length(dump))';
% toto = isfinite(Tsol_FLX_30mn_P4(dump));
% vt_15 = [ijour(toto) Tsol_FLX_30mn_P4(dump(toto))];
% fprintf(fid1,'%i\n',size(vt_15,1));
% for i=1:size(vt_15,1)
%     fprintf(fid1,'%i %6.3f\n',vt_15(i,:));
% end
% fclose(fid1);
% 
%% Preparation de la condition initiale 2016
%Ecriture 2016_2018
disp('  - Ecriture du fichier de condition initiale Humidité à 80cm sur la période 2016-2018');
file_out='Simulations/Input/CI80_2016_2018.txt';
fid1 = fopen(file_out,'w+');
nb_noeud = (1:1:130)';
moeud_obs = [1 28 60 92 130];
v = datevec(idate_SMC_FLX_30mn);
dump = find((v(:,1) == 2016) & (v(:,2) == 1) & (v(:,3) == 1) & (v(:,4) == 0) & (v(:,5) == 0));
SMC_obs = [SMCq3_FLX_30mn_P1(dump) SMCq3_FLX_30mn_P1(dump) SMCq3_FLX_30mn_P2(dump) SMCq3_FLX_30mn_P3(dump) SMCq3_FLX_30mn_P4(dump)];
v = datevec(idate_SMC_FLX_30mn);
dump = find((v(:,1) == 2016) & (v(:,2) == 1) & (v(:,3) == 1) & (v(:,4) == 0) & (v(:,5) == 0));
T_obs   = [Tsol_FLX_30mn_P1(dump) Tsol_FLX_30mn_P1(dump) Tsol_FLX_30mn_P2(dump) Tsol_FLX_30mn_P3(dump) Tsol_FLX_30mn_P4(dump)];
SMC_ini = interp1(moeud_obs,SMC_obs,nb_noeud');
T_ini   = interp1(moeud_obs,T_obs,nb_noeud');
info_ci = [2 4 1]; 
fprintf(fid1,'%i %i %i\n',info_ci(:));
for i=1:length(T_ini)
        fprintf(fid1,'%6.2f\n',T_ini(i));
end
for i=1:length(SMC_ini)
        fprintf(fid1,'%6.3f\n',SMC_ini(i)/100);
end
fclose(fid1);

disp('  - Ecriture du fichier de condition initiale Humidité à 80cm sur la période 2015-2018');
file_out='Simulations/Input/CI80_2015_2018.txt';
fid1 = fopen(file_out,'w+');
nb_noeud = (1:1:130)';
moeud_obs = [1 28 60 92 130];
SMC_obs = [0.40 0.38 0.36 0.36 SMCq3_FLX_30mn_P4(dump)/100+0.2];
v = datevec(idate_SMC_FLX_30mn);
dump = find((v(:,1) == 2016) & (v(:,2) == 1) & (v(:,3) == 1) & (v(:,4) == 0) & (v(:,5) == 0));
T_obs   = [Tsol_FLX_30mn_P1(dump) Tsol_FLX_30mn_P1(dump) Tsol_FLX_30mn_P2(dump) Tsol_FLX_30mn_P3(dump) Tsol_FLX_30mn_P4(dump)];
SMC_ini = interp1(moeud_obs,SMC_obs,nb_noeud');
T_ini   = interp1(moeud_obs,T_obs,nb_noeud');
info_ci = [2 4 1]; 
fprintf(fid1,'%i %i %i\n',info_ci(:));
for i=1:length(T_ini)
        fprintf(fid1,'%6.2f\n',T_ini(i));
end
for i=1:length(SMC_ini)
        fprintf(fid1,'%6.3f\n',SMC_ini(i));
end
fclose(fid1);





