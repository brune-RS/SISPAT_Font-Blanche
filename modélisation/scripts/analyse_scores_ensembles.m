%% Analyse des scores d'ensembles de simulation
clear;
close all;
clc;
%addpath('/home/brune/SISPAT/Fonctions matlab/')
addpath('D:/Fonctions matlab/SAFE/RSA/')
addpath('D:/Fonctions matlab/SAFE/visualization/')
addpath ('D:\Fonctions matlab')
%quelle année considérée?
annee=2014;
% annee=2017;
Lv=2.484*10^6;
dt=1800;
% afficher figures oui=1   non=0
fig_ET=0;

%% % fichiers à charger: structures de scores dans chaque dossier 
 %path_out='/media/brune/Ultra Touch/SISPAT/out5';
 path_out='D:/SISPAT/out6';
load([path_out '/' 'scores_ensemble6.mat']);
load([path_out '/' 'scores_ET6.mat']);

%matparam  matrice valeur des paramètres pour chque simu
% indice à récupérer 
root0='/home/brune/SISPAT/execution/';
root1= [root0 'param/ensemble/ensemble_6/matparam.mat'];
root1='D:/ensemble_6/matparam.mat';
mat_param=load(root1);
mat_param.matparam(:,13)=exp(mat_param.matparam(:,13).*log(10));


% 1) prof fractures min  prof fractures max  2) fraction pluie reprise fraction ruissellement   3) porosité   4)Wsat  5)hg 
% % 6)q   7)Ksat   8) beta  9)ctherm  10)FDR max 11)%pmr 12)res stomatique
variables=[" prof. fractures " "fraction reprise "   "porosité"   "Wsat"  "hg" "q"   "Ksat"   "beta"  "ctherm"  "FDR max" "%pmr" 'res.st ' 'res.tot.'];
numsim=2000;
% choisir espace de scores à afficher 
% selon les paramètres changés 

% abscisse valeur du paramètre (indice num simu et matparam) 
% ordonnée valeur du score de chaque simu 

% quel indice stat ? rmse, r2, coeff de correlation...
var=["LE", "RN","H","W5","W25","W50","drainage"];
stat_temp=[0];
indic_stat=["kge" "rmse" "cor2" "mbias"] ;
for i=1:length(indic_stat)
    for isim=1:numsim
stat_temp(isim,1)=scores.LE.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,2)=scores.RN.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,3)=scores.H.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,4)=scores.W5.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,5)=scores.W25.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,6)=scores.W50.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
stat_temp(isim,7)=scores.drainage.(['sim_' num2str(isim)]);
ET_jour_temp(:,isim)=scores.ET_j.(['sim_' num2str(isim)]).([string(indic_stat(i))]);
ET_mean_mens(:,isim)=ET_scores.mois_obs_sim_mensuel.(['sim_' num2str(isim)]).mean;
ET_max_mens(:,isim)=ET_scores.mois_obs_sim_mensuel.(['sim_' num2str(isim)]).max;
ET_min_mens(:,isim)=ET_scores.mois_obs_sim_mensuel.(['sim_' num2str(isim)]).min;
end 
stats.([string(indic_stat(i))])=stat_temp;
ET_jour.([string(indic_stat(i))])=ET_jour_temp;
end 
%% 
is=2; %kge rmse corélation
j=7; %"LE", "RN","H","W5","W25","W50", "drainage"

% LE 
figure()
scatter3(mat_param.matparam(:,1),mat_param.matparam(:,13),mat_param.matparam(:,2),[],stats.([string(indic_stat(is))])(:,j),'filled')
colorbar
xlabel(variables(1))
ylabel(variables(13))
zlabel(variables(2))
title([indic_stat(is) var(j)])
set(gca, 'YScale', 'log')
colormap turbo

% ET journalier
figure()
scatter3(mat_param.matparam(:,13),mat_param.matparam(:,2),mat_param.matparam(:,1),[],ET_jour.([string(indic_stat(is))]),'filled')
colorbar
xlabel(variables(13))
ylabel(variables(2))
zlabel(variables(1))
title([indic_stat(is) ' ET jour'])
set(gca, 'XScale', 'log')
colormap turbo

% LE 2d
figure()
%surf(mat_param.matparam(:,10),mat_param.matparam(:,13),stat_temp(:,1),stat_temp(:,1))
scatter3(mat_param.matparam(:,1),mat_param.matparam(:,13),stats.([string(indic_stat(is))])(:,j),[],stats.([string(indic_stat(is))])(:,j),'filled')
colorbar
xlabel(variables(1))
ylabel(variables(13))
%zlabel([indic_stat var(j)])
title([indic_stat(is) var(j)])
set(gca, 'YScale', 'log')
colormap turbo
%____________________________________________________________________

% une variable Vs le score 
figure()
scatter(mat_param.matparam(:,1),stat_temp(:,1))
ylabel([indic_stat '_LE'])
xlabel(variables(1))

% une variable VS une autre 
figure()
scatter(mat_param.matparam(:,10),mat_param.matparam(:,12),[],stat_temp(:,1),'filled')
colorbar
xlabel(variables(10))
ylabel(variables(12))
title([indic_stat var])




% ET journalier
figure()
scatter3(mat_param.matparam(:,13),mat_param.matparam(:,2),mat_param.matparam(:,10),'filled')
colorbar
xlabel(variables(13))
ylabel(variables(2))
zlabel(variables(10))
title([indic_stat(is) ' ET jour'])
set(gca, 'XScale', 'log')
colormap turbo





% H
figure()
scatter3(mat_param.matparam(:,10),mat_param.matparam(:,13),stat_temp(:,3),'filled')
colorbar
xlabel(variables(10))
ylabel(variables(13))
zlabel(variables(3))
%zlabel([indic_stat 'H'])
set(gca, 'YScale', 'log')
colormap turbo


%
%% DIFFERENCE ET MENSUEL 
is=1; %kge rmse corélation
j=7; %"LE", "RN","H","W5","W25","W50", "drainage"

%month=ET_scores.mois_obs_sim_mensuel.(['sim_' num2str(isim)]);
month=["Jan."	"Feb." " Mar."  "Apr."  "May" "Jun." "Jul." "Aug." "Sep." "Oct." "Nov." "Dec."]
ax=figure()
title('ET obs-sim mensuel moyenne')
for jj=1:12
    subplot(3,4,jj)
    scatter(mat_param.matparam(:,2),ET_mean_mens(jj,:),[],somme_kge,'filled')
    
    %scatter(mat_param.matparam(list,13),ET_mean_mens(jj,list),[],stats.([string(indic_stat(is))])(list,7),'filled')
%set(gca, 'XScale', 'log')
yline(0)
ylim([-40 30])
%xlabel(variables(13))
%ylabel(['ET'])
xl=xlim;
yl=ylim;
title(month(jj),'Position',[mean(xl) yl(2)-8])
grid on
end
colormap turbo
h = axes(ax,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
xlabel(h,variables(2),'FontWeight','bold');
ylabel(h,'delta ET mensuel','FontWeight','bold');
title(h,'somme kge');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
%ylabel(c,[indic_stat(is) ' ET '],'FontSize',12,'FontWeight','bold', 'Rotation',0,'Position',[c.Position(1) max(stat_temp(1,:))+2])
colormap(c,'jet')
%caxis(h,[min(stats.([string(indic_stat(is))])(:,j)),max(stats.([string(indic_stat(is))])(:,j))]);             % set colorbar limits
caxis(h,[min(somme_kge),max(somme_kge)]);             % set colorbar limits






%% SAFE toolbox test (Pianosi)Perform Regional Sensitivity Analysis (RSA) based on grouping
var=["LE", "RN","H","W5","W25","W50","drainage"];
indic_stat=["kge" "rmse" "cor2" "mbias"] ;


params= [" prof. fractures " "fraction reprise "   "FDR max"  'res.st ' 'res.tot.'];
X(:,1)=mat_param.matparam(:,1);
X(:,2)=mat_param.matparam(:,2);
X(:,3)=mat_param.matparam(:,10);
X(:,4)=mat_param.matparam(:,12);
X(:,5)=mat_param.matparam(:,13);

% Here the following output metrics are used :

 Y(:,1)=stats.([string(indic_stat(2))])(:,1); % RMSE LE;
 Y(:,2)=stats.([string(indic_stat(4))])(:,1); % bias LE;
 Y(:,3)=stats.([string(indic_stat(1))])(:,1); % KGE LE;

 Y(:,4)=stats.([string(indic_stat(2))])(:,5); % RMSE W25;
 Y(:,5)=stats.([string(indic_stat(4))])(:,5); % bias W25;
 Y(:,6)=stats.([string(indic_stat(1))])(:,5); % KGE W25

 Y(:,7)=stats.([string(indic_stat(2))])(:,6); % RMSE W50
Y(:,8)=stats.([string(indic_stat(4))])(:,6); % biais W50
 Y(:,9)=stats.([string(indic_stat(1))])(:,6); % KGE W50

Y(:,10)=stats.([string(indic_stat(2))])(:,2); % RMSE RN
Y(:,11)=stats.([string(indic_stat(4))])(:,2); % biais RN
Y(:,12)=stats.([string(indic_stat(1))])(:,2); % KGE RN

Y(:,10)=stat_temp(:,7); %drainage
%Y(:,14)=ET_mean_mens(5,:);% mai
Y(:,11)=ET_mean_mens(6,:);%juin
Y(:,12)=ET_mean_mens(7,:);%juillet

islist=[2 4 1 2 4 1 2 4 1 2 4 1];
varlist=[1 1 1 5 5 5 6 6 6 2 2 2];
varlist2=["Drainage" "D ET juin" "D ET juil"];


Y = Y(:,:); % Define the output metrics considered among those available in hymod_MulOut
n_out = size(Y,2); % Number of output metrics considered

flag = 2;  % Select the statistic to compute the sensitivity indices (here the maximum vertical distance or Kolmogorov-Smirnov statistic between CDFs) 
ngroup = 10; % Number of groups into which the output values are divided

for ii=1:n_out
    if ii<=12
     v=varlist(ii);
    xlabels(ii)=strcat(indic_stat(islist(ii)),var(v));  %légende
    else
    xlabels(ii)=varlist2(ii-12);
    end
    [Srsa(ii,:),idx(:,ii),Yk(:,ii)] = RSA_indices_groups(X,Y(:,ii),ngroup,flag); % Calculate the sensitivity indices for each output metric

end
xlabels_cell=cellstr(xlabels);
params_cell=cellstr(params);

% Plot the sensitivity indices:
figure_handle=figure()
 for ii = 1:n_out
    subplot(5,3,ii)
    title(xlabels(ii))
    boxplot1(Srsa(ii,:)) 

     % x-axis label:
    if ii >= 13 % last line
       set(gca,'XTickLabel',params_cell,'XGrid','On','FontSize',10,'FontWeight','bold')
      else
        set(gca,'XTickLabel',[],'XGrid','On')
    end     

    % y-axis label:
    aa=([1, 4, 7 ,10,13]);  %numéro subplots sur la 1ère colonne
    if ismember(ii,aa) % first column
        set(gca,'YtickLabel',[0 0.5 1],'YGrid','On','FontSize',10,'FontWeight','bold')
        ylabel('Sensitivity','FontWeight','bold')
    else
        set(gca,'YTickLabel',[],'YGrid','On')
    end
 end
 

    
%% Now we look for an explanation of the results obtained by plotting the input factors' CDFs

% Plot input factors' CDFs:
figure_handle=figure() % This figure is shown in Figure 2
n=12;
for ii = 1:n
   % nexttile
    RSA_plot_groups_b(X,idx(:,ii),Yk(:,ii),n,ii,params_cell); % Figure 

pos = get(subplot(n,5,(5*(ii-1))+5),'Position');
pos(1)=0.83;
axes('Position', pos, 'Visible', 'off');
for i=1:ngroup+1; ctick{i}=sprintf('%1.2g',Yk(i,ii)); end
colormap(feval('jet',ngroup))
c=colorbar('YTickLabel',ctick,'Fontsize',10);
ylabel(c,xlabels(ii),'Fontsize',10)
hold on
end

%%  visualisation   parallel coordinate plots

% quelq critères utiliser ? 
cas=1;

      
     
      
if cas==1    % tous les critères echantillonnés 
idx = result_seuils_kge(:,1:4);
criteres = cell(1, size(idx,2)); % Préallocation du cell array
for iii = 1:size(idx,2)
    dd = strcat("KGE ", var(iii)); % Utilisation de string pour éviter les cellules
    criteres{iii} = char(dd); % Conversion en chaîne de caractères si nécessaire
end
figure()
parcoor(X,params_cell,[],idx,criteres)
title("critères d'evaluation, les 20% des meilleures simulations")
end 


if cas==2  % différentes combinaisons de critères
figure()
parcoor(X,params_cell,[],idx)
end  

if cas ==3   % avec des différents pourcentages 
figure()
parcoor(X,params_cell)
end



%% seuils à définir
% les différents indicateurs stats et les variables 
% définir des seuils pour détereminer les simulations les plus performantes
% et ensuite réduire l'espace de paramètres 
is=2; %kge rmse corélation
j=4; %"LE", "RN","H","W5","W25","W50"
% LE 

% somme kge des variables %"LE", "H","W25","W50" et               ET jour
somme_kge=[];
for ii=1:isim
somme_kge(ii)=sum(stats.kge(ii,1:5))-stats.kge(ii,3);
end

%somme kge
figure()
scatter3(mat_param.matparam(:,13),mat_param.matparam(:,2),mat_param.matparam(:,1),[],somme_kge,'filled')
colorbar
xlabel(variables(13))
ylabel(variables(2))
zlabel(variables(1))
title('somme KGE')
set(gca, 'XScale', 'log')
colormap turbo



%% selectionner les X meilleurs simulations 
k=round(numsim/5);
%maxk(,k) 
result_seuils_kge=[];
result_seuils_rmse=[];
 for s=1:length(var)
     [as,idx_kge] = sort(stats.kge(:,s),'descend');
     result_seuils_kge(:,s) = idx_kge(1:k);
     [as,idx_rmse] = sort(stats.rmse(:,s),'ascend');
     result_seuils_rmse(:,s) = idx_rmse(1:k);
 end 
% rajouter le classement des meilleures sommes de kge
 [as,idx_kges] = sort(somme_kge,'descend');
 result_seuils_kge(:,8) = idx_kges(1:k);

% différences ET mensuel mai juin juillet
[as,idx_ETmai] = sort(abs(ET_mean_mens(5,:)),'ascend');
 result_seuils_ETmai(:) = idx_ETmai(1:k);
[as,idx_ETjuin] = sort(abs(ET_mean_mens(6,:)),'ascend');
 result_seuils_ETjuin(:) = idx_ETjuin(1:k);
 [as,idx_ETjul] = sort(abs(ET_mean_mens(7,:)),'ascend');
 result_seuils_ETjul(:) = idx_ETjul(1:k);

 % chaque critère ou une combinaison de plusieurs critères    liste variables%"LE", "RN","H","W5","W25","W50"
list=mintersect(result_seuils_kge(:,8),result_seuils_ETjuin,result_seuils_rmse(:,1));
%list=mintersect(result_seuils_rmse(:,3),result_seuils_rmse(:,1));


% ou sont situées les meilleures simulation dans l'espace de paramètres
figure()
%scatter3(mat_param.matparam(list,10),mat_param.matparam(list,13),mat_param.matparam(list,12),[],stat_temp(list,7),'filled')
scatter3(mat_param.matparam(:,1),mat_param.matparam(:,13),mat_param.matparam(:,2),'filled','r')
hold on
set(gca, 'YScale', 'log')
%scatter3(mat_param.matparam(list,10),mat_param.matparam(list,13),mat_param.matparam(list,12),[],stats.([string(indic_stat(is))])(list,j),'filled')
scatter3(mat_param.matparam(list,1),mat_param.matparam(list,13),mat_param.matparam(list,2),'filled','b')
%a=colorbar;
xlabel(variables(1))
ylabel(variables(13))
zlabel(variables(2))
a.Label.String = [indic_stat(is) var(j)];
%a.Label.String='somme kge LE RN H W5 W25 W50';

colormap turbo  %tracer dans ensemble des paramètres


%% seuils fixés manuellement
for ii=1:isim
    % LE kge 
    if stats.kge(ii,1) > 0.6  % LE kge 
    seuils(ii,1)=1;
    else
    seuils(ii,1)=0;
    end
    % LE RMSE
     if stats.rmse(ii,1) <35 % LE RMSE
    seuils(ii,2)=1;
    else
    seuils(ii,2)=0;
    end
 % RN RMSE
     if stats.rmse(ii,2) <21 % % RN RMSE
    seuils(ii,3)=1;
    else
    seuils(ii,3)=0;
     end
     % H RMSE
     if stats.rmse(ii,3) <58 % % H RMSE
    seuils(ii,4)=1;
    else
    seuils(ii,4)=0;
    end
% ET_jour RMSE
     if ET_jour.rmse(ii) <0.6 % % H RMSE
    seuils(ii,5)=1;
    else
    seuils(ii,5)=0;
    end
% ET_jour kge
     if ET_jour.rmse(ii) >0.5 % % H RMSE
    seuils(ii,6)=1;
    else
    seuils(ii,6)=0;
     end
 % ET_jour mars rmse
     if ET_mean_mens(3,ii) >-5 & ET_mean_mens(3,ii) <5% % H RMSE
    seuils(ii,7)=1;
    else
    seuils(ii,7)=0;
     end

  
location = knnsearch(ET_mean_mens(3,ii),0);

  % ET_jour avril rmse
     if ET_mean_mens(4,ii) >-5 & ET_mean_mens(4,ii) <5% % H RMSE
    seuils(ii,8)=1;
    else
    seuils(ii,8)=0;
     end

  % ET_jour mai rmse
     if ET_mean_mens(5,ii) >-10 & ET_mean_mens(5,ii) <10% % H RMSE
    seuils(ii,9)=1;
    else
    seuils(ii,9)=0;
     end

    % ET_jour juin rmse
     if ET_mean_mens(6,ii) >-10 
    seuils(ii,10)=1;
    else
    seuils(ii,10)=0;
     end

     % ET_jour juillet rmse
     if ET_mean_mens(7,ii) >-10 & ET_mean_mens(7,ii) <10% % H RMSE
    seuils(ii,11)=1;
    else
    seuils(ii,11)=0;
     end

  % ET_jour aout rmse
     if ET_mean_mens(8,ii) >-7 & ET_mean_mens(8,ii) <10% % H RMSE
    seuils(ii,12)=1;
    else
    seuils(ii,12)=0;
     end

  
  % ET_jour sep rmse
     if ET_mean_mens(9,ii) >-5 & ET_mean_mens(9,ii) <5% % H RMSE
    seuils(ii,13)=1;
    else
    seuils(ii,13)=0;
     end

  % Humidité 5cm kge 
    if stats.kge(ii,4) > 0.5  % LE kge 
    seuils(ii,14)=1;
    else
    seuils(ii,14)=0;
    end

  % Humidité 50cm kge 
    if stats.kge(ii,6) > 0.5  % LE kge 
    seuils(ii,15)=1;
    else
    seuils(ii,15)=0;
    end
end
seuils=array2table(seuils);
seuils.Properties.VariableNames = {'LE_kge' 'LE_rmse' 'RN RMSE' 'H RMSE' 'ET_j_RMSE' 'ET_j kge' 'ET_j_mars' 'ET_j_avril' 'ET_j_mai' 'ET_j_juin' 'ET_j_juil'  'ET_j_ug' 'ET_j_sep' 'hum_5cm_kge' 'hum_50cm_kge'};


%%

% liste des numeros de simulation qui respectent 
% chaque critère ou une combinaison de plusieurs critères 
list=mintersect(find(seuils.LE_kge),find(seuils.LE_rmse),find(seuils.ET_j_RMSE),find(seuils.ET_j_avril), find(seuils.ET_j_juin), find(seuils.ET_j_juil),  find(seuils.ET_j_sep))  ;

list=mintersect(find(seuils.LE_kge),find(seuils.LE_rmse),find(seuils.ET_j_RMSE),find(seuils.hum_5cm_kge),find(seuils.hum_5cm_kge));

% tracer dans ensemble des paramètres
figure()
scatter3(mat_param.matparam(list,10),mat_param.matparam(list,13),stats.([string(indic_stat(is))])(list,j),[],stats.([string(indic_stat(is))])(list,j),'filled')
colorbar
xlabel(variables(10))
ylabel(variables(13))
zlabel([indic_stat(is) var(j)])
set(gca, 'YScale', 'log')
colormap turbo

% somme kge 
figure()
scatter3(mat_param.matparam(:,10),mat_param.matparam(:,13),mat_param.matparam(:,12),[],somme_kge,'filled')
colorbar
xlabel(variables(10))
ylabel(variables(13))
zlabel(variables(12))
set(gca, 'YScale', 'log')
colormap turbo
title('somme kge LE RN H W5 W25 W50')


find(somme_kge==max(somme_kge))
seuils_min_max=[];

seuils_min_max(1)=min(mat_param.matparam(list,13));
seuils_min_max(2)=max(mat_param.matparam(list,13));
seuils_min_max(3)=min(mat_param.matparam(list,12));
seuils_min_max(4)=max(mat_param.matparam(list,12));
seuils_min_max(5)=min(mat_param.matparam(list,10));
seuils_min_max(6)=max(mat_param.matparam(list,10));
seuils_min_max=array2table(seuils_min_max);
seuils_min_max.Properties.VariableNames = {'min_list_restot' 'max_list_restot' 'min_list_resstom' 'max_list_resstom RMSE' 'min_list_dens_rac' 'max_list_dens_rac'};
