%%  COURBES DE VAN GENUCHTEN ET BROOKS & COREY 

% valeurs min et max 
h=[0.01 :0.01 :1000];
t=[0.1:0.01:1];
n1=2.2;
n2=2.8;
tetr1=0.05;
tetr2=0.2;
tets1=0.4;
tets2=0.7;
m1=1-(2/n1);
m2=1-(2/n2);
hg1=0.5;
hg2=5;
bet1=5;
bet2=15;
 ks1=5*10.^(-8);
 ks2=5*10.^(-6);

% VG
teta1=tetr1+(tets1-tetr1)./((1+(h/(hg1)).^n2).^m2);
teta2=tetr2+(tets1-tetr2)./((1+(h/(hg1)).^n2).^m2);
teta3=tetr2+(tets2-tetr2)./((1+(h/(hg1)).^n2).^m2);

figure()
semilogy(teta1,h,'LineWidth',1.5)
hold on 
semilogy(teta2,h,'LineWidth',1.5)
hold on 
semilogy(teta3,h,'LineWidth',1.5)
xlim([0.02,0.8]) 
ylim([0.09,1000])
legend('\theta r min \theta s min', '\theta r max \theta s min','\theta r max \theta s max')
xlabel('\theta (m3/m3)')
ylabel('|h| (m)')
grid on 
ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


% BC 

K1=ks2*((t./tets1).^bet1);
K2=ks2*((t./tets1).^bet1);
K3=ks2*((t./tets2).^bet1);

figure()
semilogy(t,K1,'LineWidth',2.5)
hold on 
semilogy(t,K2,'LineWidth',1.5)
hold on 
semilogy(t,K3,'LineWidth',1.5)
xlim([0.02,0.72]) 
%ylim([10^(-14),10^(-5)])
legend('\theta r min \theta s min', '\theta r max \theta s min','\theta r max \theta s max')
%legend('\beta min','\beta max')
xlabel('\theta (m3/m3)')
ylabel('K (m/s)')
grid on 
ax=gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 


% rescaling des teneurs en eau 
c= [0.85,0.5,0.4,0.35,0.3,0.2,0.15,0.1,0.1];  
sat=[0.11 0.080 0.060 0.054 0.078 0.072 0.072 0.072 0.072];
res=[0.025 0.025 0.016 0.014 0.024 0.019 0.019 0.019 0.019 ];
sat2=sat/0.2;
sat3=sat2.*c;
res2=res/0.2
res3=res2.*c
%%  Calculs de l'eau manquante sans l'infiltration préférentielle, à quel point les racines le rattrape ou pas   

surplus=ETR_day_sim2-ETR_day_obs;
dd=1330+60-365*3;
df=2060-365*4;
surplus=surplus(dd:df);
i_surplus_= find(surplus > 0.5);
 i_manque_= find(surplus < (-0.2));
 surplus_=surplus(i_surplus_);
 surplus_s=sum(surplus_)
manque_= surplus(i_manque_);
manque_s= sum(manque_)

figure()
    plot(days_sim2(dd:df),ETR_day_sim2(dd:df),'r-','linewidth',2)
   
  hold on 
    plot(days_obs(dd:df),ETR_day_obs(dd:df),'k-','linewidth',2)
    grid on
    ylabel('mm')
     
legend('Shallow (2m)','Deep (20m)','Karst (20m+cracks)')
    dynamicDateTicks([],[],'dd/mm')
    datetick('x','mmm-yy','keepticks')
       ax = gca;
ax.FontSize = 10; 
ax.FontWeight = 'bold'; 