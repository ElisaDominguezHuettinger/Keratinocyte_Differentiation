close all
clear all 
clc



%% %% Parameters

NFkB=0; %Effect: degrades Np63
IL4=0; %Effect: increases Stat3 producioton rate
    
v_Np63 =2;% 2; % 1st Positive feedback: effect of Np63 on Stat3 production
d_Stat3=1;%1; %Stat3 degradation rate
    
va_Np63=10; %10, 1.5; % 2nd Positive feedback: max effect of Np63 on Np63 produciton
v_Stat3=1; %1, 3rd positive feedback: effect of Stat3 on Np63 produciton
% Together they give rise to alpha

k_Np63 =1; %1, AC50 fpr the 2nd positive feedback
n_H    =3;  %3, Hill coefficient for the 2nd positive feedback

d_Np63=6; %6, Nominal Np63 degradation rate
d_PKC=1;%1, .1; %PKC-mediated degradation rate of Np63


nu=0;
pIL4=0;
basal=0;
basalnp=0;


Np63_Nullcline=@(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB)(Np63_t.*d_Np63 - basalnp + NFkB.*Np63_t - (Np63_t.^n_H.*va_Np63)./(Np63_t.^n_H + k_Np63.^n_H) + Np63_t.*PKC.*d_PKC)./v_Stat3;


Np63_t=0:0.1:100;



%% Plot nullclines for different values of PKC
figure;

for PKC=0

    beta=(d_Np63+NFkB+PKC*d_PKC)/v_Stat3; % Note that it depends on NfKB and PKC
    alpha=va_Np63/v_Stat3;

plot(Np63_t,Np63_Nullcline(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB),'LineWidth', 2,'Color','m', 'DisplayName', 'Np63 nullcline');
hold on
plot(Np63_t,(beta*Np63_t), 'LineWidth', 2,'Color','m', 'LineStyle','--');
hold on
plot(Np63_t,((beta*Np63_t)-alpha), 'LineWidth', 2,'Color','m', 'LineStyle','--');

for PKC=4

    beta=(d_Np63+NFkB+PKC*d_PKC)/v_Stat3; % Note that it depends on NfKB and PKC
    alpha=va_Np63/v_Stat3;
hold on
plot(Np63_t,Np63_Nullcline(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB),'LineWidth', 2,'Color','c', 'DisplayName', 'Np63 nullcline');
hold on
plot(Np63_t,(beta*Np63_t), 'LineWidth', 2,'Color','c', 'LineStyle','--');
hold on
plot(Np63_t,((beta*Np63_t)-alpha), 'LineWidth', 2,'Color','c', 'LineStyle','--');
end
end 

ylim([0,10])
xlim([0, 3])

xlabel('Np63 (t)')
ylabel('Np63 nullcline')
axis square


