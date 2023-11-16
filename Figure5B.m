close all
clear all 
clc


%% Figure 6A iv
% P+ as a function of NFkB

% Set the input parameters
IL4=0; % nomial: 0; increases to 100

v_Np63 =2;% 2; % 1st Positive feedback: effect of Np63 on Stat3 production
d_Stat3=1;%1; %Stat3 degradation rate    
va_Np63=10; %10, 1.5; % 2nd Positive feedback: max effect of Np63 on Np63 produciton
v_Stat3=1; %1, 3rd positive feedback: effect of Stat3 on Np63 produciton
% Together they give rise to alpha
k_Np63 =1.35; %1, AC50 fpr the 2nd positive feedback
n_H    =3;  %3, Hill coefficient for the 2nd positive feedback
d_Np63=6; %6, Nominal Np63 degradation rate
d_PKC=0.5;%1, .1; %PKC-mediated degradation rate of Np63

%For new Toufighi mean with the raw AMPs %FinalSep13
a_EDC=153.2607;
i_EDC=500.3689;
d_EDC=0.1029705;
aux_1= 477.7428;
aux_2=451.8896;
PKC=2.510846;

% the following two parameters are "part of the equations " for historical
% reasons (I.e. explorations of other model versions), but it turns out
% they are not needed, hence we set them to 0 always.
nu=0;
pIL4=0;
basal=0;
basalnp=0;

%% Search for P+ on PKC=0-3 while NFkB=0-1

PplusVect= [];
ii=1;
sensor=0;

for NFkB=0:0.01:0.5
for PKC=0:0.01:3
    [ SS, number_SS ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2); 
    if number_SS==2
        sensor=1;
    end
    if sensor==1 && number_SS==1
break 
    end
PplusVect(ii)=PKC;
end
ii=ii+1;
sensor=0;
end

%% Plot it

NFkBVect=0:0.01:0.5;

p=polyfit(NFkBVect, PplusVect,1); % p =   -1.7750    1.7550
y1 = polyval(p,NFkBVect);

figure(2)

plot (NFkBVect, PplusVect,'o','MarkerSize',8,'MarkerFaceColor','b');
hold on
axis square
xlabel('NFkB')
ylabel('P+')
%xlim([0,0.5])
hold on
plot(NFkBVect,y1,'--','Color',[0,0,0], LineWidth=1)
