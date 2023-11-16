close all
clear all 
clc


%% Figure 6B ii
% EDC as a function of IL4

% Set the input parameters
NFkB=0; % nominal: 0, increases to 0.5

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

%% Get EDCss on PKC=2 while IL4=0-100

EDCVect= [];
ii=1;
PKC=2;

for IL4=0:10:1000
[ SS, number_SS ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);
%EDC_ss=(((a_EDC.*(aux_1/aux_2).*SS(:,2))./(1+i_EDC.*SS(:,1)+IL4))./d_EDC);
 EDCVect(ii)=SS(2,3);
ii=ii+1;
    end

%% Plot it

IL4Vect=0:10:1000;

p=polyfit(IL4Vect, EDCVect,2); %p =   0.0000   -0.0082    2.1447

y1 = polyval(p,IL4Vect);

figure(3)

plot (IL4Vect, EDCVect,'o','MarkerSize',5);
hold on
axis square
xlabel('IL4')
ylabel('EDC_s_s/PKC=2')
hold on
plot(IL4Vect,y1,'--','Color',[0,0,0], LineWidth=1)

