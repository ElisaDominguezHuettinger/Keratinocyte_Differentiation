close all
clear all 
clc


%% Figure 6A i-iii
% Set the input parameters

IL4=1000; % nomial: 0; increases to 100 %just move here, used 0,500 and 1000
NFkB=0; % nominal: 0, increases to 0.25. %just move here, used 0,.1 and 0.25
 
    if NFkB==0 && IL4==0
        colour='k';
    elseif NFkB>0 && IL4==0
        colour='b';
    elseif NFkB==0 && IL4>0
        colour='r';
    else
        colour='m';
    end

    markersize=15; %and here, used 5,15,35 and 65


%% parameters previously determined (fixed)
NFkB=0; %Effect: degrades Np63
IL4=0; %Effect: increases Stat3 producioton rate    
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



%% Create bifurcation diagram



figure(1)

subplot(1,3,1)
hold on
axis square
xlabel('PKC')
ylabel('Stat3_s_s')


subplot(1,3,2)
hold on
axis square
xlabel('PKC')
ylabel('Np63_s_s')


subplot(1,3,3)
hold on
axis square
xlabel('PKC')
ylabel('EDC_s_s')


%%
for PKC=0:0.01:3
    
    
    % get the stable steady states
    [ SS, number_SS ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);

    for ii=1:1:3
        subplot(1,3,ii)
 scatter(PKC.*ones(1,number_SS),SS(:,ii),markersize,  colour, 'c', 'filled')  
    end


    if number_SS>1 % if there is more than one SS
    
    US=Sepparatrix_Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);
 for ii=1:1:3
        subplot(1,3,ii)
    scatter(PKC,US(:,ii),markersize./10,   colour,  'o')    
 end


   end
    
   

end


