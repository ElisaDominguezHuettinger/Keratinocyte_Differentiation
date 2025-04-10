close all
clear all 
clc



%% %% Parameters - first estimate. With this parameter set the model displays bistability

% Initial Guess for parameters
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


%% First, we want to make sure that under these parameters the model displays bistability
  %For Toufighi mean with the raw AMPs %FinalSep13
a_EDC=153.2607;
i_EDC=500.3689;
d_EDC=0.1029705;
aux_1= 477.7428;
aux_2=451.8896;
PKC=2.510846;

  %For Borowiec FLG data
% a_EDC=465.91537;
% i_EDC=492.1363;
% d_EDC=0.0173;
% aux_1= 0.0570;
% aux_2=306.048;
% PKC=3.3147;

% the following parameters are "part of the equations " for historical
% reasons (I.e. explorations of other model versions), but it turns out
% they are not needed, hence we set them to 0 always.
nu=0;
pIL4=0;
basal=0;
basalnp=0;

%% compute the nullcline Stat3
 % syms PKC IL4 v_Np63 Np63_t d_Stat3 Stat3_t basal pIL4 nu
 % dStat3= PKC+basal+IL4+v_Np63*Np63_t-d_Stat3*Stat3_t;  %dStat3_t_dt
 % dStat3=
 % PKC+basal+pIL4*IL4+v_Np63*Np63_t-d_Stat3*Stat3_t-nu*Stat3_t*Np63_t; %this is the newer version but its equivalent to equation on line 45 bc pIL4 and nu = 0
 % solve(dStat3==0, Stat3_t)

%%
  % syms va_Np63 Np63_t n_H k_Np63 v_Stat3 Stat3_t d_Np63 PKC d_PKC basalnp NFkB
  % dNp63= basalnp+va_Np63*(Np63_t^n_H/(k_Np63^n_H+Np63_t^n_H)) +v_Stat3*Stat3_t-d_Np63*Np63_t-NFkB*Np63_t-PKC*d_PKC*Np63_t;
  % solve(dNp63==0, Stat3_t)

% ans   (Np63_t*d_Np63 - basalnp + NFkB*Np63_t - (Np63_t^n_H*va_Np63)/(Np63_t^n_H + k_Np63^n_H) + Np63_t*PKC*d_PKC)/v_Stat3;
%% declare functions for the nullclines
Stat_Nullcline=@(basal,IL4, PKC, Np63_t, v_Np63, d_Stat3)(IL4 + PKC + basal + Np63_t.*v_Np63)./d_Stat3;

% it is easy to see that this is a first order  of the form
%b=(PKC+basal)/dStat3;
%m=v_Np63/d_Stat3;


Np63_Nullcline=@(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB)(Np63_t.*d_Np63 - basalnp + NFkB.*Np63_t - (Np63_t.^n_H.*va_Np63)./(Np63_t.^n_H + k_Np63.^n_H) + Np63_t.*PKC.*d_PKC)./v_Stat3;


Np63_t=0:0.1:100;



%% Plot nullclines for different values of PKC
figure;

for PKC=0:2:4

    beta=(d_Np63+NFkB+PKC*d_PKC)/v_Stat3; % Note that it depends on NfKB and PKC

%figure
n1=plot(Np63_t,Stat_Nullcline(basal,IL4, PKC, Np63_t, v_Np63, d_Stat3), 'LineWidth', 2,'Color','b', 'DisplayName', 'Stat3 nullcline');
hold on
n2=plot(Np63_t,Np63_Nullcline(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB),'LineWidth', 2,'Color','m', 'DisplayName', 'Np63 nullcline');
hold on
%plot(Np63_t,(beta*Np63_t), 'LineWidth', 2,'Color','g', 'LineStyle','--');

[ stable_positive_real_solution_matrix , ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);


scatter(stable_positive_real_solution_matrix(:,2), stable_positive_real_solution_matrix(:,1), 100, 'k', 'filled')
end

ylim([0,9])
xlim([0, 2.5])

%legend([n1,n2],'AutoUpdate','off')

xlabel('Np63')
ylabel('Stat3')
axis square

title('Qualitative changes as PKC increases from 0 to 4')

%% Create bifurcation diagram
figure

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
xline(3.3147)
xline(2.510846)

% como en lo que tu usaste pero hay que modificar el codigo de 
for PKC=0:0.09:5
    
    
    % get the stable steady states
    [ SS, number_SS ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);

    for ii=1:1:3
        subplot(1,3,ii)
 scatter(PKC.*ones(1,number_SS),SS(:,ii),40,  'r')  
    end

    if number_SS>1 % if there is more than one SS
    
    US=Sepparatrix_Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);
 for ii=1:1:3
        subplot(1,3,ii)
    scatter(PKC,US(:,ii),40,  'b')    
 end
    end
    
    

end


%% Focus on PKC =2

figure

PKC=2;

alpha=va_Np63/v_Stat3;
beta=(d_Np63+NFkB+PKC*d_PKC)/v_Stat3; % Note that it depends on NfKB and PKC

n1=plot(Np63_t,Stat_Nullcline(basal,IL4, PKC, Np63_t, v_Np63, d_Stat3), 'LineWidth', 1+0.5*PKC,'Color','k', 'DisplayName', 'Stat3 nullcline');
hold on
n2=plot(Np63_t,Np63_Nullcline(va_Np63, Np63_t, n_H, k_Np63, v_Stat3, d_Np63, PKC, d_PKC, basalnp, NFkB),'LineWidth', 1+0.5*PKC,'Color','r', 'DisplayName', 'Np63 nullcline');

[ stable_positive_real_solution_matrix , ]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);

scatter(stable_positive_real_solution_matrix(:,2), stable_positive_real_solution_matrix(:,1), 150, 'c', 'filled')

ylim([0,10])
xlim([0, 4])

%legend([n1,n2],'AutoUpdate','off')

xlabel('Np63')
ylabel('Stat3')
axis square

% for the vector plot
[Np63_t_v, Stat3_t_v] = meshgrid(0:0.6:4,0:1:10);

dy1dt= basalnp+NFkB+va_Np63.*(Np63_t_v.^n_H./(k_Np63.^n_H+Np63_t_v.^n_H)) +v_Stat3.*Stat3_t_v-d_Np63.*Np63_t_v-PKC.*d_PKC.*Np63_t_v; %dNp63_t_dt
dx1dt= PKC+basal+IL4+v_Np63.*Np63_t_v-d_Stat3.*Stat3_t_v;  %dStat3_t_dt
   
quiver(Np63_t_v,Stat3_t_v,dy1dt,dx1dt,1, 'b')


%%
 unstable_sol=Sepparatrix_Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2);
%%
scatter(unstable_sol(1,2),unstable_sol(1,1), 100, 'm')
%%

% to characterize the sepparatrix, we will integrate back from the unstable
% solution


for epsilon=[-0.1, 0.1]
y0=[unstable_sol(1,1)+epsilon unstable_sol(1,2) 3 1] ;
%y0=[8.1061    6.2438    0.7462    6.2438]; %condiciones inicales optimizadas para fit final de Toufighi
tspan=[0:-0.001:-100];
SolOPT=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2),tspan,y0);

plot(SolOPT.y(2,:),SolOPT.y(1,:), ['--' 'g'],'LineWidth',2)

end