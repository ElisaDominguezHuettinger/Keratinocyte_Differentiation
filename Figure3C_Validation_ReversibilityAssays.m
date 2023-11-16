close all
clear all 
clc

%% %% OPTIMAL PARAMETERS
NFkB=0; %Effect: degrades Np63
IL4=0; %Effect: increases Stat3 producioton rate
    
v_Np63 =2;% 2; % 1st Positive feedback: effect of Np63 on Stat3 production
d_Stat3=1;%1; %Stat3 degradation rate
    
va_Np63=10; %1.5; % 2nd Positive feedback: max effect of Np63 on Np63 produciton
v_Stat3=1; % 3rd positive feedback: effect of Stat3 on Np63 produciton
% Together they give rise to alpha

k_Np63 =1.35; % AC50 fpr the 2nd positive feedback
n_H    =3;  % Hill coefficient for the 2nd positive feedback

d_Np63=6; % Nominal Np63 degradation rate
d_PKC=0.5; %PKC-mediated degradation rate of Np63

%For Toufighi final fit
a_EDC=153.2607;
i_EDC=500.3689;
d_EDC=0.1029705;
aux_1= 477.7428;
aux_2=451.8896;
PKC=2.510846;

%% Experimental data
CSTAraw=[12.6983333300000	13.0196666700000	13.4086666700000	13.5166666700000	13.6813333300000	13.8566666700000	13.8786666700000	13.9650000000000	13.8273333300000	13.8446666700000];
CXCL10raw=[6.30966666700000	8.72166666700000	7.25433333300000	6.49000000000000	6.43733333300000	6.88733333300000	6.76133333300000	6.50033333300000	6.42233333300000	6.86833333300000];
CXCL11raw=[6.46200000000000	8.17833333300000	7.74033333300000	7.05566666700000	6.94833333300000	7.41300000000000	7.19766666700000	6.89900000000000	6.97966666700000	6.99133333300000];
CSTA=[12.6983333300000	13.0196666700000	13.4086666700000	13.5166666700000	13.6813333300000	13.8566666700000	13.8786666700000	13.9650000000000	13.8273333300000	13.8446666700000]./max(CSTAraw);
CXCL10=[6.30966666700000	8.72166666700000	7.25433333300000	6.49000000000000	6.43733333300000	6.88733333300000	6.76133333300000	6.50033333300000	6.42233333300000	6.86833333300000]./max(CXCL10raw);
CXCL11=[6.46200000000000	8.17833333300000	7.74033333300000	7.05566666700000	6.94833333300000	7.41300000000000	7.19766666700000	6.89900000000000	6.97966666700000	6.99133333300000]./max(CXCL11raw);
RNASE7raw= [7.57800000000000	7.61066666700000	9.43300000000000	10.2956666700000	10.5896666700000	10.5196666700000	9.84100000000000	9.25266666700000	8.94600000000000	8.41533333300000];
S100A7raw= [7.43800000000000	8.10833333300000	8.97933333300000	9.75766666700000	11.0156666700000	11.9436666700000	12.1403333300000	12.2493333300000	12.2816666700000	12.2896666700000];
SLPIraw= [9.96866666700000	11.3613333300000	12.5953333300000	13.5186666700000	14.3783333300000	14.8773333300000	14.8920000000000	14.8293333300000	14.8000000000000	14.7493333300000];
RNASE7= [7.57800000000000  7.61066666700000	9.43300000000000	10.2956666700000	10.5896666700000	10.5196666700000	9.84100000000000	9.25266666700000	8.94600000000000	8.41533333300000]./max(RNASE7raw);
S100A7= [7.43800000000000  8.10833333300000	8.97933333300000	9.75766666700000	11.0156666700000	11.9436666700000	12.1403333300000	12.2493333300000	12.2816666700000	12.2896666700000]./max(S100A7raw);
SLPI= [9.96866666700000  11.3613333300000	12.5953333300000	13.5186666700000	14.3783333300000	14.8773333300000	14.8920000000000	14.8293333300000	14.8000000000000	14.749333330000]./max(SLPIraw);
FLG=[0.8100    0.8009    0.8626    0.9293    0.9534    0.9794    1.0000    0.9847    0.9553    0.9400];
IVL=[0.7462    0.7756    0.8323    0.8751    0.9141    0.9632    0.9754    0.9930    0.9969    1.0000];

M=[CSTA;CXCL10;CXCL11;RNASE7;S100A7;SLPI;FLG;IVL];

se = std(M) / sqrt(size(M,1));

EDC_exp=mean([CSTA;CXCL10;CXCL11;RNASE7;S100A7;SLPI;FLG;IVL],1);
t_exp =[    0     5    10    15    20    25    30    35    40    45];

%% Plot exp data
figure
hold on
plot(t_exp, CSTA,'o' ,'MarkerEdgeColor','[0.4471    0.7176    0.9882]','MarkerSize',15,'LineWidth',2);
plot(t_exp, CXCL10,'d' ,'MarkerEdgeColor','[0    0.8157    0.9412]','MarkerSize',15,'LineWidth',2);
plot(t_exp, CXCL11,'o' ,'MarkerEdgeColor','[ 0.1216    0.4902    0.2196]','MarkerSize',15,'LineWidth',2);
plot(t_exp, RNASE7,'d' ,'MarkerEdgeColor','[0.5922    1.0000    0.3216]','MarkerSize',15,'LineWidth',2);
plot(t_exp, S100A7,'d' ,'MarkerEdgeColor','[1.0000    0.8824    0.1216]','MarkerSize',15,'LineWidth',2);
plot(t_exp, SLPI,'d' ,'MarkerEdgeColor','[1.0000    0.5255    0.1098]','MarkerSize',15,'LineWidth',2);
plot(t_exp, FLG,'o' ,'MarkerEdgeColor','[0.4510    0.2078    0.9412]','MarkerSize',15,'LineWidth',2);
plot(t_exp, IVL,'o' ,'MarkerEdgeColor','[0.6706         0         0]','MarkerSize',15,'LineWidth',2);
plot(t_exp, EDC_exp,['pentagram' 'k'],'MarkerSize',19,'LineWidth',2,'MarkerFaceColor','k');


%% Integrate the model
t_exp =[    0     5    10    15    20    25    30    35    40    45];

EDC_exp=mean([CSTA;CXCL10;CXCL11;RNASE7;S100A7;SLPI;FLG;IVL],1);
tspan=[0:0.1:150];
y0= [8.1061    6.2438    0.7462    6.2438]; %Toufighi final fit
SolOPT=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2),tspan,y0);

figure
hold on
% plot(t_exp, RNASE7,['s' 'k'],'MarkerSize',10,'LineWidth',2);
% plot(t_exp, S100A7,['s' 'k'],'MarkerSize',10,'LineWidth',2);
% plot(t_exp, SLPI,['s' 'k'],'MarkerSize',10,'LineWidth',2);
% plot(t_exp, FLG,['o' 'r'],'MarkerSize',10,'LineWidth',2);
% plot(t_exp, IVL,['d' 'b'],'MarkerSize',10,'LineWidth',2);
errorbar(t_exp,EDC_exp,se,"both","o","LineWidth",1.5, Color='k')
hold on
plot(SolOPT.x,SolOPT.y(3,:),['--' 'm'],'LineWidth',2)
yline(0.9354,'LineWidth',2) %% for toufighi Final

hold on
axis square
xlabel('time (Days post-Ca challenge')
ylabel('EDC, mean over individual "dynamic"')
xlim([0 50])

%% Check SS
[stable_positive_real_solution_matrix , counter_stable_positive_real_solution]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2)

%% reversibility experiments
%(Jadali and Ghazizadeh, 2010) Fig 1, expression levels (PCR), Webplot Digitizer, Mouse keratinocytes

time_R=[	0	72	120	168	240]; %d
Ca_R=[	0.05	1.2	0.05	0.05	0.05]; %mM
p63_R=[	105.35	42.45	24.67	42.73	153.33]; %AU
INV_R=[	4.14	104.41	135.54	88.24	91.12]; %AU
FIL_R=[	10.88	102.71	120.98	70.32	22.18];
LOR_R=[	7.71	107.13	95.25	45.01	23.23];

   
mean_EDC=mean([INV_R; LOR_R; FIL_R]);
mean_EDC=mean_EDC./max(mean_EDC);
p63_R=p63_R./(max(p63_R));

figure(2);
subplot(3,2,1)
scatter(time_R, Ca_R, 'ok', 'filled')
hold on
line([time_R(1) time_R(2)], [Ca_R(2) Ca_R(2)],'Color', 'k')
line([time_R(1) time_R(1)], [Ca_R(1) Ca_R(2)],'Color', 'k')
line([time_R(2) time_R(2)], [Ca_R(2) Ca_R(1)],'Color', 'k')
line([time_R(2) time_R(end)], [Ca_R(1) Ca_R(1)],'Color', 'k')
axis square
ylabel('input (Ca mM)')
xlabel('time (days)')


subplot(3,2,3)
scatter(time_R, p63_R, 'rd', 'filled')
hold on
line([time_R(1) time_R(2)], [1 1],'Color', 'k')
line([time_R(1) time_R(1)], [Ca_R(1) 1],'Color', 'k')
line([time_R(2) time_R(2)], [1 Ca_R(1)],'Color', 'k')
line([time_R(2) time_R(end)], [Ca_R(1) Ca_R(1)],'Color', 'k')
axis square
ylabel('p63 (AU)')
xlabel('time (days)')


subplot(3,2,5)
scatter(time_R, INV_R./max(INV_R), 'bs', 'filled')
hold on
scatter(time_R, FIL_R/max(FIL_R), 'bd', 'filled')
scatter(time_R, LOR_R/max(LOR_R), 'bo', 'filled')
scatter(time_R, mean_EDC, 'bp', 'filled')
line([time_R(1) time_R(2)], [1 1],'Color', 'r')
line([time_R(1) time_R(1)], [Ca_R(1) 1],'Color', 'r')
line([time_R(2) time_R(2)], [1 Ca_R(1)],'Color', 'r')
line([time_R(2) time_R(end)], [Ca_R(1) Ca_R(1)],'Color', 'r')
axis square
legend('INV', 'LOC', 'FLG', 'mean')
ylabel('differentiation markers (AU)')
xlabel('time (days)')


sgtitle('Jadali and Ghazizadeh, 2010, reversibility experiments, Fig 1')

%(Jadali and Ghazizadeh, 2010) Fig 1, expression levels (PCR) from Mouse keratinocytes, measured by Webplot Digitizer.

% Consistent with previous data (eg. Toufigli, lena), upon addition of Ca, p63 decreases and 
% the differentiation markers increase.
% 
% Then, after removal of Ca, the epidermal differntiation 
% markers decrease again.
% Interestingly, INV and LOC keep on oncreasing for two days
% after decreasing Calcium levels before decreasing to lower (INV) and
% basal (LOC) levels.
% 
% This experimental dataset suggests that epidermal differentiation is reversible
% 
% Note that,  under low Ca conditions, EDC decrease while p63, its activator, increasaes, possibly because
% Stat3 is increasing. (this is very speculative).


%% Declare the initial conditions

%y0= [7.9794    6.2438    0.7694    6.2438]; %optimized
y0= [0.9794    0.2438    mean_EDC(1)   6.2438];%experimental

tspanREV1= [0 72];

PKC=3.0145;
% now we integrate the model
SolREV1=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2),tspanREV1,y0);

% now we widraw the Calcium and simulate again:
PKC=0.01;

y0_int_2=SolREV1.y(:,end);
%%
tspanREV2= [0 168];


SolREV2=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2),tspanREV2,y0_int_2);
%%

timeREV=[SolREV1.x SolREV2.x+72];
Stat3REV=[SolREV1.y(1,:) SolREV2.y(1,:)];
p63REV=[SolREV1.y(2,:) SolREV2.y(2,:)];

EDCREV=[SolREV1.y(3,:) SolREV2.y(3,:)];

%% plot besides the reversibility experimental data

subplot(3,2,2)
plot(timeREV, Stat3REV,'LineWidth', 2, 'Color', 'g')
hold on

line([time_R(1) time_R(2)], [0.9*max(Stat3REV) 0.9*max(Stat3REV)],'Color', 'k')
line([time_R(1) time_R(1)], [0.1*max(Stat3REV) 0.9*max(Stat3REV)],'Color', 'k')
line([time_R(2) time_R(2)], [0.9*max(Stat3REV) 0.1*max(Stat3REV)],'Color', 'k')
line([time_R(2) time_R(end)], [0.1*max(Stat3REV) 0.1*max(Stat3REV)],'Color', 'k')

axis square
ylabel('Stat3 (AU)')
xlabel('time (days)')


subplot(3,2,4)
plot(timeREV, p63REV, 'LineWidth', 2, 'Color', 'r')
hold on

line([time_R(1) time_R(2)], [0.9*max(p63REV) 0.9*max(p63REV) ],'Color', 'k')
line([time_R(1) time_R(1)], [0.1*max(p63REV)  0.9*max(p63REV) ],'Color', 'k')
line([time_R(2) time_R(2)], [0.9*max(p63REV)  0.1*max(p63REV) ],'Color', 'k')
line([time_R(2) time_R(end)], [0.1*max(p63REV)  0.1*max(p63REV) ],'Color', 'k')

axis square
ylabel('p63 (AU)')
xlabel('time (days)')


subplot(3,2,6)
plot(timeREV, EDCREV, 'LineWidth', 2, 'Color', 'b ')
hold on
scatter(time_R, INV_R./max(INV_R), 'bs', 'filled')
scatter(time_R, FIL_R/max(FIL_R), 'bd', 'filled')
scatter(time_R, LOR_R/max(LOR_R), 'bo', 'filled')
scatter(time_R, mean_EDC, 'bp', 'filled')

line([time_R(1) time_R(2)], [0.9*max(EDCREV)  0.9*max(EDCREV)],'Color', 'r')
line([time_R(1) time_R(1)], [0.1*max(EDCREV) 0.9*max(EDCREV)],'Color', 'r')
line([time_R(2) time_R(2)], [0.9*max(EDCREV) 0.1*max(EDCREV)],'Color', 'r')
line([time_R(2) time_R(end)], [ 0.1*max(EDCREV)  0.1*max(EDCREV)],'Color', 'r')

axis square
ylabel('differentiation markers (AU)')
xlabel('time (days)')

%% Plot only EDC

figure
plot(timeREV, EDCREV, 'LineWidth', 2, 'Color', 'b ')
hold on
scatter(time_R, INV_R./max(INV_R), 'bs', 'filled')
scatter(time_R, FIL_R/max(FIL_R), 'bd', 'filled')
scatter(time_R, LOR_R/max(LOR_R), 'bo', 'filled')
scatter(time_R, mean_EDC, 'bp', 'filled')

line([time_R(1) time_R(2)], [0.9*max(EDCREV)  0.9*max(EDCREV)],'Color', 'r')
line([time_R(1) time_R(1)], [0.1*max(EDCREV) 0.9*max(EDCREV)],'Color', 'r')
line([time_R(2) time_R(2)], [0.9*max(EDCREV) 0.1*max(EDCREV)],'Color', 'r')
line([time_R(2) time_R(end)], [ 0.1*max(EDCREV)  0.1*max(EDCREV)],'Color', 'r')

axis square
ylabel('differentiation markers (AU)')
xlabel('time (days)')
xlim([0 240])

