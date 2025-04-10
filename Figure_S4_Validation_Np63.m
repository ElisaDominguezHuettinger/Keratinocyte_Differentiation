close all
clear all 
clc


%% Input the experimental data - it was taken from this publication: 
%Toufighi, K. et al. (2015) ‘Dissecting the Calcium-Induced Differentiation of Human Primary Keratinocytes Stem Cells by Integrative and Structural Network Analyses’, PLoS Computational Biology, 11(5), pp. 1–27. doi: 10.1371/journal.pcbi.1004256.
Toufighi_t =[    0     5    10    15    20    25    30    35    40    45];%./24; %time scale: Days

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


% the Terminal differentiation markers corresond the mean value of the
% epidermal differentiation complex (EDC) genes above.
EDC_exp=mean([CSTA;CXCL10;CXCL11;RNASE7;S100A7;SLPI;FLG;IVL],1);

% normalized to the fist timepoint ( the max)
Np63_exp=[7.721333333	7.14	6.675	6.946333333	7.060333333	6.731666667	6.723666667	6.757	6.810666667	6.883333333]./7.721333333;
SE_Np63_exp=[0.16264105	0.083032122	0.045177428	0.046997636	0.297690182	0.052650842	0.064957251	0.062268237	0.065547269	0.043811465];

% normalized to the max
Stat3_exp=[10.32066667	10.96533333	11.15966667	11.16066667	11.167	11.006	10.75666667	10.59766667	10.35066667	10.505]./11.1670;
SE_Stat3_exp=[0.070638359	0.063125096	0.152096826	0.026206445	0.15981969	0.235867194	0.267181794	0.110354479	0.256115425	0.161599299];

%% Borowiec et al., 2013 FLG
BorowiecFLG_t=[0  24  48  72  96  124];
BorowiecFLG_EDC=[0.605700685 1  0.656454003  0.362359611  0.300950191  0.249947883];

%% declare the data we will use for the optimization
data=EDC_exp;
t_exp=Toufighi_t;

 % data=BorowiecFLG_EDC;
 % t_exp=BorowiecFLG_t;

%% Integration time
tspan=[0:0.1:150];

%%  Parameters -  the ones that were optimized
%Final Stat3 and Np63 parameter values

NFkB=0; %Effect: degrades Np63
IL4=0; %Effect: increases Stat3 producioton rate   
v_Np63 =2; %1st Positive feedback: effect of Np63 on Stat3 production
d_Stat3=1; %Stat3 degradation rate   
va_Np63=10; %2nd Positive feedback: max effect of Np63 on Np63 produciton
v_Stat3=1; %3rd positive feedback: effect of Stat3 on Np63 produciton
k_Np63 =1.35; %AC50 fpr the 2nd positive feedback
n_H=3;  %Hill coefficient for the 2nd positive feedback
d_Np63=6; %Nominal Np63 degradation rate
d_PKC=0.5; %PKC-mediated degradation rate of Np63

%EDC values for Toufighi fit

a_EDC=153.2607;
i_EDC=500.3689;
d_EDC=0.1029705;
aux_1= 477.7428;
aux_2=451.8896;
PKC=2.510846;

nu=.045;


y0_int=[8.106053    6.2438    0.7462    6.2438];
y0=[8.106053    6.2438    0.7462    6.2438];
%% check how the model works with an initial guess
SolOPT=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux_1, aux_2),tspan,y0_int);

figure
plot(Toufighi_t, EDC_exp,['o' 'k'],'MarkerSize',10,'LineWidth',2);
hold on
%plot(BorowiecFLG_t, BorowiecFLG_EDC,['s' 'b'],'MarkerSize',10,'LineWidth',2);
plot(SolOPT.x,SolOPT.y(3,:),['--' 'm'],'LineWidth',2)
%yline(0.0975,'LineWidth',2) %% this is the EDC_ss for Borowiec Final

 
hold on
axis square
xlabel('time (Days post-Ca challenge')
ylabel('EDC, mean over individual "dynamic"')
xlim([0 50])


x0=[a_EDC, i_EDC, d_EDC, aux_1, aux_2, PKC, y0_int(1)];
Cost0=CostFunction_Keratinocytes_INT(x0,t_exp,data, y0, tspan, PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63,d_PKC, aux_1, aux_2)
title([' Initial guess, cost ' num2str(Cost0) 'x=' num2str(x0) ' PKC= ' num2str(PKC) ])


%% Stat3 - for this the experimental data is "unresolved", ie. we no not see any real trend, therefore we will skip it
figure;
subplot(2,1,1)
plot(Toufighi_t, Stat3_exp./max(Stat3_exp),['o' 'k'],'MarkerSize',10,'LineWidth',2);
hold on
axis square
xlabel('time (days post-Ca challenge (toufighi data)')
ylabel('Stat3')
plot(Toufighi_t, Stat3_exp./max(Stat3_exp)+SE_Stat3_exp, 'Color', 'k', 'LineStyle', '--', 'LineWidth',2);
plot(Toufighi_t, Stat3_exp./max(Stat3_exp)-SE_Stat3_exp, 'Color', 'k', 'LineStyle', '--', 'LineWidth',2);

subplot(2,1,2)
plot(SolOPT.x,SolOPT.y(1,:),['--' 'm'],'LineWidth',2)
axis square
xlabel('time (days post-Ca challenge')
ylabel('Stat3')
title(' "unresolved" ')
xlim([0,  45])

%% Np63
figure
subplot(1,3,2)
plot(Toufighi_t, Np63_exp,['o' 'r'],'MarkerSize',10,'LineWidth',2);
hold on
plot(Toufighi_t, Np63_exp+SE_Np63_exp, 'Color', 'r', 'LineStyle', '--', 'LineWidth',2);
plot(Toufighi_t, Np63_exp-SE_Np63_exp, 'Color', 'r', 'LineStyle', '--', 'LineWidth',2);
xlabel('time (days post-Ca challenge')
ylabel('Np63 (Toufighi et al. 2015)')
xlim([0,  45])
%ylim([0, 1])
axis square
ax = gca;
ax.FontSize = 16;

subplot(1,3,1)
plot(SolOPT.x,SolOPT.y(2,:),['--' 'm'],'LineWidth',2)
%plot(Sol.x/24,Sol.y(2,:)./Sol.y(2,1),['--' 'm'],'LineWidth',2)
%plot(Sol2.x,Sol2.y(2,:),['--' 'm'],'LineWidth',2)
xlabel('time (days post-Ca challenge')
ax = gca;
ax.FontSize = 16; 
ylabel('Np63 (Model)')
xlim([0,  45])
axis square

%% Add the lena data, simply plot it


%(Lena et al., 2008) Fig 1a bottom panel
Lena_tExp=[0	24	72	120];
Lena_Np63=[0.702965272	0.540229576	0.263155137	0.094690465];


subplot(1,3,3)
plot(Lena_tExp, Lena_Np63,['s' 'r'],'MarkerSize',10,'LineWidth',2);

xlabel('time (days post-Ca challenge')
ylabel('Np63 ((Lena et al., 2008) Fig 1a bottom panel')
xlim([0,  120])
%ylim([0, 1])
axis square
ax = gca;
ax.FontSize = 16;


