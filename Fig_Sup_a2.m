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
%Toufighi, K. et al. (2015) ‘Dissecting the Calcium-Induced Differentiation of Human Primary Keratinocytes Stem Cells by Integrative and Structural Network Analyses’, PLoS Computational Biology, 11(5), pp. 1–27. doi: 10.1371/journal.pcbi.1004256.

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
t_exp =[    0     5    10    15    20    25    30    35    40    45]; %days

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
xlabel('Time post Ca challenge (days)')
ylabel('Terminal differentation markers')
