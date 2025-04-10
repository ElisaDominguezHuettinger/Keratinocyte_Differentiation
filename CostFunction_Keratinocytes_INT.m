function Cost=CostFunction_Keratinocytes_Sep2023_INT(x,t_exp,data, y0, tspan, PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63,d_PKC, aux_1, aux_2)

%Parameters of EDC

a_EDC=x(1);
i_EDC=x(2);
d_EDC=x(3);
aux_1=x(4);
aux_2=x(5);

PKC=x(6);

%if length(x)==5
%    y0_opt=y0;
%else
y01_opt=x(7);
% y02_opt=x(8);
% y03_opt=x(9);
% y04_opt=x(10);
% y0_opt=[y01_opt y02_opt y03_opt y04_opt];
y0_opt=y0;
y0_opt(1)=y01_opt;
%end


% We integrate the model
Sol=ode15s(@(t,y)Keratinocyte_Differentiation_ODE_Model_Integral_Np63(t,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC,aux_1, aux_2),tspan,y0_opt);
    
%%
t1=Sol.x;
%Stat3_t=Sol.y(1,:);
%Np63_t =Sol.y(2,:);
EDC_t  =Sol.y(3,:); 
    
% Predictions

% interpolate those values corresponding to the experiments
EDC_pre_Interpol= interp1(t1,EDC_t,t_exp);

% Calculate the cost of the predition vs. the experimental data
Cost=(sum((EDC_pre_Interpol-data).^2))./length(data);

% Initial conditions - in the absence of PKC they should correspond to the
% experimental data


end