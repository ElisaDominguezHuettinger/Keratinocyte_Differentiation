function dydt=Keratinocyte_Differentiation_ODE_Model_Integral_Np63(~,y,PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux1, aux2)

    dydt = zeros(4,1);
    
    Stat3_t=y(1);
    Np63_t =y(2);
    EDC_t =y(3); % this is the output
    aux_t=y(4);

    dydt(1) = PKC+IL4+v_Np63*Np63_t-d_Stat3*Stat3_t;  %dStat3_t_dt
    dydt(2) = va_Np63*(Np63_t^n_H/(k_Np63^n_H+Np63_t^n_H)) +v_Stat3*Stat3_t-d_Np63*Np63_t-NFkB*Np63_t-PKC*d_PKC*Np63_t; %dNp63_t_dt
    dydt(3) = a_EDC*aux_t/(1+i_EDC*Stat3_t)-d_EDC*EDC_t; %dEDC_t_dt 
    dydt(4)=  aux1*Np63_t-aux2*aux_t;
    
    
end
