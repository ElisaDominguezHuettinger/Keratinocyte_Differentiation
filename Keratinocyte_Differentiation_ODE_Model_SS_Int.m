function [ stable_positive_real_solution_matrix , counter_stable_positive_real_solution]=Keratinocyte_Differentiation_ODE_Model_SS_Int(PKC,NFkB, IL4, v_Np63,va_Np63, d_Stat3, k_Np63, n_H, v_Stat3, d_Np63, d_PKC, a_EDC, i_EDC, d_EDC, aux1, aux2)

    %% declare the model variables as symbolic variables
    syms Stat3_t Np63_t EDC_t aux_t
    
    
    %% write down the equations, equating them to 0, because we are looking for steady states
    dydt1= 0== PKC+v_Np63*Np63_t-d_Stat3*Stat3_t;  %dStat3_t_dt
    dydt2= 0== NFkB+va_Np63*(Np63_t^n_H/(k_Np63^n_H+Np63_t^n_H)) +v_Stat3*Stat3_t-d_Np63*Np63_t-PKC*d_PKC*Np63_t; %dNp63_t_dt
    dydt3= 0== a_EDC*aux_t/(1+i_EDC*Stat3_t+IL4)-d_EDC*EDC_t; %dEDC_t_dt 
    dydt4= 0== aux1*Np63_t-aux2*aux_t;
    
    
    %% put them into a vector
     equations = [dydt1 dydt2 dydt3 dydt4];
      % the same with the variables
    
      vars=[Stat3_t Np63_t EDC_t aux_t];
     %% ¡Solve the system! compute all the solutions in an open range :D

      range = [NaN NaN; NaN NaN;NaN NaN;NaN NaN];
     sol = vpasolve(equations, vars, range);
     
     
     
     %% now we will select only the stable positive real steady staes

     %% For this we need an expression for the Jacobian matrix
    dStat3= PKC+v_Np63*Np63_t-d_Stat3*Stat3_t;  %dStat3_t_dt
    dNp63 = NFkB+va_Np63*(Np63_t^n_H/(k_Np63^n_H+Np63_t^n_H)) +v_Stat3*Stat3_t-d_Np63*Np63_t-PKC*d_PKC*Np63_t; %dNp63_t_dt
    dEDC = a_EDC*aux_t/(1+i_EDC*Stat3_t+IL4)-d_EDC*EDC_t; %dEDC_t_dt 
    daux = aux1*Np63_t-aux2*aux_t;
    
    %% compute the jacobian matrix
    J=jacobian([dStat3 dNp63 dEDC daux], vars);

    %% declare a empty matrix where we will store the stable positive real solutions 
     stable_positive_real_solution_matrix=[];
     counter_stable_positive_real_solution=0;
     
    %% now we will loop over all the solutions and select only those that are positive, real and stable
         for sol_num=1:1:length(sol.Stat3_t)


             %% Check if solutions are real
             if isreal([sol.Stat3_t(sol_num) sol.Np63_t(sol_num) sol.EDC_t(sol_num) sol.aux_t(sol_num)])
                 %% check if solutions are positive (and real)
                 if sum(double(([sol.Stat3_t(sol_num) sol.Np63_t(sol_num) sol.EDC_t(sol_num) sol.aux_t(sol_num)])>=0))==4

                     %% Check if solutions are stable (and positive and real)
                     %% evaluate the jacobian matrix at the solution
                     Jeval=subs(J, vars, [sol.Stat3_t(sol_num) sol.Np63_t(sol_num) sol.EDC_t(sol_num) sol.aux_t(sol_num)]);

                     %% print these eigenvalues
                     eigenvals=eig(Jeval) ;

                     %% for a steady state to be stable, we need all the eigenvalues to be < 0
                     if (sum(double(eigenvals)<0))==4
                      %   disp('solution is stable :)')
                         %% we are collecting these solutions in a matrix
                         % matrix concatenates solutions
                         stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix;
                             double( [sol.Stat3_t(sol_num) sol.Np63_t(sol_num) sol.EDC_t(sol_num) sol.aux_t(sol_num)])];

                         counter_stable_positive_real_solution=counter_stable_positive_real_solution+1;


                     else
                     %    disp('solution is unstable :(')
                     end

                     %% print the solution number
                 %    [sol.Stat3_t(sol_num) sol.Np63_t(sol_num) sol.EDC_t(sol_num)]
                     
                 else
                   %  disp('solution has negative elements, choose another value')
                 end

             else
                % disp('solution is complex, we move to another value')
             end

         end

         %% total number of stable positive real solutions is:
        % counter_stable_positive_real_solution
         % = numero de renglones en matriz)


    
end