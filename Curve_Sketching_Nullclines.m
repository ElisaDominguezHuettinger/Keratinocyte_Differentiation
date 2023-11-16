close all
clear all
clc

%%
% symbolic analysis of nullclines


syms IL4 PKC Np63_t v_Np63 d_Stat3 nu
syms alpha beta n_H k_Np63

Dispay_symbolic_solutions=0;
%%

if Dispay_symbolic_solutions==1
%% obtain nullclines

syms Stat3_t
dStat3_t_dt=PKC+IL4+v_Np63*Np63_t-d_Stat3*Stat3_t-nu*Np63_t*Stat3_t;  %dStat3_t_dt
solve(dStat3_t_dt==0, Stat3_t)

%%
syms va_Np63 v_Stat3 d_Np63 NFkB PKC d_PKC
dNp63_t_dt=va_Np63*(Np63_t^n_H/(k_Np63^n_H+Np63_t^n_H)) +v_Stat3*Stat3_t-d_Np63*Np63_t-NFkB*Np63_t-PKC*d_PKC*Np63_t; %dNp63_t_dt

solve(dNp63_t_dt==0, Stat3_t)
pretty((Np63_t*d_Np63 + NFkB*Np63_t - (Np63_t^n_H*va_Np63)/(Np63_t^n_H + k_Np63^n_H) + Np63_t*PKC*d_PKC)/v_Stat3)

%% to see that this is equivalent to the expression with alpha and beta, we collect terms
collect((Np63_t*d_Np63 + NFkB*Np63_t - (Np63_t^n_H*va_Np63)/(Np63_t^n_H + k_Np63^n_H) + Np63_t*PKC*d_PKC)/v_Stat3,Np63_t)

pretty(((NFkB + d_Np63 + PKC*d_PKC)/v_Stat3)*Np63_t - (Np63_t^n_H*va_Np63)/(v_Stat3*(Np63_t^n_H + k_Np63^n_H)))
    
%% declare functions for the nullclines
% Stat_Nullcline=(IL4 + PKC + Np63_t.*v_Np63)./(d_Stat3+nu.*Np63_t);

%Np63_Nullcline=-alpha.*(Np63_t.^n_H./(k_Np63.^n_H+Np63_t.^n_H))+Np63_t.*beta;

Np63_Nullcline(Np63_t)=-alpha.*(Np63_t.^n_H./(k_Np63.^n_H+Np63_t.^n_H))+Np63_t.*beta
diff(Np63_Nullcline,Np63_t)
solve(diff(Np63_Nullcline,Np63_t)==0, Np63_t)

diff(diff(Np63_Nullcline,Np63_t), Np63_t)
solve(diff(diff(Np63_Nullcline,Np63_t))==0, Np63_t)

subs(diff(diff(Np63_Nullcline,Np63_t), Np63_t), Np63_t,k_Np63 )
%% kurvendiscussion
subs(Stat_Nullcline,Np63_t,0)
%%
limit(Stat_Nullcline,Np63_t,Inf)
%%
solve(diff(Stat_Nullcline, Np63_t)==0, Np63_t) % empty set; this is to show that the nullcine is strictly monotone
end


%% Declare nullclines as functions

Np63_t_Nullcline=@(beta,alpha,n_H,k_Np63,Np63_t)beta.*Np63_t-alpha.*(Np63_t.^n_H./(k_Np63.^n_H+Np63_t.^n_H));
Stat_Nullcline=@(IL4, PKC, Np63_t, v_Np63, d_Stat3, nu)(IL4 + PKC + Np63_t.*v_Np63)./(d_Stat3+nu.*Np63_t);

%% we create the Np63 vector
Np63_t=0:0.01:15; 
%% First the Stat3 nullcline

IL4=0;
PKC=0;
v_Np63=1;
d_Stat3=1;
nu=1;

%%
(IL4 + PKC)/d_Stat3 % interseccion en 0 subs(Stat_Nullcline,Np63_t,0)
v_Np63/nu % limit(Stat_Nullcline,Np63_t,Inf)



%%
figure
plot(Np63_t, Stat_Nullcline(IL4, PKC, Np63_t, v_Np63, d_Stat3, nu), 'Color','k','LineWidth',3)
xlabel('Np63(t)');
ylabel('Stat3 Nullcine');
hold on
axis square

line([0, max(Np63_t)], [v_Np63/nu, v_Np63/nu], 'Color','c','LineStyle', ':','LineWidth',2)
text(max(Np63_t), v_Np63/nu, ' v_N_p_6_3/\nu', 'Color','c');%, 'FontSize',20)

scatter(0, (IL4 + PKC)/d_Stat3, 100, 'm','Filled')
text(0, (IL4 + PKC)/d_Stat3, ' (IL4 + PKC)/d_S_t_a_t_3', 'Color','m');%, 'FontSize',20)

fontsize(gcf,scale=2)

%% Next, we analyze the Np63 nullcline

% we lump parameters together: 

%alpha=va_Np63/v_Stat3;
%beta=(d_Np63+NFkB+PKC*d_PKC)/v_Stat3; % Note that it depends on NfKB and PKC

alpha=1;
n_H=10;
k_Np63=5;
%%
beta=0;

%%
figure;
%subplot(1,2,1)
hold on
plot(Np63_t, Np63_t_Nullcline(beta, alpha,n_H,k_Np63,Np63_t), 'Color','k','LineWidth',2)
xlabel('Np63_t');
ylabel('Np63(t) nullcline with \beta=0');
hold on
axis square

line([k_Np63, k_Np63], [-1,0], 'Color','m','LineStyle', '--','LineWidth',2)
line([0, max(Np63_t)], [-alpha,-alpha],'Color','r','LineStyle', ':','LineWidth',2) %linea arriba
line([0, max(Np63_t)], [-.5*alpha,-.5*alpha],'Color','y','LineStyle', '--','LineWidth',2) %linea arriba

text(k_Np63, -.4*alpha, ' k_N_p_6_3', 'Color','m');%, 'FontSize',20)
text(max(Np63_t),-alpha, ' -\alpha', 'Color','r');%, 'FontSize',20)
fontsize(gcf,scale=2)
%%
%subplot(1,2,2)
figure;

for  beta=[0, 0.1, .5, 1, 2]

hold on
plot(Np63_t, Np63_t_Nullcline(beta, alpha,n_H,k_Np63,Np63_t), 'Color','k','LineWidth',3*beta+1)
xlabel('Np63(t)');
ylabel('Np63(t) nullcline');
hold on
axis square
ylim([-alpha, 10])
line([0, max(Np63_t)], [0-alpha, beta*max(Np63_t)-alpha], 'Color','g','LineStyle', ':','LineWidth',2)
line([0, max(Np63_t)], [0, 0], 'Color','y','LineStyle', '--','LineWidth',2)
line([0, max(Np63_t)], [0, beta*max(Np63_t)], 'Color','b','LineStyle', ':','LineWidth',2)

if beta==.5
text(0.5*max(Np63_t), beta*0.5*max(Np63_t)-alpha, ' \beta \times Np63(t)-\alpha', 'Color','g');%, 'FontSize',20)
end
%fontsize(gcf,scale=2)
end
line([k_Np63, k_Np63], [-1,10], 'Color','m','LineStyle', '--','LineWidth',2)

text(k_Np63, 10, ' k_N_p_6_3', 'Color','m');%, 'FontSize',20)

%sgtitle('Np63 nullcline')
fontsize(gcf,scale=2)
%%
%beta=5
%plot(Np63_t, Np63_t_Nullcline(beta, alpha,n_H,k_Np63,Np63_t), 'Color','k','LineWidth',3)


%% difference Np63 nullcline and beta*max(Np63_t)-alpha

%% we want to see where the inflexion point is; here we show that

%{
syms IL4 PKC Np63_t v_Np63 d_Stat3 nu

syms alpha beta n_H k_Np63

Np63_Nullcline=-alpha*(Np63_t^n_H/(k_Np63^n_H+Np63_t.^n_H))+Np63_t*beta;

asymptote=beta*Np63_t-alpha

%subs(Np63_Nullcline, Np63_t, 0)-asymptote
%%
simplify(subs(Np63_Nullcline-asymptote, Np63_t, 0)) %alpha
%%
limit(Np63_Nullcline-asymptote,Np63_t,Inf) % the difference is 0

%%
subs(Np63_Nullcline-asymptote, Np63_t, k_Np63)
%}

%%
beta=.1
figure;

n_H=300

hold on
plot(Np63_t, Np63_t_Nullcline(beta, alpha,n_H,k_Np63,Np63_t), 'Color','k','LineWidth',3*beta+1)
xlabel('Np63(t)');
ylabel('Np63(t) nullcline');
hold on
axis square
ylim([-alpha,1])
xlim([4, 6])

line([0, max(Np63_t)], [0-alpha, beta*max(Np63_t)-alpha], 'Color','g','LineStyle', ':','LineWidth',2)
text(0.5*max(Np63_t), beta*0.5*max(Np63_t)-alpha, ' \beta \times Np63(t)-\alpha', 'Color','g');%, 'FontSize',20)

line([0, max(Np63_t)], [-alpha,-alpha], 'Color','r','LineStyle', '--','LineWidth',2)
text(0,-alpha, ' -\alpha', 'Color','r');%, 'FontSize',20)

%
line([k_Np63, k_Np63], [-1,10], 'Color','m','LineStyle', '--','LineWidth',2)
text(k_Np63, -alpha+1, ' k_N_p_6_3', 'Color','m');%, 'FontSize',20)

%
line([0, max(Np63_t)], [beta*k_Np63-.5*alpha,beta*k_Np63-.5*alpha],'Color','c','LineStyle', '--','LineWidth',2) %linea arriba
text(10, beta*k_Np63-.5*alpha, ' \beta*k_N_p_6_3-\alpha/2', 'Color','c');%, 'FontSize',20)

%
line([0, max(Np63_t)], [0, beta*max(Np63_t)], 'Color',[0.9290 0.6940 0.1250],'LineStyle', ':','LineWidth',2)
text(0,1, ' \beta \times Np63(t)', 'Color',[0.9290 0.6940 0.1250]);%, 'FontSize',20)

%
fontsize(gcf,scale=2)


%%
%perpendicular line
perpendicular=@(x)-1/beta*x+k_Np63*(beta+1/beta)-alpha/2


line(Np63_t, -1/beta.*Np63_t, 'Color','b','LineStyle', ':','LineWidth',2)

line(Np63_t, perpendicular(Np63_t), 'Color','b','LineStyle', ':','LineWidth',2)

%line([0, alpha*beta], [0, -alpha], 'Color','b','LineStyle', ':','LineWidth',2)

%%
mu=k_Np63*(beta+1/beta)-alpha/2;
scatter(-mu/(-beta-1/beta),-beta*mu/(-beta-1/beta), 50, 'k')
text(-mu/(-beta-1/beta),-beta*mu/(-beta-1/beta), ' A', 'Color','k', 'FontSize',20)

%syms beta x k_Np63 alpha
scatter((-mu-alpha)/(-beta-1/beta),(alpha/beta-beta*mu)/(-beta-1/beta), 50, 'k')
text((-mu-alpha)/(-beta-1/beta),(alpha/beta-beta*mu)/(-beta-1/beta), ' B', 'Color','k', 'FontSize',20)

%-1/beta*x+k_Np63*(beta+1/beta)-alpha/2

%%
syms beta x k_Np63 alpha mu
%mu=k_Np63*(beta+1/beta)-alpha/2;
Ax=-mu/(-beta-1/beta)
Ay=-beta*mu/(-beta-1/beta)


Bx=(-mu-alpha)/(-beta-1/beta)
By=(alpha/beta-beta*mu)/(-beta-1/beta)


pretty([Ax, Ay])


%% eucclidean distance

dist=sqrt((Bx-Ax)^2+(By-Ay)^2)



simplify((((beta*mu - alpha/beta)/(beta + 1/beta) - (beta*mu)/(beta + 1/beta))^2 + ((alpha + mu)/(beta + 1/beta) - mu/(beta + 1/beta))^2)^(1/2))

(alpha^2/(beta^2 + 1))^(1/2)