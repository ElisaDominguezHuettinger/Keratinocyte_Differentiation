close all
clearvars
clc

%% Load data
Names_of_genes=load('Names_of_genes.mat');
Names_of_genes=Names_of_genes.Names_of_genes;
selected_genes_normalized=load('selected_genes_normalized.mat');
selected_genes_normalized=selected_genes_normalized.selected_genes_normalized;
selectedgenesraw=load('selected_genes_raw.mat');
selectedgenesraw=selectedgenesraw.selectedgenesraw;
t_exp=load('times_in_hours.mat');
t_exp=t_exp.times_in_hours;

%% Variables of interest
Ap = {'CSTA', 'CXCL10', 'CXCL11', 'RNASE7', 'S100A7','SLPI'};
Bp = {'FLG','IVL'};

    
%% Plot each gene available


color=['r', 'b', 'g', 'm', 'c', 'y', 'k', [0.4940 0.1840 0.5560]];

%% Plot Ap 

figure('Name','Ap1')
for i=1:length(Ap)
    %Set parameters and data
    index_i=find(strcmp(Names_of_genes,Ap(i)));
    data=selectedgenesraw(index_i,:);
    
    data=data./max(data);
    
    % Plot data 
    hold on 
    scatter(t_exp, data,150, color(i), 'd') 
end
%%
index_Ai=[];
    for i=1:length(Ap)
        % Optimize over Ap
        index_Ai=[index_Ai find(strcmp(Names_of_genes,Ap(i)))];
    
    end
mean_AMP=mean(selectedgenesraw([ index_Ai],:));

mean_AMP=mean_AMP./max(mean_AMP);
scatter(t_exp, mean_AMP,250, 'd', 'k', 'filled') 

%%




axis square
xlim([0, 45]);
xlabel('time (hours post-Ca challenge')
legend('-DynamicLegend');
legend('show');
drawnow;
%hold off


%% Plot Bp 
%figure('Name','Bp1')
for i=1:length(Bp)
    %Set parameters and data
    index_i=find(strcmp(Names_of_genes,Bp(i)));
    data=selectedgenesraw(index_i,:);
    data=data./max(data);
    
    hold on 
    
    if i==1;
    scatter(t_exp, data,250, color(i), 'o', 'filled') 
    else
        scatter(t_exp, data,250, color(i), 's', 'filled') 
    end;
end

legend([Ap 'mean AMP' Bp] )

axis square
xlim([0, 45]);
xlabel('time (hours post-Ca challenge')
ylabel('Terminal differentiation markers (normalized to max)')
legend('-DynamicLegend');
legend('show');
drawnow;
hold off
%%
title('Toufighi et al., 2015; FLG, from 0.05 to 1.2 mM CaCl2')


