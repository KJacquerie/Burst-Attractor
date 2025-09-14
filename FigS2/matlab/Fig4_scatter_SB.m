clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "PairBased"; 
region = "HPC"; 
bound_type = "SB";
plot_w=0; 
Fig_num = "Fig4_Model"; 


Vbin = load(sprintf('%s/%s/data/%s/%s/%s/Vbin.dat', directory_name,Fig_num, experiment_name, region, bound_type));
Vspk = load(sprintf('%s/%s/data/%s/%s/%s/Vspk.dat', directory_name,Fig_num, experiment_name, region, bound_type));
if(plot_w==1)
    w = load(sprintf('%s/%s/data/%s/%s/%s/w.dat',  directory_name,Fig_num, experiment_name, region, bound_type));
end
w_BACK = load(sprintf('%s/%s/data/%s/%s/%s/w_BACK.dat',  directory_name,Fig_num, experiment_name, region, bound_type));
tmp =  load('idx_w.mat'); % this was saved from Fig 1
idx_w = tmp.idx_w;
points = [1 19 32 47 80] ;


%%



Ncycles=2; 
Tcycle=40000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;


color_tonic = [253 236 214]./255; 
color_burst = [157 195 236]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_hist = [0.7 0.7 0.7]; 
% color rainbow
color_p = [200 52 93; 
            255 193 0;
            112 173 71;
            237 125 49
            91 155 213; 
            ]./255;

        
ptx_long = 15; 
ptx_short = ptx_long/2; 
pty = 2.3; 
ptx = 2.7;
ptx_scatter=3;

%% w

Mks = 4.5; 

smax_tonic = max(max(w_BACK(:,1:2:3)));
smin_tonic = min(min(w_BACK(:,1:2:3)));
    
    
figure
hold on 
plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
plot((w_BACK(:,1)'-smin_tonic)/(smax_tonic-smin_tonic), (w_BACK(:,3)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', 1,'MarkerEdgeColor',color_hist,'MarkerFaceColor',color_hist); 
for idx_p=1:1:length(points)
    plot((w_BACK(idx_w(points(idx_p)),1)'-smin_tonic)/(smax_tonic-smin_tonic), (w_BACK(idx_w(points(idx_p)),3)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
end
xlim([-0.1 1.1])
ylim([-0.1 1.1])
xticks([0 1])
xticklabels({'',''})
yticks([0 1])
yticklabels({'',''})
box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_tonic',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')

%%

smax_burst = max(max(w_BACK(:,2:2:end)));
smin_burst = min(min(w_BACK(:,2:2:end)));
figure
hold on 
plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
plot((w_BACK(:,2)'-smin_burst)/(smax_burst-smin_burst), (w_BACK(:,4)'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', 1,'MarkerEdgeColor',color_hist,'MarkerFaceColor',color_hist); 
for idx_p=1:1:length(points)
    plot((w_BACK(idx_w(points(idx_p)),2)'-smin_burst)/(smax_burst-smin_burst), (w_BACK(idx_w(points(idx_p)),4)'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
end
xlim([-0.1 1.1])
ylim([-0.1 1.1])
xticks([0 1])
xticklabels({'',''})
yticks([0 1])
yticklabels({'',''})
box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_burst',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')


slim = [smin_tonic smin_burst; smax_tonic smax_burst]; 
csvwrite(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_slim',Fig_num, experiment_name, region, bound_type), slim)
