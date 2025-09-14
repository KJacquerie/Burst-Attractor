clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/GitHub/Burst-Attractor"; 
experiment_name = "Graupner2012_TUNE"; 
bound_type = "SB";
fMax=50; 
new_curves=1; 
Fig_num = "Fig3"; 

Ncycles=1; 
Tcycle=80000; 
Nstates=4; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;

%tmp =  load(sprintf('%s/%s/data/%s/idx_w.mat', directory_name, Fig_num, experiment_name));
%idx_w = tmp.idx_w;
idx_w = load('idx_wl.dat');


  


  w_rand= load(sprintf('%s/%s/data/%s/w.dat', directory_name, Fig_num, experiment_name));
  %wtemp = w_rand(idx_wl, :);
  %w_rand =wtemp; 
  w_STATE= load(sprintf('%s/%s/data/%s/w_BACK.dat', directory_name, Fig_num, experiment_name));
  
 %%
Vbin = load(sprintf('%s/%s/data/%s/Vbin.dat',directory_name, Fig_num, experiment_name));
Vspk = load(sprintf('%s/%s/data/%s/Vspk.dat',directory_name, Fig_num, experiment_name));

if(bound_type=="SB")
    %V1 = load(sprintf('%s/%s/data/%s/V1.dat',directory_name, Fig_num, experiment_name));
    %V2 = load(sprintf('%s/%s/data/%s/V2.dat',directory_name, Fig_num, experiment_name));
    LFP= load(sprintf('%s/%s/data/%s/LFP.dat',directory_name, Fig_num, experiment_name));
end
gion= load(sprintf('%s/%s/data/%s/gion.dat',directory_name, Fig_num, experiment_name));
gsyn= load(sprintf('%s/%s/data/%s/gsyn.dat',directory_name, Fig_num, experiment_name));

SPB_tonic= load(sprintf('%s/%s/data/%s/SPB_tonic.dat',directory_name, Fig_num, experiment_name));
SPB_burst= load(sprintf('%s/%s/data/%s/SPB_burst.dat',directory_name, Fig_num, experiment_name));
PER_tonic= load(sprintf('%s/%s/data/%s/PER_tonic.dat',directory_name, Fig_num, experiment_name));
PER_burst= load(sprintf('%s/%s/data/%s/PER_burst.dat',directory_name, Fig_num, experiment_name));
IBF_tonic= load(sprintf('%s/%s/data/%s/IBF_tonic.dat',directory_name, Fig_num, experiment_name));
IBF_burst= load(sprintf('%s/%s/data/%s/IBF_burst.dat',directory_name, Fig_num, experiment_name));
DC_tonic= load(sprintf('%s/%s/data/%s/DC_tonic.dat',directory_name, Fig_num, experiment_name));
DC_burst= load(sprintf('%s/%s/data/%s/DC_burst.dat',directory_name, Fig_num, experiment_name));
freq_tonic= load(sprintf('%s/%s/data/%s/freq_tonic.dat',directory_name, Fig_num, experiment_name));
freq_burst= load(sprintf('%s/%s/data/%s/freq_burst.dat',directory_name, Fig_num, experiment_name));

%%


color_tonic = [253 236 214]./255; 
color_burst = [157 195 236]./255; 
color_gray = [ 0.3 0.3 0.3]; 

% couleur rainbow
        color_p = [200 52 93; 
            255 193 0;
            112 173 71;
            237 125 49
            91 155 213; 
            ]./255;


ptx_long = 15; 
ptx_short = ptx_long/2; 
pty = 4; 

ptx_hist = 3; 
pty_hist= 3;

color_hist = [0.7 0.7 0.7]; 


%% w
points = [3 19 32 47 100] ;
figure
hold on
for idx_t=1:1:Nstates
    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end
    

%xlim([T1 T2])
mm = min(min(w_STATE)); 
MM = max(max(w_STATE)); 
ylim([mm*0.9 MM*1.1]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('%s/%s/fig/%s/w',directory_name,  Fig_num, experiment_name), '-depsc', '-painters')


mm_MM = [mm MM]; 
csvwrite(sprintf('%s/%s/fig/%s/w_mm_MM',directory_name,  Fig_num, experiment_name), mm_MM)

%% Raster

Tstate = 20000; 


for idx_state=1:2:8
T1 = idx_state*Tstate/2-2000; 
T2 = idx_state*Tstate/2+2000;
figure
count=1;
idx_t = find(Vbin(1,:)==1); 
hold on
for idx=2:1:size(Vbin,1) % start at 2 to remove the inhibitor
    idx_t = find(Vbin(idx,:)==1); 
    
    plot(idx_t, Vbin(idx,idx_t)+count,'o','color', 'k', 'MarkerSize', 0.5)
    count=count-0.01; 
end

xlim([T1 T2])
ylim([1 2])
xticks([0 1000])
xticklabels({'',''})
yticks([0 1])
yticklabels({'',''})

axis off
box off


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long 2.5]);
print(sprintf('%s/%s/fig/%s/Raster%d', directory_name, Fig_num, experiment_name, idx_state), '-depsc', '-painters')



interv = T1/dt:T2/dt;

yy=smooth(LFP,1500);

figure
hold on
plot([0 interv(end)], [0.1 0.1], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.2 0.2], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.3 0.3], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.4 0.4], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 

plot(t(interv)*1e-3, yy(interv), 'color', color_gray, 'linewidth', 0.75) 
ylim([-0.05 0.45])
xlim([(T1-100)/1000 (T2+100)/1000])
yticks([0 0.1])
yticklabels({'',''})
xticks(0)
xticklabels({''})

box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long 1.5]);
print(sprintf('%s/%s/fig/%s/LFP%d', directory_name, Fig_num, experiment_name, idx_state), '-depsc')

end

%% -----------------------------------------
%  
%       ANALYSIS OF THE FIRING ACTIVITY
%  
%  -----------------------------------------

%  ______ SPB
figure
for idx=1:1:size(SPB_burst,2)
    histogram(SPB_burst(:,idx),'BinEdges', 0.5:1:5.5, 'FaceColor', color_hist)
    xticks([1 2 3 4 5])
    xticklabels({'','','','',''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([1 6])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/SPB%d',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end


%  ______ IBF

figure
for idx=1:1:size(IBF_burst,2)
    histogram(IBF_burst(:,idx),'BinEdges', 0:10:100,  'FaceColor', color_hist)
    xticks([0 100])
    xticklabels({'',''})
    yticks([0 100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([0 100])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/IBF%d',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end

%  ______ PER

figure
for idx=1:1:size(PER_burst,2)
    histogram(1000./PER_burst(:,idx),'BinEdges', 3:0.05:4, 'FaceColor', color_hist)
    xticks([0])
    xticklabels({''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([3 4])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/PER%d',directory_name, Fig_num,experiment_name, idx), '-depsc', '-painters')
end

%  ______ freq

figure
for idx=1:1:size(freq_tonic,2)
    histogram(freq_tonic(:,idx),'BinEdges', 0:10:70, 'FaceColor', color_hist)
    xticks([0 70])
    xticklabels({'',''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([0 70])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/freq%d',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end

%% -----------------------------------------
%           
%               HISTOGRAM-SCATTER
%                   w STATE
%
%  -----------------------------------------

edges = [0:0.025:1];
for idx=1:1:size(w_STATE,2)
    figure
    histogram(w_STATE(:,idx)','BinEdges', edges,'FaceColor', color_hist)
    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1000])
    yticklabels({'',''})
    ylim([-100 1300])
    xlim([0 1])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/whist%d',directory_name, Fig_num, experiment_name,  idx), '-depsc', '-painters')
end


%% Comparison w SCATTER


ptx_scatter=3;

Mks = 5;
count=1;%% -----------------------------------------
%           
%               HISTOGRAM w STATE
%                   normalized
%
%  -----------------------------------------



ptx_scatter=3;

Mks = 5;
count=1;
edges = [0:0.025:1];
for idx=1:2:size(w_STATE,2)
    smax_tonic = max(max(w_STATE(:,idx)));
    smin_tonic = min(min(w_STATE(:,idx)));
    figure
    histogram((w_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 300 ])
    yticklabels({'',''})
    ylim([-100 300])
    xlim([-0.1 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_whist%d',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end
 
for idx=2:2:size(w_STATE,2)
    figure
    smax_burst = max(max(w_STATE(:,idx)));
    smin_burst = min(min(w_STATE(:,idx)));
    histogram((w_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 300])
    yticklabels({'','', ''})
    ylim([-100 300])
    xlim([-0.1, 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_whist%d',directory_name, Fig_num,experiment_name,  idx), '-depsc', '-painters')
end



%% -----------------------------------------
%           
%               SCATTER w STATE
%                   normalized
%
%  -----------------------------------------
idx=5; 
smax_tonic = max(max([w_STATE(:,idx), w_STATE(:,idx+2)]));
smin_tonic = min(min([w_STATE(:,idx), w_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((w_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic), (w_STATE(:,idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Color', color_hist,'Markersize', 1)
    for idx_p=1:1:length(points)
        plot((w_STATE(idx_w(points(idx_p)),idx)'-smin_tonic)/(smax_tonic-smin_tonic), (w_STATE(idx_w(points(idx_p)),idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
    end

    xticks([0 1])
    xticklabels({'', ''})
    yticks([0 1])
    yticklabels({'',''})
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    box off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_scatter*1.5 ptx_scatter*1.5]);
print(sprintf('%s/%s/fig/%s/NORM_ScatterTonic_Article',directory_name, Fig_num,experiment_name), '-depsc', '-painters')



%%
idx=6; 
smax_burst = max(max([w_STATE(:,idx), w_STATE(:,idx+2)]));
smin_burst = min(min([w_STATE(:,idx), w_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((w_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst), (w_STATE(:,idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Color', color_hist,'Markersize', 1)
    for idx_p=1:1:length(points)
        plot((w_STATE(idx_w(points(idx_p)),idx)'-smin_burst)/(smax_burst-smin_burst), (w_STATE(idx_w(points(idx_p)),idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
    end

    xticks([0 1])
    xticklabels({'', ''})
    yticks([0 1])
    yticklabels({'',''})
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    box off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_scatter*1.5 ptx_scatter*1.5]);

print(sprintf('%s/%s/fig/%s/NORM_ScatterBurst_Article',directory_name, Fig_num, experiment_name), '-depsc', '-painters')


slim = [smin_tonic smin_burst; smax_tonic smax_burst]; 
csvwrite(sprintf('%s/%s/fig/%s/slim',directory_name, Fig_num,experiment_name), slim)


