clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/GitHub/Burst-Attractor"; 
experiment_name = "Graupner2012_RESET"; 
bound_type = "SB";
fMax=50; 
new_curves=0; 
Fig_num = "Fig1"; 

Ncycles=4; 
Tcycle=40000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;

%tmp =  load(sprintf('%s/%s/data/%s/idx_w.mat', directory_name, Fig_num, experiment_name));
%idx_w = tmp.idx_w;
idx_w = load('idx_wl.dat');


  

  w_rand= load(sprintf('%s/%s/data/%s/w.dat', directory_name, Fig_num, experiment_name));
  w_STATE= load(sprintf('%s/%s/data/%s/w_BACK.dat', directory_name, Fig_num, experiment_name));
  
  l_rand= load(sprintf('%s/%s/data/%s/late_w.dat', directory_name, Fig_num, experiment_name));
  l_STATE= load(sprintf('%s/%s/data/%s/g_BACK.dat', directory_name, Fig_num, experiment_name));
 
  wl_rand = l_rand.*w_rand;
  wl_STATE = w_STATE.*l_STATE;


%%
Vbin = load(sprintf('%s/%s/data/%s/Vbin.dat',directory_name, Fig_num, experiment_name));
Vspk = load(sprintf('%s/%s/data/%s/Vspk.dat',directory_name, Fig_num, experiment_name));

if(bound_type=="SB")
    V1 = load(sprintf('%s/%s/data/%s/V1.dat',directory_name, Fig_num, experiment_name));
    V2 = load(sprintf('%s/%s/data/%s/V2.dat',directory_name, Fig_num, experiment_name));
    LFP= load(sprintf('%s/%s/data/%s/LFP.dat',directory_name, Fig_num, experiment_name));
end
gion= load(sprintf('%s/%s/data/%s/gion.dat',directory_name, Fig_num, experiment_name));
gsyn= load(sprintf('%s/%s/data/%s/gsyn.dat',directory_name, Fig_num, experiment_name));

%%
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
points = [3 19 32 47 100];%[1 19 32 47 80] ;
figure
hold on
for idx_t=1:1:2*Ncycles
    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end
    

%xlim([T1 T2])
mm = min(min(w_STATE)); 
MM = max(max(w_STATE)); 
ylim([mm MM]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
%print(sprintf('%s/%s/fig/%s/w',directory_name,  Fig_num, experiment_name), '-depsc', '-painters')


%mm_MM = [mm MM]; 
%csvwrite(sprintf('%s/%s/fig/%s/w_mm_MM',directory_name,  Fig_num, experiment_name), mm_MM)


%%

figure
hold on
for idx_t=1:1:2*Ncycles
    plot([idx_t*20e3 idx_t*20e3], [0 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
plot(l_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( l_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end
%end
%xlim([T1 T2])
mm = min(min(l_STATE)); 
MM = max(max(l_STATE)); 
ylim([mm-0.00005 MM+0.00005]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
if(mm==MM)
    yticks([mm])
    yticklabels({''})
else
    yticks([mm MM])
    yticklabels({'', ''})
end



set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
%print(sprintf('%s/%s/fig/%s/l',directory_name,  Fig_num, experiment_name), '-depsc', '-painters')

%mm_MM = [mm MM]; 
%csvwrite(sprintf('%s/%s/fig/%s/l_mm_MM',directory_name,  Fig_num, experiment_name), mm_MM)



%%


figure
hold on
for idx_t=1:1:2*Ncycles
    plot([idx_t*20e3 idx_t*20e3], [0 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
plot(wl_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( wl_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end
%end
%xlim([T1 T2])
mm = min(min(wl_STATE)); 
MM = max(max(wl_STATE)); 
ylim([mm MM]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
%print(sprintf('%s/%s/fig/%s/wl',directory_name, Fig_num, experiment_name), '-depsc', '-painters')



%mm_MM = [mm MM]; 
%csvwrite(sprintf('%s/%s/fig/%s/wl_mm_MM',directory_name,  Fig_num, experiment_name), mm_MM)


%% Raster
T1 = Tcycle/2-1500; 
T2 = Tcycle/2+2500;

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
%print(sprintf('%s/%s/fig/%s/Raster', directory_name, Fig_num, experiment_name), '-depsc', '-painters')

%% LFP

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
%print(sprintf('%s/%s/fig/%s/LFP', directory_name, Fig_num, experiment_name), '-depsc')



%% Plot V

    tshort = T1:T2;
    TV1= 500; 
    TV2 = 4500;


    NVplot = 5; 

    figure
    plot(V1(1,TV1:TV2)', 'color', color_gray, 'linewidth', 0.5)
    xticks(0)
    xticklabels({''})
    yticks([-50 0])
    yticklabels({'',''})
    ylim([-90 25])
    xlim([1 TV2-TV1])
    axis off
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_long 0.6]);
    print(sprintf('%s/%s/fig/%s/VInh',directory_name, Fig_num, experiment_name),'-depsc')


    idxV = find(freq_tonic<10 & freq_tonic>0); 
    for idx=1:1:4
    figure(idx+1)
        plot(V1(idxV(idx),TV1:TV2)', 'color', color_gray, 'linewidth', 0.5)
        xticks(0)
        xticklabels({''})
        yticks([-50 0])
        yticklabels({'',''})
        ylim([-90 25])
        xlim([1 TV2-TV1])
        axis off
        box off
        set(gcf,'PaperPositionMode','auto');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 ptx_long 0.6]);
        print(sprintf('%s/%s/fig/%s/V%d',directory_name, Fig_num, experiment_name, idx), '-depsc')
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
    ylim([-100 1200])
    xlim([0 1])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    %print(sprintf('%s/%s/fig/%s/whist%d',directory_name, Fig_num, experiment_name,  idx), '-depsc', '-painters')
end



%% Comparison w SCATTER


% ptx_scatter=3;
% 
% Mks = 5;
% count=1;
% figure
% for idx=1:2:5
%     subplot(1,3,count)
%     hold on 
%     plot([0 1], [0 1], 'color', [0.5 0.5 0.5])
%     plot(w_STATE(:,idx)', w_STATE(:,idx+2)','o', 'Color', color_hist,'Markersize', 1)
%     for idx_p=1:1:length(points)
%         plot(w_STATE(idx_w(points(idx_p)),idx)', w_STATE(idx_w(points(idx_p)),idx+2)','o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
%     end
% 
%     xticks([0 1])
%     xticklabels({'',''})
%     yticks([0 1])
%     yticklabels({'',''})
%     xlim([0 1])
%     ylim([0 1])
%     box off
%     count=count+1;
% end
% 
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx_scatter*4.5 ptx_scatter]);
% print(sprintf('%s/%s/fig/ScatterTonic_subplot',directory_name, Fig_num), '-depsc', '-painters')
% 
% 
% %%
% count=1;
% figure
% for idx=2:2:6
%     subplot(1,3,count)
%     hold on 
%     plot([0 1], [0 1], 'color', [0.5 0.5 0.5])
%     plot(w_STATE(:,idx)', w_STATE(:,idx+2)','o', 'Color', color_hist,'Markersize', 1)
%     for idx_p=1:1:length(points)
%         plot(w_STATE(idx_w(points(idx_p)),idx)', w_STATE(idx_w(points(idx_p)),idx+2)','o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
%     end
% 
%     xticks([0 1])
%     xticklabels({'',''})
%     yticks([0 1])
%     yticklabels({'',''})
%     xlim([0 1])
%     ylim([0 1])
%     box off
%     count=count+1;
% end
% 
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx_scatter*4.5 ptx_scatter]);
% print(sprintf('%s/%s/fig/ScatterBurst_subplot',directory_name, Fig_num), '-depsc', '-painters')
% 
% 
% 
% %%
% figure
% idx=5; 
%     hold on 
%     plot([0 1], [0 1], 'color', [0.5 0.5 0.5])
%     plot(w_STATE(:,idx)', w_STATE(:,idx+2)','o', 'Color', color_hist,'Markersize', 1)
%     for idx_p=1:1:length(points)
%         plot(w_STATE(idx_w(points(idx_p)),idx)', w_STATE(idx_w(points(idx_p)),idx+2)','o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
%     end
% 
%     xticks([0 min(w_STATE(:,idx)) max(w_STATE(:,idx))  1])
%     xticklabels({'', '','', ''})
%     yticks([0 min(w_STATE(:,idx+2)) max(w_STATE(:,idx+2)) 1])
%     yticklabels({'','','',''})
%     xlim([0 1])
%     ylim([0 1])
%     box off
% 
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx_scatter*1.5 ptx_scatter*1.5]);
% print(sprintf('%s/%s/fig/ScatterTonic_Article',directory_name, Fig_num), '-depsc', '-painters')
% 
% 
% 
% 
% %%
% 
% figure
% idx=6; 
%     hold on 
%     plot([0 1], [0 1], 'color', [0.5 0.5 0.5])
%     plot(w_STATE(:,idx)', w_STATE(:,idx+2)','o', 'Color', color_hist,'Markersize', 1)
%     for idx_p=1:1:length(points)
%         plot(w_STATE(idx_w(points(idx_p)),idx)', w_STATE(idx_w(points(idx_p)),idx+2)','o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
%     end
% 
%     xticks([0 min(w_STATE(:,idx)) max(w_STATE(:,idx))  1])
%     xticklabels({'', '','', ''})
%     yticks([0 min(w_STATE(:,idx+2)) max(w_STATE(:,idx+2)) 1])
%     yticklabels({'','','',''})
%     xlim([0 1])
%     ylim([0 1])
%     box off
% 
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx_scatter*1.5 ptx_scatter*1.5]);
% print(sprintf('%s/%s/fig/ScatterBurst_Article',directory_name, Fig_num), '-depsc', '-painters')
% 


%% -----------------------------------------
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
    yticks([0 500 ])
    yticklabels({'',''})
    ylim([-100 500])
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
    yticks([0 500])
    yticklabels({'','', ''})
    ylim([-100 500])
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
idx=1; 
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
idx=2; 
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


%% -----------------------------------------
%           
%               HISTOGRAM wl STATE
%                   normalized
%
%  -----------------------------------------



ptx_scatter=3;

Mks = 5;
count=1;
edges = [0:0.025:1];
for idx=1:2:size(w_STATE,2)
    smax_tonic = max(max(wl_STATE(:,idx)));
    smin_tonic = min(min(wl_STATE(:,idx)));
    figure
    histogram((wl_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 500 1000])
    yticklabels({'',''})
    ylim([-100 500])
    xlim([-0.1 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_wlhist%d',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end
 
for idx=2:2:size(wl_STATE,2)
    figure
    smax_burst = max(max(wl_STATE(:,idx)));
    smin_burst = min(min(wl_STATE(:,idx)));
    histogram((wl_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 500 1000])
    yticklabels({'','', ''})
    ylim([-100 500])
    xlim([-0.1, 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_wlhist%d',directory_name, Fig_num,experiment_name,  idx), '-depsc', '-painters')
end




%% -----------------------------------------
%           
%               SCATTER wl STATE
%                   normalized
%
%  -----------------------------------------



clear gca

idx=5; 
smax_tonic = max(max([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
smin_tonic = min(min([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((wl_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic), (wl_STATE(:,idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Color', color_hist,'Markersize', 1)
    for idx_p=1:1:length(points)
        plot((wl_STATE(idx_w(points(idx_p)),idx)'-smin_tonic)/(smax_tonic-smin_tonic), (wl_STATE(idx_w(points(idx_p)),idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
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
print(sprintf('%s/%s/fig/%s/NORM_ScatterTonic_Article_wl',directory_name, Fig_num,experiment_name), '-depsc', '-painters')




idx=6; 
smax_burst = max(max([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
smin_burst = min(min([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((wl_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst), (wl_STATE(:,idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Color', color_hist,'Markersize', 0.7)
    for idx_p=1:1:length(points)
        plot((wl_STATE(idx_w(points(idx_p)),idx)'-smin_burst)/(smax_burst-smin_burst), (wl_STATE(idx_w(points(idx_p)),idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', Mks*0.7,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
    end

    xticks([0 1])
    xticklabels({'', ''})
    yticks([0 1])
    yticklabels({'',''})
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    box off
    axis off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 2.5 2.5]);

print(sprintf('%s/%s/fig/%s/NORM_ScatterBurst_Article_wl',directory_name, Fig_num, experiment_name), '-depsc', '-painters')


slim = [smin_tonic smin_burst; smax_tonic smax_burst]; 
csvwrite(sprintf('%s/%s/fig/%s/slim_wl',directory_name, Fig_num,experiment_name), slim)


%% -----------------------------------------
%           
%               HISTOGRAM wl STATE
%                   normalized
%
%  -----------------------------------------



ptx_scatter=3;

Mks = 5;
count=1;
edges = [0:0.025:1];
smax_tonic = MM;%max(max(wl_STATE(:,idx)));
smin_tonic = mm;%min(min(wl_STATE(:,idx)));
    
for idx=1:2:size(w_STATE,2)
    
    figure
    histogram((wl_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 600 1000])
    yticklabels({'',''})
    ylim([-100 600])
    xlim([-0.1 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_wlhist%d_1',directory_name, Fig_num, experiment_name, idx), '-depsc', '-painters')
end
 
for idx=2:2:size(wl_STATE,2)
    figure
    smax_burst = MM;%max(max(wl_STATE(:,idx)));
    smin_burst = mm;%min(min(wl_STATE(:,idx)));
    histogram((wl_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst),'BinEdges', edges,'FaceColor', color_hist)
    xticks([0, 1])
    xticklabels({'',''})
    yticks([0 600 1000])
    yticklabels({'','', ''})
    ylim([-100 600])
    xlim([-0.1, 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/NORM_wlhist%d_1',directory_name, Fig_num,experiment_name,  idx), '-depsc', '-painters')
end




%% -----------------------------------------
%           
%               SCATTER wl STATE
%                   normalized
% between the range defined by tonic 
%
%  -----------------------------------------
% FOR THE PAPER I CHOOSE THIS VERSION  



idx=3; 
smax_tonic = MM;%max(max([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
smin_tonic = mm;%min(min([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((wl_STATE(:,idx)'-smin_tonic)/(smax_tonic-smin_tonic), (wl_STATE(:,idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Color', color_hist,'Markersize', 1)
    for idx_p=1:1:length(points)
        plot((wl_STATE(idx_w(points(idx_p)),idx)'-smin_tonic)/(smax_tonic-smin_tonic), (wl_STATE(idx_w(points(idx_p)),idx+2)'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
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
print(sprintf('%s/%s/fig/%s/NORM_ScatterTonic_Article_wl_%d',directory_name, Fig_num,experiment_name, idx), '-dsvg', '-painters')




idx=4; 
smax_burst = smax_tonic; %max(max([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
smin_burst = smin_tonic; % ;%min(min([wl_STATE(:,idx), wl_STATE(:,idx+2)]));
figure
    hold on 
    
    %plot([smin_burst smax_burst], [smin_burst smax_burst], 'color', [0.5 0.5 0.5])
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    
    plot((wl_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst), (wl_STATE(:,idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Color', color_hist,'Markersize', 1)
    for idx_p=1:1:length(points)
        plot((wl_STATE(idx_w(points(idx_p)),idx)'-smin_burst)/(smax_burst-smin_burst), (wl_STATE(idx_w(points(idx_p)),idx+2)'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
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

print(sprintf('%s/%s/fig/%s/NORM_ScatterBurst_Article_wl_%d',directory_name, Fig_num, experiment_name, idx), '-dsvg', '-painters')


slim = [smin_tonic smin_burst; smax_tonic smax_burst]; 
csvwrite(sprintf('%s/%s/fig/%s/slim_wl_1',directory_name, Fig_num,experiment_name), slim)

%% ____ corr coff ___
crcf = zeros(6,1);
rgcf = zeros(6,1);
for idx=1:1:6
wl_X =(wl_STATE(:,idx)'-smin_burst)/(smax_burst-smin_burst);  
wl_Y = (wl_STATE(:,idx+2)'-smin_burst)/(smax_burst-smin_burst);
mat = corrcoef(wl_X, wl_Y);
crcf(idx) = mat(2);

%[~,~,~,~,stats] = regress(wl_X, wl_Y); 
%rgcf(idx) = stats(1);

mdl = fitlm(wl_X,wl_Y);
ftcf(idx) = mdl.Rsquared.Ordinary;



end

%% _______ MEAN and STD _____


sum_TONIC = reshape([w_STATE(:,1) w_STATE(:,3) w_STATE(:,5) w_STATE(:,7)], 2500*4,1);
mean_TONIC  = mean(sum_TONIC);
std_TONIC  = std(sum_TONIC); 

sum_BURST = reshape([w_STATE(:,2) w_STATE(:,4) w_STATE(:,6) w_STATE(:,8)], 2500*4,1);
mean_BURST  = mean(sum_BURST);
std_BURST  = std(sum_BURST); 

mean_std = [mean_TONIC mean_BURST; std_TONIC std_BURST]; 
csvwrite(sprintf('%s/%s/fig/%s/mean_std',directory_name, Fig_num,experiment_name), slim)

%%

mean_state = mean(w_STATE, 1); 
std_state = std(w_STATE); 

figure
hold on
%for idx_t=1:1:2*Ncycles
%    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
%end
errorbar([20e3:20e3:160e3], mean_state, std_state, 'color', [0.5 0.5 0.5])

%xlim([T1 T2])
mm = min(min(w_STATE)); 
MM = max(max(w_STATE)); 
ylim([mm MM]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})
xlim([0 160000])

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('%s/%s/fig/%s/error',directory_name,  Fig_num, experiment_name), '-depsc', '-painters')

