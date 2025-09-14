clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Deperrois2020_STD"; 
region = "Cortex"; 
bound_type = "HB";
new_curves=0; 
Fig_num = "Fig4_Model"; 


Ncycles=2; 
Tcycle=40000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;
%%
tmp =  load('idx_w.mat'); % this was saved from Fig 1
idx_w = tmp.idx_w;
if(new_curves==1)
    w = load(sprintf('%s/%s/data/%s/%s/%s/w.dat',  directory_name,Fig_num, experiment_name, region, bound_type));

    w_rand = w(idx_w,:); 
    save(sprintf('%s/%s/data/%s/%s/%s/w_rand.mat', directory_name, Fig_num, experiment_name, region, bound_type), 'w_rand');

    wMAX=1;
    wMIN=0; 
    Tstate=20000; 
    
    demo_lim= zeros(size(w,1),Ncycles*2);
    for idx_state = 1:2:Ncycles*2
        for idw = 1:1:size(w,1)
            if(w(idw,idx_state*Tstate)>=wMAX*0.99)
                idx_win = find(w(idw, (idx_state-1)*Tstate+100:idx_state*Tstate)>=wMAX*0.99); 
                demo_lim(idw, idx_state)=(idx_state-1)*Tstate+100+idx_win(1);
            else
                if (w(idw,idx_state*Tstate)<=wMIN*0.99)
                    idx_win = find(w(idw, (idx_state-1)*Tstate+100:idx_state*Tstate)<=wMIN*0.99); 
                    demo_lim(idw, idx_state)=(idx_state-1)*Tstate+100+idx_win(1);
                else
                    demo_lim(idw,idx_state)=idx_state*Tstate; 
                end
            end 
            if(w(idw,idx_state*Tstate)==w(idw, (idx_state-1)*Tstate+100))
                demo_lim(idw,idx_state)=0; 
            end
        end
    end


    for idx_state = 2:2:Ncycles*2
        for idw = 1:1:size(w,1)
            if(w(idw,idx_state*Tstate)>=wMAX*0.99)
                idx_win = find(w(idw, (idx_state-1)*Tstate+1000:idx_state*Tstate)>=wMAX*0.99); 
                demo_lim(idw, idx_state)=(idx_state-1)*Tstate+1000+idx_win(1);
            else
                if (w(idw,idx_state*Tstate)<=wMIN*0.99)
                    idx_win = find(w(idw, (idx_state-1)*Tstate+1000:idx_state*Tstate)<=wMIN*0.99); 
                    demo_lim(idw, idx_state)=(idx_state-1)*Tstate+1000+idx_win(1);
                else
                    demo_lim(idw,idx_state)=idx_state*Tstate; 
                end
            end 
            if(w(idw,idx_state*Tstate)==w(idw, (idx_state-1)*Tstate+500))
                demo_lim(idw,idx_state)=0; 
            end
            
          
        end
    end


    save(sprintf('%s/%s/data/%s/%s/%s/demo_lim.mat', directory_name, Fig_num, experiment_name, region, bound_type), 'demo_lim');

    
else
  demo_lim = load(sprintf('%s/%s/data/%s/%s/%s/demo_lim.mat',directory_name, Fig_num, experiment_name, region, bound_type));
  demo_lim=demo_lim.demo_lim;   
  %w_rand = load(sprintf('%s/%s/data/%s/%s/%s/w_rand.dat',  directory_name,Fig_num, experiment_name, region, bound_type));

end

w_BACK = load(sprintf('%s/%s/data/%s/%s/%s/w_BACK.dat',  directory_name,Fig_num, experiment_name, region, bound_type));

points = [1 19 32 47 80] ;


%%


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

idx_points = idx_w(points); 


%% SLOPE NUMERIQUE

slope = zeros(size(w_BACK,1),Ncycles*2);
nPre= 50; 
nPost=50; 
Tstate =Tcycle/2; 
w0 = 0.5*ones(nPre*nPost, 1);
for idx=1:1:Ncycles*2
    for idw = 1:1:size(slope,1)
        %if(demo_lim(idw,idx)==idx*Tstate)
            %slope(idw,idx) = 0; 
        %else
            slope(idw, idx) = (w_BACK(idw, idx) - w0(idw))./ (demo_lim(idw, idx) - (idx-1)*Tstate);
        %end
        if(demo_lim(idw,idx)==0)
            slope(idw,idx) =0; 
        end
    end
    w0 = w_BACK(:, idx);
end

%% reconstruction

if(new_curves==1)
    figure
    hold on
    plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)


    for idx_p=1:1:length(points)
        plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
        w0 = 0.5*ones(nPre*nPost, 1);
        for idx_state=1:1:Ncycles*2
            plot( [(idx_state-1)*Tstate, demo_lim(idx_w(points(idx_p)),idx_state)], [w0(idx_w(points(idx_p))) w0(idx_w(points(idx_p)))+slope(idx_w(points(idx_p)),idx_state).*(demo_lim(idx_w(points(idx_p)),idx_state)-(idx_state-1)*Tstate)], '-o','color', color_p(idx_p, :),'linewidth', 2)
            w0 = w_BACK(:,idx_state);

        end

    end
end



%% ABSOLUTE   %%%




Mks = 5;
count=1;


idx=1; 
%smax_tonic = max(max(abs(slope(:,1:2:end))));
%smin_tonic = min(min(abs(slope(:,1:2:end))));

slopeNZ1=abs(slope(:,idx)); 
slopeNZ1(slopeNZ1==0)=inf; 

slopeNZ2=abs(slope(:,idx+2)); 
slopeNZ2(slopeNZ2==0)=inf; 


smax_tonic = max(max([abs(slope(:,idx)), abs(slope(:,idx+2))]));
smin_tonic = min(min([slopeNZ1, slopeNZ2]));

NP_tonic =0; 
figure
    hold on 
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    for idw=1:1:size(slope,1)
        if(slope(idw,idx) ~=0 && slope(idw,idx+2)~=0)
            plot((abs(slope(idw,idx))'-smin_tonic)/(smax_tonic-smin_tonic), (abs(slope(idw,idx+2))'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Color', color_hist,'Markersize', 1)
            NP_tonic=NP_tonic+1; 
        end
    end
    for idx_p=1:1:length(points)
        if((abs(slope(idx_w(points(idx_p)),idx)) ~=0 && abs(slope(idx_w(points(idx_p)),idx+2))~=0))
            plot((abs(slope(idx_w(points(idx_p)),idx))'-smin_tonic)/(smax_tonic-smin_tonic), (abs(slope(idx_w(points(idx_p)),idx+2))'-smin_tonic)/(smax_tonic-smin_tonic),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
        end
    end

    xticks([0 1])
    xticklabels({'',''})
    yticks([ 0 1])
    yticklabels({'', '',''})
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    box off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_tonic',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')




%%


idx=2; 
%smax_burst = max(max(abs(slope(:,2:2:end))));
%smin_burst = min(min(abs(slope(:,2:2:end))));

slopeNZ1=abs(slope(:,idx)); 
slopeNZ1(slopeNZ1==0)=inf; 

slopeNZ2=abs(slope(:,idx+2)); 
slopeNZ2(slopeNZ2==0)=inf; 


smax_burst = max(max([abs(slope(:,idx)), abs(slope(:,idx+2))]));
smin_burst = min(min([slopeNZ1, slopeNZ2]));

NP_burst =0; 
figure
hold on 
    plot([-2 2], [-2 2], 'color', [0.5 0.5 0.5])
    for idw=1:1:size(slope,1)
        if((slope(idw,idx) ~=0 && slope(idw,idx+2)~=0))
            plot((abs(slope(idw,idx))'-smin_burst)/(smax_burst-smin_burst), (abs(slope(idw,idx+2))'-smin_burst)/(smax_burst-smin_burst),'o', 'Color', color_hist,'Markersize', 1)
            NP_burst=NP_burst+1; 
        end
    end
    for idx_p=1:1:length(points)
        if((abs(slope(idx_w(points(idx_p)),idx)) ~=0 && abs(slope(idx_w(points(idx_p)),idx+2))~=0))
            plot((abs(slope(idx_w(points(idx_p)),idx))'-smin_burst)/(smax_burst-smin_burst), (abs(slope(idx_w(points(idx_p)),idx+2))'-smin_burst)/(smax_burst-smin_burst),'o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
        end
    end

    xticks([0 1])
    xticklabels({'',''})
    yticks([ 0 1])
    yticklabels({'', '',''})
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    box off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_burst',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')



slim = [smin_tonic smin_burst; smax_tonic smax_burst; NP_tonic NP_burst]; 
csvwrite(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/scatter/%s_%s_%s_slim',Fig_num, experiment_name, region, bound_type), slim)



