clc
close all 
clear all 

pt=11;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%

color_blue = [106 153 208]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_pink = [250/255 244/255 247/255];

color_corr = [63/255 92/255 206/255]; 
color_uncorr = [130/255 187/255 255/255]; 
color_green = [112 173 71]./255; 
%%

TAG =1; 
switch TAG
    case 0
        expm_name = "NMOD"; 
    case 1
        expm_name = "TAG";
end
w = load(sprintf('/Users/kathleen/Documents/GitHub/Burst-Attractor/Fig4/data/w_%s.dat',expm_name ));


%%

figure
count=1; 
    v = [15000 0; 30000 0;30000 1; 15000 1];
    v2 = [45000 0; 60000 0;60000 1; 45000 1];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices', v, 'FaceColor', color_pink, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', v2, 'FaceColor', color_pink, 'EdgeColor', 'none')
    
for idx=1:1:size(w,1)
    hold on
    if(w(idx,:)~=zeros(1,60000))
     if(count==1 || count==3)
         color_=color_corr; 
     end
     
     if(count>3)
         color_= color_uncorr; 
     end
     %if(count==2)
     %    color_= color_green;
     %end

     plot(w(idx,:), 'color', color_, 'linewidth', 2 ); 
     count=count+1;
    end
end

box off
axis off
%xlim([-80 80])
ylim([0 1])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8.3 2.5]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig5/w_%s', expm_name),'-depsc', '-painters')
