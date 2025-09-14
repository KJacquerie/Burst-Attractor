clear all 
close all
clc

%%
% Load data
Cpre = load('../Cpre.dat');
Cpost = load('../Cpost.dat');
w = load('../w.dat'); 
V = load('../V.dat');

% Folder for saving figures
save_folder = '../fig';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

%% Tonic
T1 = 500; 
T2 = 1.5e3;

plot_and_save(Cpre(1,T1:T2) + Cpost(1,T1:T2), 'tonic_Ca1.svg', 'k')
plot_and_save(Cpre(2,T1:T2) + Cpost(2,T1:T2), 'tonic_Ca2.svg', 'k')
plot_and_save(V(2,T1:T2), 'tonic_V2.svg', [0.3 0.3 0.3])
plot_and_save(V(3,T1:T2), 'tonic_V3.svg', [0.3 0.3 0.3])
plot_and_save(V(4,T1:T2), 'tonic_V4.svg', [0.3 0.3 0.3])
plot_spike_raster(V, 100, 11e3, 'tonic_spikes.pdf')

%% Burst
T1 = 17.43e3; 
T2 = 18.43e3;

plot_and_save(Cpre(1,T1:T2) + Cpost(1,T1:T2), 'burst_Ca1', 'k')
plot_and_save(Cpre(2,T1:T2) + Cpost(2,T1:T2), 'burst_Ca2', 'k')
plot_and_save(V(2,T1:T2), 'burst_V2', [0.3 0.3 0.3])
plot_and_save(V(3,T1:T2), 'burst_V3', [0.3 0.3 0.3])
plot_and_save(V(4,T1:T2), 'burst_V4', [0.3 0.3 0.3])
plot_spike_raster(V, 16e3, 26e3, 'burst_spikes.pdf')


function [spike_times, neuron_ids] = detect_spikes(V, T1, T2)
    spike_times = [];
    neuron_ids = [];
    for i = 1:size(V, 1)
        vtrace = V(i, T1:T2);
        spikes = find(vtrace(1:end-1) < 0 & vtrace(2:end) >= 0); % upward 0mV crossing
        spike_times = [spike_times, spikes + T1]; % relative to global time
        neuron_ids = [neuron_ids, i * ones(1, numel(spikes))];
    end
end

% Plot spike raster and save
function plot_spike_raster(V, T1, T2, fname)
    [spike_times, neuron_ids] = detect_spikes(V, T1, T2);
    figure
    plot(spike_times, neuron_ids, 'k.', 'MarkerSize', 4)
    xlabel('Time (ms)')
    ylabel('Neuron ID')
    format_fig
    print(gcf, fullfile('../fig', fname), '-dpdf', '-painters')
end



% Define figure formatting function
function format_fig()
    set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1)
    box off
    set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 11, 2.5]) % width = 6cm, height = 2cm
    set(gca, 'TickDir', 'out')
end

% Plot and save helper
function plot_and_save(data, filename, color)
    figure
    plot(data, 'Color', color, 'linewidth', 1)
        ylim([-0.2 3.8])
        xlim([1 1000])
    format_fig
    print(gcf, fullfile('../fig', filename), '-dsvg', '-painters')
end

