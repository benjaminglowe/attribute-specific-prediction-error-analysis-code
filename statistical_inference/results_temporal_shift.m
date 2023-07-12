% Code written to generate Fig. 4B.
%
% author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

%% Loading in analysis parameters
clear
parameters_path = "/Volumes/DATA 1/PhD/EEG_decode_v2/code/statistical_inference";                                                       % path to parameters.mat
cd(parameters_path)
load parameters.mat

%% Defining important variables
path = "/Users/benjaminlowe/Datasets/PhD/EEG_decode_v2/derivatives/results/MVPA";
cd(path) 
data_files = dir('sub*.mat');                                               % subject files
N = size(data_files, 1);                                                    % number of subjects
boot_perms = 100000;                                                        % number of boot permutations
N1_win = 50:73;                                                             % N170 window
dissociation_window = 73:99;                                               % dissociation window
windows = {N1_win; dissociation_window};
TL_fields = ["ori_TL_siz"; "siz_TL_ori"];                                   % transfer learning data structure field names 
x_axis_lim = 20;
subplot_titles = ['Train: Orientation Violation', newline, 'Generalise: Size Violation'; ...
                  'Train: Size Violation', newline, 'Generalise: Orientation Violation'];
time_points = 141;
titles = {'120-260 ms'; '260-400 ms'};
smooth_n = 15;

%% Aggregating data
TGMs = struct();
TGMs.siz_TL_ori = zeros(N, time_points, time_points);
TGMs.ori_TL_siz = zeros(N, time_points, time_points);

for n = 1:N
    load(data_files(n).name)
    
    siz_TL_ori = squeeze(mean(TL(:, 1, 2, :, :), 1)).*percentage_scaler;
    siz_TL_ori = filter2(vert_conv, siz_TL_ori, 'same');
    siz_TL_ori = filter2(hori_conv, siz_TL_ori, 'same');
    siz_TL_ori = siz_TL_ori(TGM_crop, TGM_crop);

    ori_TL_siz = squeeze(mean(TL(:, 2, 1, :, :), 1)).*percentage_scaler;
    ori_TL_siz = filter2(vert_conv, ori_TL_siz, 'same');
    ori_TL_siz = filter2(hori_conv, ori_TL_siz, 'same');
    ori_TL_siz = ori_TL_siz(TGM_crop, TGM_crop);


    TGMs.siz_TL_ori(n, :, :) = siz_TL_ori;
    TGMs.ori_TL_siz(n, :, :) = ori_TL_siz;
end

%% Calculating temporal biases
data = zeros(N, size(TL_fields, 1), size(windows, 1));
t_shifts = struct();
for f = 1:size(TL_fields, 1)
    these_t_shifts = zeros(N, size(windows, 1));% , boot_perms);
        these_TGMs = TGMs.(TL_fields(f, 1));
        for n = 1:N
            part_TGM = squeeze(these_TGMs(n, :, :));
            for w = 1:size(windows, 1)
                this_win = windows{w, 1};
                perm_window = part_TGM(this_win, this_win);
                upper = mean(triu(perm_window, 1), 'all');
                lower = mean(tril(perm_window, -1), 'all');
                difference = upper-lower;
                data(n, f, w) = difference;
            end      
        end
    t_shifts.(TL_fields(f, 1)) = squeeze(mean(these_t_shifts, 1));
end

%% Organising data
data_1 = squeeze(data(:, 1, 1));
data_2 = squeeze(data(:, 1, 2));
data_3 = squeeze(data(:, 2, 1));
data_4 = squeeze(data(:, 2, 2));

all_data = [data_1, data_2, data_3, data_4]';
BFs = bayesfactor(all_data);
errorbars = zeros(1, 4);
for i = 1:4
    win_data = all_data(i, :);
    out = std(win_data);
    errorbars(1, i) = out/sqrt(N);
end

%% Plotting data
close all
figure('Position', [10, 10, 200, 720])
inds = [1, 2; 3, 4];
tlo = tiledlayout(2, 1);

colours = [0.5, 0.5, 0.5; 1, 0, 0];

annotation('textbox', [0.04, 0.98, 0, 0], 'string', 'B', 'FontSize', 25)


lims = [-0.22, 0.075; -0.075, 0.22];

for i = 1:size(inds, 1)
    
    h(i) = nexttile(tlo);

    ind = inds(i, :);
    plot_data = all_data(ind, :);
    b = bar(mean(plot_data, 2), 'FaceColor', 'flat');
    b.CData(1, :) = [0.5, 0.5, 0.5];
    b.CData(2, :) = [1 0 0];
    hold on
    e = errorbar(mean(plot_data, 2), errorbars(:, ind), '-k', ...
        'LineWidth', 1.2);
    e.LineStyle = 'none';
    hold on
    ylabel('Temporal Bias (%)')
    
    these_BFs = BFs(ind, 1);
    for j = 1:size(these_BFs,1)
        disp(these_BFs(j, 1))
    end

    xticklabels(["150 to 265 ms", "265 to 395 ms"])

    set(findobj(gcf,'type','axes'), 'FontSize', 12)
  
    ylim(lims(i, :))
    
    xtickangle(h(i), 90)

    box off
end
