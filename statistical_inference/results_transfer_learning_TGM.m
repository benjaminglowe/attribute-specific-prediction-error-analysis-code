% Code written to generate Fig. 4A. Reviewer request cells collectively
% generate the TGM contrast analysis within the Supplementary Material.
% Commented lines are basically parameters for generating variants of the
% figure which constitute its supplements.
%
% author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

%% Loading in analysis parameters
clear
parameters_path = "";                                                       % path to parameters.mat
cd(parameters_path)
load parameters.mat

%% Defining important variables
path = "";                                                                  % path to MVPA data within derivatives
cd(path)   
data_files = dir('sub*.mat');                                               % subject files
N = size(data_files, 1);                                                    % number of subjects                                                          % 
TGM_fields = ["ori_TL_siz"; "siz_TL_ori"];
time_points = 141;
subplot_titles = ["Train: Orientation Expectation", "Generalise: Size Expectation"; ...
                  "Train: Size Expectation", "Generalise: Orientation Expectation"];
N1_onset = 50;                                                              % control window onset
N1_offset = 73;                                                             % control window offset

s_onset = 73;                                                               % dissociation window onset
s_offset = 99;                                                              % dissociation window offset
BF_crit = 6;

axis_titles = {'Orientation Expectation', 'Size Expectation'; ...
               'Size Expectation', 'Orientation Expectation'};

%% Aggregating data
TGMs = struct();

for n = 1:N
    load(data_files(n).name)
    
    siz_TL_ori = squeeze(mean(TL(:, 1, 2, :, :), 1)).*percentage_scaler;
    ori_TL_siz = squeeze(mean(TL(:, 2, 1, :, :), 1)).*percentage_scaler;
    
    diff = ori_TL_siz - siz_TL_ori;

    siz_TL_ori = filter2(vert_conv, siz_TL_ori, 'same');
    siz_TL_ori = filter2(hori_conv, siz_TL_ori, 'same');

    ori_TL_siz = filter2(vert_conv, ori_TL_siz, 'same');
    ori_TL_siz = filter2(hori_conv, ori_TL_siz, 'same');

    diff = filter2(vert_conv, diff, 'same');
    diff = filter2(hori_conv, diff, 'same');

    siz_TL_ori = siz_TL_ori(TGM_crop, TGM_crop);
    ori_TL_siz = ori_TL_siz(TGM_crop, TGM_crop);
    diff = diff(TGM_crop, TGM_crop);

    TGMs.siz_TL_ori(n, :, :) = siz_TL_ori;
    TGMs.ori_TL_siz(n, :, :) = ori_TL_siz;
    TGMs.diff(n, :, :) = diff;
end

%% Calculating Bayes Factors
BFs = struct();
tic
for i = 1:size(TGM_fields, 1)
    this_BF = TGMs.(TGM_fields(i, 1));
    this_BF = reshape(this_BF, N, time_points^2);
    this_BF = bayesfactor(this_BF', 'nullvalue', chance, ...
                                    'interval', [-Inf, Inf], ...
                                    'rscale', rscale);
    this_BF = reshape(this_BF, time_points, time_points);
    BFs.(TGM_fields(i, 1)) = this_BF;
end
toc

save BF_TGMs BFs

%% Extracting TGM cells meeting alternative support criterion
load BF_TGMs 
BF_alt = struct();
BF_nul = struct();
for i = 1:size(TGM_fields, 1)
    this_BF = BFs.(TGM_fields(i, 1));
    this_BF = reshape(this_BF, 1, size(TGM_crop, 2)^2);
    alt_output = zeros(size(this_BF));
    nul_output = zeros(size(this_BF));
    for j = 1:size(this_BF, 2)
        if this_BF(1, j) > BF_crit
            alt_output(1, j) = 1;
        elseif this_BF(1, j) < 1/BF_crit
            nul_output(1, j) = 1;
        end
    end
    alt_output = reshape(alt_output, size(TGM_crop, 2), size(TGM_crop, 2));
    nul_output = reshape(nul_output, size(TGM_crop, 2), size(TGM_crop, 2));
    BF_alt.(TGM_fields(i, :)) = alt_output;
    BF_nul.(TGM_fields(i, :)) = nul_output;
end

%% Initialising Figure
close all
figure('Position', [10 10 400 720]);
% figure('Position', [10, 10, 720, 400])
tlo = tiledlayout(2, 1);
% tlo = tiledlayout(1, 2);

annotation('textbox', [0.04, 0.98, 0, 0], 'string', 'A', 'FontSize', 25)

% TGM subplots
for i = 1:size(TGM_fields, 1)
    h(i) = nexttile(tlo);
    TGM_plot = squeeze(mean(TGMs.(TGM_fields(i, 1)), 1));
    % TGM_plot = BFs.(TGM_fields(i, 1));
    clims = [49, 51];
    % clims = [48.5, 51.5];
    % clims = [];
    % clims = [1/100, 100];
    imagesc(h(i), TGM_plot, clims); hold on
    xticks(0:20:size(TGM_crop, 2))                                                % x-axis ticks
    yticks(0:20:size(TGM_crop, 2))                                                % y-axis ticks
    xticklabels(-100:100:600)                                               % x-axis tick labels
    yticklabels(-100:100:600)                                               % y-axis tick labels
    %these_axis_titles = axis_titles{i, :};
    xlabel(['Generalisation Time (ms):', newline, axis_titles{i, 2}], ...  % x-axis label
            'FontSize', 14)                      
    ylabel(['Training Time (ms):', newline, axis_titles{i, 1}], ... % y-axis label
            'FontSize', 14)                            
    colormap(colorcet('D1A'))
    
    set(gca,'ColorScale','log')

%     title(subplot_titles(i, :), 'FontSize', 14)
    
    % Outline clusters
    this_BF_alt = BF_alt.(TGM_fields(i, 1));
    [x, y, ~] = find(this_BF_alt==1);
    for j = 1:size(x, 1)    
        x_ind = x(j, 1);
        y_ind = y(j, 1);
        try
            boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'N');
        catch
        end
        try
            boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'S');
        catch
        end
        try 
            boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'E');
        catch
        end
        try
            boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'W');
        catch
        end
        plot(boundary(:, 2), boundary(:, 1), '-k', 'LineWidth', 1.2)
    end

    this_BF_nul = BF_nul.(TGM_fields(i, 1));
    [x, y, ~] = find(this_BF_nul==1);
    for j = 1:size(x, 1)    
        x_ind = x(j, 1);
        y_ind = y(j, 1);
        try
            boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'N');
        catch
        end
        try
            boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'S');
        catch
        end
        try 
            boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'E');
        catch
        end
        try
            boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'W');
        catch
        end
        plot(boundary(:, 2), boundary(:, 1), '-b', 'LineWidth', 1.2)
    end
    
    plot(1:size(TGM_crop, 2), '-k', 'LineWidth', 1.2)                       % plotting diagonal inside the TGM
    patch([s_offset, s_onset, s_onset, s_offset], ...
          [s_onset, s_onset, s_offset, s_offset], 'red', ...
          'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2); hold on
    patch([N1_offset, N1_onset, N1_onset, N1_offset], ...
          [N1_onset, N1_onset, N1_offset, N1_offset], 'red', ...
          'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2);

end

cbh = colorbar;
cbh.Layout.Tile = 'west'; 
% cbh.Layout.Tile = 'south'; 
% cbh.Label.String = 'Classification Accuracy (%)';
cbh.Label.String = 'Bayes Factor';
cbh.Label.FontSize = 14;

set(findobj(gcf,'type','axes'), 'FontSize', 12)

%% Calculating Bayes factors (reviewer request)
% TGMs.diff = reshape(TGMs.diff, N, time_points^2);
% BFs.diff = bayesfactor(TGMs.diff', 'nullvalue', 0, ...
%                                    'interval', [-Inf, Inf], ...
%                                    'rscale', rscale);
% BFs.diff = reshape(BFs.diff, time_points, time_points);
% % save BF_TGMs BFs
% 
% TGMs.diff = reshape(TGMs.diff, N, time_points, time_points);

%% Categorising BFs (reviewer request)
% % load BF_TGMs
% BF_alt.diff = nan(size(TGM_crop,1),size(TGM_crop,1));
% BF_nul.diff = nan(size(TGM_crop,1),size(TGM_crop,1));
% 
% for i = 1%:size(TGM_fields, 1)
%     this_BF = BFs.diff;
%     this_BF = reshape(this_BF, 1, size(TGM_crop, 2)^2);
%     alt_output = zeros(size(this_BF));
%     nul_output = zeros(size(this_BF));
%     for j = 1:size(this_BF, 2)
%         if this_BF(1, j) > BF_crit
%             alt_output(1, j) = 1;
%         elseif this_BF(1, j) < 1/BF_crit
%             nul_output(1, j) = 1;
%         end
%     end
%     alt_output = reshape(alt_output, size(TGM_crop, 2), size(TGM_crop, 2));
%     nul_output = reshape(nul_output, size(TGM_crop, 2), size(TGM_crop, 2));
%     BF_alt.diff = alt_output;
%     BF_nul.diff = nul_output;
% end


%% Plotting (reviewer request)
% close all
% figure()
% TGM_plot = squeeze(mean(TGMs.diff, 1));
% clims = [-1.2, 1.2];
% imagesc(TGM_plot, clims); hold on
% xticks(0:20:size(TGM_crop, 2))                                              % x-axis ticks
% yticks(0:20:size(TGM_crop, 2))                                              % y-axis ticks
% xticklabels(-100:100:600)                                                   % x-axis tick labels
% yticklabels(-100:100:600)  
% colormap(colorcet('D1A'))
% 
% this_BF_alt = BF_alt.diff;
% [x, y, ~] = find(this_BF_alt==1);
% for j = 1:size(x, 1)    
%     x_ind = x(j, 1);
%     y_ind = y(j, 1);
%     try
%         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'N');
%     catch
%     end
%     try
%         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'S');
%     catch
%     end
%     try 
%         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'E');
%     catch
%     end
%     try
%         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'W');
%     catch
%     end
%     plot(boundary(:, 2), boundary(:, 1), '-k', 'LineWidth', 1.2)
% end
% 
% 
% this_BF_nul = BF_nul.diff;
% [x, y, ~] = find(this_BF_nul==1);
% for j = 1:size(x, 1)    
%     x_ind = x(j, 1);
%     y_ind = y(j, 1);
%     try
%         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'N');
%     catch
%     end
%     try
%         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'S');
%     catch
%     end
%     try 
%         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'E');
%     catch
%     end
%     try
%         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'W');
%     catch
%     end
%     plot(boundary(:, 2), boundary(:, 1), '-b', 'LineWidth', 1.2)
% end
% 
% annotation('textbox', [0.04, 0.98, 0, 0], 'string', 'A', 'FontSize', 25)
% 
% plot(1:size(TGM_crop, 2), '-k', 'LineWidth', 1.2)                           % plotting diagonal inside the TGM
% patch([s_offset, s_onset, s_onset, s_offset], ...
%       [s_onset, s_onset, s_offset, s_offset], 'red', ...
%       'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2); hold on
% patch([N1_offset, N1_onset, N1_onset, N1_offset], ...
%       [N1_onset, N1_onset, N1_offset, N1_offset], 'red', ...
%       'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2);
% 
% cbh = colorbar;
% 
% 
% xlabel('Generalisation Time (ms)','FontSize', 14)                      
% ylabel('Training Time (ms)', 'FontSize', 14)  
% 
% cbh.Label.String = 'Difference in Classification Accuracy (%)';
% cbh.Label.FontSize = 14;
% 
% set(findobj(gcf,'type','axes'), 'FontSize', 12)

%% Plotting Bayes factors (reviewer request)
% figure(2)
% TGM_plot = BFs.diff;% squeeze(mean(TGMs.diff, 1));
% clims = [1/10, 10];
% imagesc(TGM_plot, clims); hold on
% xticks(0:20:size(TGM_crop, 2))                                              % x-axis ticks
% yticks(0:20:size(TGM_crop, 2))                                              % y-axis ticks
% xticklabels(-100:100:600)                                                   % x-axis tick labels
% yticklabels(-100:100:600)  
% colormap(colorcet('D1A'))
% 
% annotation('textbox', [0.04, 0.98, 0, 0], 'string', 'B', 'FontSize', 25)
% 
% % this_BF_alt = BF_alt.diff;
% % [x, y, ~] = find(this_BF_alt==1);
% % for j = 1:size(x, 1)    
% %     x_ind = x(j, 1);
% %     y_ind = y(j, 1);
% %     try
% %         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'N');
% %     catch
% %     end
% %     try
% %         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'S');
% %     catch
% %     end
% %     try 
% %         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'E');
% %     catch
% %     end
% %     try
% %         boundary = bwtraceboundary(this_BF_alt, [x_ind, y_ind], 'W');
% %     catch
% %     end
% %     plot(boundary(:, 2), boundary(:, 1), '-k', 'LineWidth', 1.2)
% % end
% 
% 
% % this_BF_nul = BF_nul.diff;
% % [x, y, ~] = find(this_BF_nul==1);
% % for j = 1:size(x, 1)    
% %     x_ind = x(j, 1);
% %     y_ind = y(j, 1);
% %     try
% %         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'N');
% %     catch
% %     end
% %     try
% %         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'S');
% %     catch
% %     end
% %     try 
% %         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'E');
% %     catch
% %     end
% %     try
% %         boundary = bwtraceboundary(this_BF_nul, [x_ind, y_ind], 'W');
% %     catch
% %     end
% %     plot(boundary(:, 2), boundary(:, 1), '-b', 'LineWidth', 1.2)
% % end
% 
% 
% plot(1:size(TGM_crop, 2), '-k', 'LineWidth', 1.2)                           % plotting diagonal inside the TGM
% patch([s_offset, s_onset, s_onset, s_offset], ...
%       [s_onset, s_onset, s_offset, s_offset], 'red', ...
%       'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2); hold on
% patch([N1_offset, N1_onset, N1_onset, N1_offset], ...
%       [N1_onset, N1_onset, N1_offset, N1_offset], 'red', ...
%       'FaceAlpha', 0, 'LineStyle', '--', 'LineWidth', 1.2);
% 
% cbh = colorbar;
% set(gca,'ColorScale','log')
% 
% xlabel('Generalisation Time (ms)','FontSize', 14)                      
% ylabel('Training Time (ms)', 'FontSize', 14)  
% 
% cbh.Label.String = 'Bayes Factor';
% cbh.Label.FontSize = 14;
% 
% set(findobj(gcf,'type','axes'), 'FontSize', 12)
