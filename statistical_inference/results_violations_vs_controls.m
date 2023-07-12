% Code written to generate Fig. 2.
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
N = size(data_files, 1);                                                    % number of subjects
onset = 40;                                                                 % time point of stimulus onset
chance = 50;                                                                % theoretical chance (%)
percentage_scaler = 100;                                                    % scaler converting DVs to a percentage
time_points = 181;                                                          % time points within each epoch
chance_vector = ones(1, time_points).*chance;                               % theoretical chance vector for plotting
interval = [0.5, Inf];                                                      % prior range
rscale = 0.7071;                                                            % prior width
smooth_win = 15;
smooth_mthd = 'movmean';                                                    % time-series smoothing method
conditions = ["siz_vio", "ori_vio", "col_vio"];                             % conditions classified against control       
BF_crit = 6;                                                                % critical Bayes factor
BF_cats = ["alt"; "nul"; "inc"];                                            % Bayes factor categories
sz = 10;                                                                    % scatter plot data point size
BF_scatter_lim = 8;                                                         % y-axis limit of BF scatter plots
ylim_values = [49.4, 52.5];                                                 % y-axis limits
subplot_titles = ["Size Violation", ""; "Orientation Violation", ""];

%% Aggregating data
time_series = struct();                                                     % struct for storing time-series data
time_series.col_vio = zeros(N, time_points);                                % col_vios vs. controls time-series
time_series.siz_vio = zeros(N, time_points);                                % siz_vios vs. controls time-series
time_series.ori_vio = zeros(N, time_points);                                % ori_vios vs. controls time-series

TGMs = struct();
TGMs.col_vio = zeros(N, size(TGM_crop, 2), size(TGM_crop, 2));
TGMs.siz_vio = zeros(N, size(TGM_crop, 2), size(TGM_crop, 2));
TGMs.ori_vio = zeros(N, size(TGM_crop, 2), size(TGM_crop, 2));

for n = 1:N
    load(data_files(n).name)
    vio_vs_cont = squeeze(mean(vio_vs_cont, 1)).*percentage_scaler;
    diagonals = zeros(size(conditions, 2), time_points);
    part_TGMs = zeros(size(conditions, 2), size(TGM_crop, 2), ...
                     size(TGM_crop, 2));
    for c = 1:size(conditions, 2)
        diagonals(c, :) = pop_diagonal(squeeze(vio_vs_cont(c, :, :)));
    end
    diagonals = smoothdata(diagonals, 2, smooth_mthd, smooth_win);
  
    time_series.col_vio(n, :) = diagonals(1, :);                             % storing time-series data
    time_series.siz_vio(n, :) = diagonals(2, :);
    time_series.ori_vio(n, :) = diagonals(3, :);

    TGMs.col_vio(n, :, :) = part_TGMs(1, :, :);                             % storing TGMs
    TGMs.siz_vio(n, :, :) = part_TGMs(2, :, :);
    TGMs.ori_vio(n, :, :) = part_TGMs(3, :, :);
end

%% Calculating Bayes factors for time-series data
BF_time_series = struct();                                                  % struct for storing time-series Bayes factors

for c = 1:size(conditions, 2)                                               % for each pairwise comparison between violations and controls
    data = time_series.(conditions(1, c))';                                 % data for this comparison
    BF_time_series.(conditions(1, c)) = bayesfactor(data, ...               % computing Bayes factors
                                                  'nullvalue', chance, ...
                                                  'interval', interval, ...
                                                  'rscale', rscale);
end

%% Categorising Bayes factors
BF_alt = struct();                                                          % structure for storing time points with evidence for the alternative
BF_nul = struct();                                                          % structure for storing time points with evidence for the null
BF_inc = struct();                                                          % structure for storing time points with inconclusive evidence

for c = 1:size(conditions, 2)
    BF_data = BF_time_series.(conditions(1, c));                            % condition specific Bayes factor values
    alt = NaN(time_points, 2);                                              % empty matrix for storing values
    nul = NaN(time_points, 2);
    inc = NaN(time_points, 2);
    for t = 1:time_points                                                   % for loop for categorising Bayes factors
        score = BF_data(t, 1);
        if score > BF_crit
            alt(t, :) = [t, score];
        elseif score < 1/BF_crit
            nul(t, :) = [t, score];
        else
            inc(t, :) = [t, score];
        end
    end
    alt(any(isnan(alt), 2), :) = [];                                        % remove NaNs for each matrix
    nul(any(isnan(nul), 2), :) = [];
    inc(any(isnan(inc), 2), :) = [];
    BF_alt.(conditions(1, c)) = alt;                                        % places Bayes factor arrays into structures defined above
    BF_nul.(conditions(1, c)) = nul;
    BF_inc.(conditions(1, c)) = inc;
end

%% Creating figure
close all
figure('Position', [10, 10, 600, 600]);

% Accuracy subplot
subplot_secs = [1, 2, 3, 4, 5, 6, 7, 8, 9, ...
                11, 12, 13, 14, 15, 16, 17, 18, 19, ...
                21, 22, 23, 24, 25, 26, 27, 28, 29];
subplot(6, 10, subplot_secs);                                               % accuracy (main) subplot
colours = ['-r'; '-b'; '-k'];                                               % colour of each time-series
ignore = plot(chance_vector, '--k'); hold on                                % plotting theoretical chance
ignore.Annotation.LegendInformation.IconDisplayStyle = 'off';               % ignoring theoretical chance for legend
for c = 1:size(conditions, 2)                                               % loop plotting time-series data
    plot_data = time_series.(conditions(1, c));                             % calling data from time-series field
    plot_data = squeeze(mean(plot_data, 1));                                % averaging across participants
    plot_SE = std(time_series.(conditions(1, c)))./sqrt(N);                 % calculating standard error
    [hl, hp] = boundedline(1:time_points, plot_data, plot_SE, ...           % plotting data
                           colours(c, :), 'alpha');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';               % removing error bars from legend
end
xticks(20:20:time_points)                                                   % x-axis ticks
xticklabels(-100:100:700)                                                      % x-axis labels
xlim([0, time_points-1])                                                    % x-axis limit
ylim(ylim_values)                                                           % y-axis limit
ylabel('Classification Accuracy (%)', 'FontSize', 14)                       % y-axis label
lgd = legend('Size Expectation', 'Orientation Expectation', ...                 % plot legend 
             'Brightness Expectation');
lgd.FontSize = 12;                                                          % legend font size
legend boxoff                                                               % removing legend boarder
vline(onset, '--k');                                                        % stimulus onset vertical line   
box off

% Bayes factors subplot
legend_titles = ["Size Exp."; "Ori. Exp."; "Bri. Exp."];
colours_BF = [1, 0, 0; 0, 0, 1; 0, 0, 0];                                   % scatter plot colours
subplot_secs = [31, 32, 33, 34, 35, 36, 37, 38, 39; ...
                41, 42, 43, 44, 45, 46, 47, 48, 49; ...
                51, 52, 53, 54, 55, 56, 57, 58, 59];
lgd_pos = [0.9, 0.455, 0, 0; ...
           0.9, 0.313, 0, 0; ...
           0.9, 0.17, 0, 0];
for c = 1:size(conditions, 2)                                               % Bayes factor scatter subplots
    sub_c = subplot(6, 10, subplot_secs(c, :));                             % defining scatter subplot
    ignore = plot(ones(1, time_points), '--k');                             % plotting baseline Bayes factor
    ignore.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on

    alt = BF_alt.(conditions(1, c));                                        % this condition's indicies for alternative support
    nul = BF_nul.(conditions(1, c));                                        % this condition's indicies for null support
    inc = BF_inc.(conditions(1, c));                                        % this condition's indicies for inclusive support
    scatter(alt(:, 1), alt(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5], ... % plotting scatter plots
                                      'MarkerFaceColor', colours_BF(c, :))
    scatter(inc(:, 1), inc(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5])
    scatter(nul(:, 1), nul(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
                                      'MarkerFaceColor', [0.5, 0.5, 0.5])
    xticks(20:20:time_points)                                               % x-axis ticks
    xticklabels(-100:100:700)                                               % x-axis labels
    xlim([0, time_points-1])                                                % x-axis limits
    set(gca, 'YScale', 'log')                                               % changing y-axis scale to log
    ylim([10^(-BF_scatter_lim),10^BF_scatter_lim])                          % y-axis limits
    yticks(10^(-BF_scatter_lim):10^BF_scatter_lim:10^BF_scatter_lim)

    if c == 2
        ylabel('BF', 'FontSize', 14)
    end
    set(findobj(gcf,'type','axes'), 'FontSize', 12, 'box', 'off')
    vline(onset, '--k');
    lgd = legend('BF > 6', '1/6 < BF < 6', '1/6 > BF');
    legend boxoff
    set(lgd, 'Position', lgd_pos(c, :))
    title(lgd, legend_titles(c, 1), 'FontWeight', 'normal')
end
xlabel('Time (ms)', 'FontSize', 14)                                         % x-axis title
