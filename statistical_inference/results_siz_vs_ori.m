% Code written to generate Fig. 3.
%
% author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

%% Defining important variables
clear 
path = "";                                                                  % path to MVPA directory
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
smooth_win = 2;                                                            % time points being smoothed across within time-series data
smooth_mthd = 'movmean';                                                    % time-series smoothing method     
BF_crit = 6;                                                                % critical Bayes factor
BF_cats = ["alt"; "nul"; "inc"];                                            % Bayes factor categories
sz = 10;                                                                    % scatter plot data point size
BF_scatter_lim = 8;                                                         % y-axis limit of BF scatter plots
ylim_values = [49.4, 52.5];                                                 % y-axis limits

%% Aggregating data
time_series = zeros(N, time_points);                                        % array for storing time-series data

for n = 1:N
    load(data_files(n).name)
    
    siz_vs_ori = squeeze(mean(siz_vs_ori, 1)).*percentage_scaler;           % averaging across runs    
    siz_vs_ori_diag = max(siz_vs_ori.*eye(time_points));                    % extracting diagonal time-series from TGM

    siz_vs_ori_diag = smoothdata(siz_vs_ori_diag, smooth_mthd, smooth_win); % smoothing time-series

    time_series(n, :) = siz_vs_ori_diag;
end

%% Calculating time-series Bayes factors
BF_time_series = bayesfactor(time_series', 'nullvalue', chance, ...
                                           'interval', interval, ...
                                           'rscale', rscale);

%% Categorising time-series Bayes factors
alt = NaN(time_points, 2);                                                  % empty matrix for storing values
nul = NaN(time_points, 2);
inc = NaN(time_points, 2);

for t = 1:time_points                                                       % for loop for categorising Bayes factors
    score = BF_time_series(t, 1);
    if score > BF_crit
        alt(t, :) = [t, score];
    elseif score < 1/BF_crit
        nul(t, :) = [t, score];
    else
        inc(t, :) = [t, score];
    end
end

alt(any(isnan(alt), 2), :) = [];                                            % remove NaNs for each matrix
nul(any(isnan(nul), 2), :) = [];
inc(any(isnan(inc), 2), :) = [];

%% Initialising figure
close all
figure('Position', [10, 10, 600, 400]);

%% Accuracy subplot
subplot_secs = [1, 2, 3, 4, 5, 6, 7, 8, 9, ...
                11, 12, 13, 14, 15, 16, 17, 18, 19, ...
                21, 22, 23, 24, 25, 26, 27, 28, 29];
subplot(4, 10, subplot_secs)
plot_data = squeeze(mean(time_series, 1));                                  % time-series mean
plot_SE = std(time_series)./sqrt(N);                                        % claculating standard error
ignore = plot(ones(1, time_points).*chance, '--k'); hold on                 % theoretical chance
ignore.Annotation.LegendInformation.IconDisplayStyle = 'off';               % ignoring theoretical chance for legend
[hl, hp] = boundedline(1:time_points, plot_data, plot_SE, '-k', ...         % plotting time-series with standard error
                       'alpha'); hold on
hp.Annotation.LegendInformation.IconDisplayStyle = 'off';                   % removing error bars from legend
ylim(ylim_values)                                                           % figure y-axis limits
vline(onset, '--k'); hold on                                                % vertical line denoting stimulus onset
xticks(20:20:time_points)                                                   % x-axis ticks
xticklabels(-100:100:700)                                                   % x-axis labels
xlim([0, time_points-1])                                                    % x-axis limit
ylabel('Classification Accuracy (%)', 'FontSize', 14)                       % y-axis label
box off

%% Bayes factors subplot
subplot_secs = [31, 32, 33, 34, 35, 36, 37, 38, 39];
subplot(4, 10, subplot_secs)
ignore = plot(ones(1, time_points), '--k'); hold on                         % horizontal line denoting BF = 1
ignore.Annotation.LegendInformation.IconDisplayStyle = 'off';
scatter(alt(:, 1), alt(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5], ...   % alternative hypothesis evidence time points
                                  'MarkerFaceColor', [1, 0, 0])
hold on
scatter(inc(:, 1), inc(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5])       % inconclusive evidence time points
hold on
scatter(nul(:, 1), nul(:, 2), sz, 'MarkerEdgeColor', [0.5, 0.5, 0.5], ...   % null hypothesis evidence time points
                                  'MarkerFaceColor', [0.5, 0.5, 0.5])
hold on
set(gca, 'YScale', 'log')
xticks(20:20:time_points)                                                    % x-axis ticks
xticklabels(-100:100:700)                                                   % x-axis labels
xlim([0, time_points-1])                                                    % x-axis limits
set(gca, 'YScale', 'log')                                                   % changing y-axis scale to log
ylim([10^(-BF_scatter_lim),10^BF_scatter_lim])                              % y-axis limits
yticks(10^(-BF_scatter_lim):10^BF_scatter_lim:10^BF_scatter_lim)            % y-axis ticks
vline(onset, '--k')                                                         % vertical line denoting stimulus onset
ylabel('BF', 'FontSize', 14)                                                % y-axis label
xlabel('Time (ms)', 'FontSize', 14)
box off
lgd = legend('BF > 6', '1/6 < BF < 6', '1/6 > BF');
set(lgd, 'Position', [0.9, 0.185, 0, 0])
legend boxoff
set(findobj(gcf,'type','axes'), 'FontSize', 12, 'box', 'off')
