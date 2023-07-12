% Code written for generating Fig. 2--Supplement 3. Requires the user to
% have already run get_erp_data.py.
% 
% author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

%% Defining important variables
clear
load chan64
root = '';                                                                  % root of BIDS directory.
cd(root)
erp_file = dir(fullfile(root,'derivatives','results','ERPs','*.mat'));
erp_data = load(fullfile(erp_file.folder,erp_file.name));
time_points = 181;
stim_onset = 40;
conditions = {'siz_vio','ori_vio','control'};
colours = {'-r','-b','-k'};

%% Plotting data
close all
channels.Left_Frontal = {'FP1','AF7','AF3','F7','F5','F3'};
channels.Central_Frontal = {'FPz','AFz','F1','Fz','F2'};
channels.Right_Frontal = {'FP2','AF4','AF8','F4','F6','F8'};
channels.Left_Central = {'FC3','FC5','FT7','C3','C5','T7','CP3','CP5','TP7'};
channels.Central_Central = {'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'};
channels.Right_Central = {'FC4','FC6','FT8','C4','C6','T8','CP4','CP6','TP8'};
channels.Left_Posterior = {'P3','P5','P7','P9','PO7'};
channels.Central_Posterior = {'P1','Pz','P2','PO3','POz','PO4','O1','Oz','O2','Iz'};
channels.Right_Posterior = {'P4','P6','P8','P10','PO8'};

rois = fieldnames(channels);
n_plots = length(rois);
figure('Position',[10,10,800,800])
tiledlayout(sqrt(n_plots),sqrt(n_plots))

for i = 1:n_plots

    nexttile
    plot_channels = channels.(rois{i,1});
    plot_channel_inds = [];

    for j = 1:length(chans)                                                 % getting channel indicies

        check_channel = chans(j).labels;

        if sum(strcmp(plot_channels,check_channel)) > 0

            plot_channel_inds = [plot_channel_inds,j];

        end
    end

    ignore = plot(zeros(1,time_points),'--k'); hold on                                
    ignore.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    hold on

    for c = 1:length(conditions)

        cond_data = erp_data.(conditions{1,c});
        cond_data = squeeze(mean(cond_data(:,plot_channel_inds,:),2));
        cond_mean = squeeze(mean(cond_data,1));
        cond_se = squeeze(std(cond_data,1))/sqrt(size(cond_data,1));

        [hl,hp] = boundedline(1:time_points,cond_mean,cond_se,colours{1,c}, ...
                      'alpha');
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
    end

    xticks(20:20:time_points)                                               % x-axis ticks
    xlim([0,time_points-1])

    if i > 6
        xticklabels(-100:100:700)                                           % x-axis labels
    else
        xticklabels([])
    end

    vline(stim_onset,'--k'); 

    if i == 1 || i == 4 || i == 7
        ylabel('\muV','FontSize',14,'rotation',0)
    end

    this_title = rois{i,1};
    this_title = split(this_title,'_')';
    this_title = [this_title{1,1},' ',this_title{1,2}];

    title(this_title,'FontSize',10)
end

lgd = legend('Size Violation','Orientation Violation','Control');
lgd.FontSize = 12;                                                          % legend font size
legend boxoff 
