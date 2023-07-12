% Code written to generate Fig. 2--Supplement 2 and Fig. 3--Supplement 2.
%
% author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

%% Defining important variables
clear
load chan64
root = '';                                                                  % root of BIDS directory
cd(root)
results_dir = fullfile(root,'derivatives','results','MVPA');
data_files = dir(fullfile(results_dir,'sub-*.mat'));
N = length(data_files);
time_points = 181;
time_bin = 10;
stim_onset = 40;
sig_onset = 70;
sig_offset = 140;
sig_onset_diff = 90;
sig_offset_diff = 120;

%% Aggregating data
clear weights
weights.siz_exp = nan(N,time_points,length(chans));
weights.ori_exp = nan(N,time_points,length(chans));
weights.diff = nan(N,time_points,length(chans)); 

for n = 1:N

    data_file = fullfile(data_files(n).folder,data_files(n).name);
    load(data_file)
    
    siz_exp = squeeze(mean(vio_vs_cont_SL(:,1,:,:),1))*-1;
    ori_exp = squeeze(mean(vio_vs_cont_SL(:,2,:,:),1))*-1;
    diff = squeeze(mean(siz_vs_ori_SL,1))*-1;

    weights.siz_exp(n,:,:) = siz_exp;
    weights.ori_exp(n,:,:) = ori_exp;
    weights.diff(n,:,:) = diff;

end

%% Plotting weight projections for size and orientation expectation (Fig. 2--Supplement 2).
close all
figure('Position',[10,10,400,720]);
comparisons = {'siz_exp';'ori_exp'};
n_tiles = (sig_offset-sig_onset)/time_bin;
t = tiledlayout(n_tiles,length(comparisons));
t.TileSpacing = 'none';
t.Padding = 'compact';
comparison_names = {
    ['Size',newline,'Expectation'], ...
    ['Orientation',newline,'Expectation']
    };

y = [.895,.774,.654,.53,.4075,.285,.164];
ind = 1;

for t = sig_onset:time_bin:sig_offset-time_bin
    
    t_data = nan(length(chans),length(comparisons));

    for c = 1:length(comparisons)
    
        nexttile

        data = squeeze(mean(weights.(comparisons{c,1}),1));
        data = squeeze(mean(data(t:t+time_bin-1,:),1));

        topoplot(data,chans,'style','map','electrodes','off');
        colormap(colorcet('D1A'))

        if t == sig_onset

            this_title = comparison_names{1,c};
            title(this_title,'FontSize',14)

        end

        t_data(:,c) = data;

    end
    
    % data = t_data(:,1) - t_data(:,2);

    % topoplot(data,chans,'style','map','electrodes','off');
    % colormap(colorcet('D1A'))
    
    R = corrcoef(t_data);
    disp(R)
    time_title = [num2str(t*5-200),' to ',num2str(t*5-200+time_bin*5),' ms'];
    annotation('textbox',[0.425,y(ind),0.5,0],'string',time_title, ...
        'FontSize',10,'EdgeColor','none')
    ind = ind+1;
end

cbh = colorbar;
cbh.Layout.Tile = 'south'; 
cbh.Label.String = 'Reconstructed Weight Projection (a.u.)';
cbh.Label.FontSize = 14;
cbh.Ticks = [];

set(findobj(gcf,'type','axes'), 'FontSize', 12)

%% Plotting weight projections for size vs. orientation violation (Fig. 3--Supplement 2).
% figure('Position',[10,10,400,800]);
% n_tiles = (sig_offset_diff-sig_onset_diff)/time_bin;
% t = tiledlayout(n_tiles,1);
% t.TileSpacing = 'none';
% t.Padding = 'compact';
% 
% y = [.825,.535,.24];
% ind = 1;
% 
% for t = sig_onset_diff:time_bin:sig_offset_diff-time_bin
% 
%     nexttile
%     data = squeeze(mean(weights.diff,1));
%     data = squeeze(mean(data(t:t+time_bin-1,:),1));
% 
%     topoplot(data,chans,'style','map','electrodes','off');
%     colormap(colorcet('D1A'))
% 
%     time_title = [num2str(t*5-200),' to ',num2str(t*5-200+time_bin*5),' ms'];
%     annotation('textbox',[0.75,y(ind),0.5,0],'string',time_title, ...
%         'FontSize',12,'EdgeColor','none')
%     ind = ind+1;
% 
%     if t == sig_onset_diff
% 
%         this_title = '  Size Vio. vs. Ori. Vio.';
%         title(this_title,'FontSize',14)
% 
%     end
% 
% end
% 
% cbh = colorbar;
% cbh.Layout.Tile = 'south'; 
% cbh.Label.String = 'Reconstructed Weight Projection (a.u.)';%['Reconstructed Weight',newline,'Projection (a.u.)'];
% cbh.Label.FontSize = 14;
% cbh.Ticks = [];
% 
% set(findobj(gcf,'type','axes'), 'FontSize', 12)
