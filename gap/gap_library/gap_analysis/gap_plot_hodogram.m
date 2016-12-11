function [figH] = gap_plot_hodogram(times, vector_ts, time_range, figure_window_title ...
    , vars, normal, max_var)
%Plots a hodogram of a vector time series in its frame of reference.
%USAGE:
%   [figH] = gap_plot_hodogram(times, vector_ts, time_range, figure_window_title ...
%                               , vars, normal, max_var)
%WHERE
%   times - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   vector_ts - is a Nx3 array with N row-vectors (e.g. B_mva)
%   time_range - (optional) 2x1 or 1x2 vector with two elements: start and 
%           end of the time span. Both should be given in the datenum 
%           format. Only this time interval is plotted.
%   figure_window_title - (optional) label of the plot window. It also
%           appears on the bottom of the figure.
%   vars - (optional) 3x1 representing eigenvalues of the MVA problem.
%           vars(1) is the smallest. This appears only as text on the plot
%           without any further change.
%   normal - (optional) appears only as a text on the plot.
%   max_var - (optional) appears only as a text on the plot.
%
%OUTPUT
%   figH   - handle to the produced figure.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
if (nargin > 2)
    [times, vector_ts] = gap_select_time_range(times, time_range, vector_ts);
end
if (size(times,1) < 3)
    error('MATLAB:gap_plot_hodogram', 'Too few points for a hodogram');
end

maxes = max(vector_ts,[],1);
mins = min(vector_ts,[],1);
range = max(maxes - mins);
center = (maxes + mins) / 2;
plot_margin = 0.13;

if (nargin > 3)
    figH = figure('name',figure_window_title);
else
    figH = figure;
end
if (nargin <= 4)
    vars = ones(1,3);
end
str_normal = '';
if (nargin > 5)
    str_normal = sprintf('[%.2f,%.2f,%.2f]', normal(1), normal(2), normal(3));
end
str_max_var = '';
if (nargin > 6)
    str_max_var = sprintf('[%.2f,%.2f,%.2f]', max_var(1), max_var(2), max_var(3));
end

subplot(6,2,[1,3,5]);  %max vs intermediate
subplot('Position',[plot_margin 0.5 0.5-plot_margin 1-plot_margin-0.5]);
plot(vector_ts(:,2), vector_ts(:,3), vector_ts(1,2), vector_ts(1,3), 'ro');
set(gca,'XAxisLocation','top');
ylabel('z'); 
ylim([center(3)-range/2, center(3)+range/2]);
xlim([center(2)-range/2, center(2)+range/2]);
ylab_str = sprintf('y  v_{max}/v_{int}=%.2f', vars(3)/vars(2));
xlabel(ylab_str); 

subplot(6,2,[2,4,6]);  %max vs min
subplot('Position',[0.5 0.5 1.0-plot_margin-0.5 1-plot_margin-0.5]);
plot(vector_ts(:,1), vector_ts(:,3), vector_ts(1,1), vector_ts(1,3), 'ro');
set(gca,'YAxisLocation','right','XAxisLocation','top');
ylabel(['z max var: ', str_max_var]); 
ylim([center(3)-range/2, center(3)+range/2]);
xlabel('x'); xlim([center(1)-range/2, center(1)+range/2]);

subplot(6,2,[7 9 11]);  %min vs intermediate
subplot('Position',[plot_margin plot_margin 0.5-plot_margin 0.5-plot_margin]);
plot(vector_ts(:,2), vector_ts(:,1), vector_ts(1,2), vector_ts(1,1), 'ro');
ylabel('x'); ylim([center(1)-range/2, center(1)+range/2]);
xlim([center(2)-range/2, center(2)+range/2]);
ylab_str = sprintf('y  v_{int}/v_{min}=%.2f\n%s', vars(2)/vars(1), str_normal);
xlabel(ylab_str); 

height_one=(0.5-plot_margin)/3;
widht_one=1.0-plot_margin-0.5;
for i=1:3
    subplot(6,2,6+i*2);
    subplot('Position',[0.5 plot_margin+(3-i)*height_one widht_one height_one]);
    plot(times, vector_ts(:,i));
    xlim([times(1) times(end)]);
    ylim([center(i)-range/2, center(i)+range/2]);
    datetick('x','MM:SS', 'keeplimits');
    set(gca,'YAxisLocation','right');
    if i<3
        set(gca,'xtick',[]);
    else
        set(gca,'XMinorTick','on');
    end
end
if (nargin > 3)
    xlabel(figure_window_title);
end

end
