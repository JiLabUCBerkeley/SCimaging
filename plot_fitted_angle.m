function plot_fitted_angle(filePath,filename)


close all %close any open figures


if nargin == 0 
    filePath = 'E:\data_summaries\yajies_project\combined_data\clockws';
    filename = 'clockws';
end

% Create directory to save stuff
result = 'test';
if ~exist(fullfile(filePath,result))
    mkdir(fullfile(filePath,result));
end

% Load the data spreadsheet 
file = fullfile(filePath,filename);
load(file);

% Get the OS ROIs and the fitted angles
for n = 1:length(xlsHead)
    if strcmp(xlsHead{n},'isOS?'), os_rois = data(:,n); end
    if strcmp(xlsHead{n},'fittedAngle'), ind_fit = n; end
end    



% Make the 0's into NaN to ignore the rois that aren't OS, DS, etc
% (otherwise get a bunch of 0s in the plot that don't mean anything)
os_rois(os_rois==0) = NaN;

% Fitted angle for OS ROIs
data_fitAng = data(:,ind_fit);
data_fitAng_os_rois = data_fitAng.*os_rois;

% Put all values on 0-180 scale
data_fitAng_os_rois(data_fitAng_os_rois>=180) = data_fitAng_os_rois(data_fitAng_os_rois>=180)-180;





% Plot fitted angle distribution

fig = figure;
hold on;
edges = 0:10:180;
h = histogram(data_fitAng_os_rois,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 45 90 135 180]);
hold off;
      
figtitle = 'Fitted Angle';
savefig(fig,fullfile(filePath,result,figtitle))
saveas(fig,fullfile(filePath,result,figtitle),'png')
saveas(fig,fullfile(filePath,result,figtitle),'pdf')


