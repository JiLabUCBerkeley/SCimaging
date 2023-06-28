% Make pairwise comparisons of fitted angle and preferred direction across
% sessions, for OS ROIs

% Create a folder within the stack folder, e.g. "session_compare"
% Copy the allPar files into this folder and rename according to type of
% stimulus, e.g. grating, ortho drifting bars, etc.
close all

filePath = 'E:\data_summaries\yajies_project\combined_data\all_sessions';
identifier = 'all_combined'; 

filename_grating = 'grating';
filename_counterclk = 'counterclk';
filename_clockws = 'clockws';
filename_ortho = 'ortho';

result = 'test';
% Create directory to save stuff
if ~exist(fullfile(filePath,result))
    mkdir(fullfile(filePath,result));
end

% Load the data for each file
% Drifting grating:
file_grating = fullfile(filePath,filename_grating);
load(file_grating);
data_grating = data;

% Orthogonal bars:
file_ortho = fullfile(filePath,filename_ortho);
load(file_ortho);
data_ortho = data;

% Clockwise bars:
file_clockws = fullfile(filePath,filename_clockws);
load(file_clockws);
data_clockws = data;

% Counterclockwise bars:
file_counterclk = fullfile(filePath,filename_counterclk);
load(file_counterclk);
data_counterclk = data;



% Get data for OS, preferred direction, and fitted angle
load(file_grating, 'xlsHead'); %headings for the spreadsheet
% Find indices for the relevant columns
for n = 1:length(xlsHead)
    if strcmp(xlsHead{n},'isOS?'), ind_os = n; end
    if strcmp(xlsHead{n},'theta(pref_direction)'), ind_prefDir = n; end
    if strcmp(xlsHead{n},'fittedAngle'), ind_fitAng = n; end
end  

% Use these indices to pull the data from each dataset
grating_os = data_grating(:,ind_os);
grating_prefDir = data_grating(:,ind_prefDir);
grating_fitAng = data_grating(:,ind_fitAng);
grating_fitAng(grating_fitAng>=180) = grating_fitAng(grating_fitAng>=180)-180; %scale fitAng to 0-180

% Orthogonal bars:
ortho_os = data_ortho(:,ind_os); 
ortho_prefDir = data_ortho(:,ind_prefDir); 
ortho_fitAng = data_ortho(:,ind_fitAng); 
ortho_fitAng(ortho_fitAng>=180) = ortho_fitAng(ortho_fitAng>=180)-180; %scale fitAng to 0-180

% Clockwise bars:
clockws_os = data_clockws(:,ind_os); 
clockws_prefDir = data_clockws(:,ind_prefDir); 
clockws_fitAng = data_clockws(:,ind_fitAng); 
clockws_fitAng(clockws_fitAng>=180) = clockws_fitAng(clockws_fitAng>=180)-180; %scale fitAng to 0-180

% Counterclockwise bars:
counterclk_os = data_counterclk(:,ind_os); 
counterclk_prefDir = data_counterclk(:,ind_prefDir); 
counterclk_fitAng = data_counterclk(:,ind_fitAng); 
counterclk_fitAng(counterclk_fitAng>=180) = counterclk_fitAng(counterclk_fitAng>=180)-180; %scale fitAng to 0-180

% For each ROI, compare pairwise between the 4 conditions. For each pair,
% if the ROI is OS under both conditions, take the difference in fitAng and
% prefDir and add to a matrix.
fitAng_grating_vs_ortho = [];
fitAng_grating_vs_clockws = [];
fitAng_grating_vs_counterclk = [];
fitAng_ortho_vs_clockws = [];
fitAng_ortho_vs_counterclk = [];
fitAng_clockws_vs_counterclk = [];
prefDir_grating_vs_ortho = [];
prefDir_grating_vs_clockws = [];
prefDir_grating_vs_counterclk = [];
prefDir_ortho_vs_clockws = [];
prefDir_ortho_vs_counterclk = [];
prefDir_clockws_vs_counterclk = [];

% Also, keep track of how many ROIs are tuned under both conditions, for
% each pairwise comparison
grating_ortho_os = [];
grating_clockws_os = [];
grating_counterclk_os = [];
ortho_clockws_os = [];
ortho_counterclk_os = [];
clockws_counterclk_os = [];

nROI = size(data_grating,1);

% Calculate difference in fitAng and prefDir
for ii = 1:nROI
    if grating_os(ii) == 1 && ortho_os(ii) == 1 
        d_fitAng = abs(grating_fitAng(ii) - ortho_fitAng(ii));
        fitAng_grating_vs_ortho = [fitAng_grating_vs_ortho;d_fitAng];
        d_prefDir = abs(grating_prefDir(ii) - ortho_prefDir(ii));
        prefDir_grating_vs_ortho = [prefDir_grating_vs_ortho;d_prefDir];
        grating_ortho_os = [grating_ortho_os;1];
    end
    
    if grating_os(ii) == 1 && clockws_os(ii) == 1 
        d_fitAng = abs(grating_fitAng(ii) - clockws_fitAng(ii));
        fitAng_grating_vs_clockws = [fitAng_grating_vs_clockws;d_fitAng];
        d_prefDir = abs(grating_prefDir(ii) - clockws_prefDir(ii));
        prefDir_grating_vs_clockws = [prefDir_grating_vs_clockws;d_prefDir];
        grating_clockws_os = [grating_clockws_os;1];
    end
    
    if grating_os(ii) == 1 && counterclk_os(ii) == 1 
        d_fitAng = abs(grating_fitAng(ii) - counterclk_fitAng(ii));
        fitAng_grating_vs_counterclk = [fitAng_grating_vs_counterclk;d_fitAng];
        d_prefDir = abs(grating_prefDir(ii) - counterclk_prefDir(ii));
        prefDir_grating_vs_counterclk = [prefDir_grating_vs_counterclk;d_prefDir];
        grating_counterclk_os = [grating_counterclk_os;1];
    end
    
    if ortho_os(ii) == 1 && clockws_os(ii) == 1 
        d_fitAng = abs(ortho_fitAng(ii) - clockws_fitAng(ii));
        fitAng_ortho_vs_clockws = [fitAng_ortho_vs_clockws;d_fitAng];
        d_prefDir = abs(ortho_prefDir(ii) - clockws_prefDir(ii));
        prefDir_ortho_vs_clockws = [prefDir_ortho_vs_clockws;d_prefDir];
        ortho_clockws_os = [ortho_clockws_os;1];
    end
    
    if ortho_os(ii) == 1 && counterclk_os(ii) == 1 
        d_fitAng = abs(ortho_fitAng(ii) - counterclk_fitAng(ii));
        fitAng_ortho_vs_counterclk = [fitAng_ortho_vs_counterclk;d_fitAng];
        d_prefDir = abs(ortho_prefDir(ii) - counterclk_prefDir(ii));
        prefDir_ortho_vs_counterclk = [prefDir_ortho_vs_counterclk;d_prefDir];
        ortho_counterclk_os = [ortho_counterclk_os;1];
    end
    
    if clockws_os(ii) == 1 && counterclk_os(ii) == 1 
        d_fitAng = abs(clockws_fitAng(ii) - counterclk_fitAng(ii));
        fitAng_clockws_vs_counterclk = [fitAng_clockws_vs_counterclk;d_fitAng];
        d_prefDir = abs(clockws_prefDir(ii) - counterclk_prefDir(ii));
        prefDir_clockws_vs_counterclk = [prefDir_clockws_vs_counterclk;d_prefDir];
        clockws_counterclk_os = [clockws_counterclk_os;1];
    end
end

% Put all the fitAng differences on a 0-90 scale, and prefDir differences
% on a 0-180 scale
fitAng_grating_vs_ortho(fitAng_grating_vs_ortho>=90) = fitAng_grating_vs_ortho(fitAng_grating_vs_ortho>=90)-90; 
fitAng_grating_vs_clockws(fitAng_grating_vs_clockws>=90) = fitAng_grating_vs_clockws(fitAng_grating_vs_clockws>=90)-90;
fitAng_grating_vs_counterclk(fitAng_grating_vs_counterclk>=90) = fitAng_grating_vs_counterclk(fitAng_grating_vs_counterclk>=90)-90;
fitAng_ortho_vs_clockws(fitAng_ortho_vs_clockws>=90) = fitAng_ortho_vs_clockws(fitAng_ortho_vs_clockws>=90)-90;
fitAng_ortho_vs_counterclk(fitAng_ortho_vs_counterclk>=90) = fitAng_ortho_vs_counterclk(fitAng_ortho_vs_counterclk>=90)-90;
fitAng_clockws_vs_counterclk(fitAng_clockws_vs_counterclk>=90) = fitAng_clockws_vs_counterclk(fitAng_clockws_vs_counterclk>=90)-90;
prefDir_grating_vs_ortho(prefDir_grating_vs_ortho>=180) = prefDir_grating_vs_ortho(prefDir_grating_vs_ortho>=180)-180; 
prefDir_grating_vs_clockws(prefDir_grating_vs_clockws>=180) = prefDir_grating_vs_clockws(prefDir_grating_vs_clockws>=180)-180;
prefDir_grating_vs_counterclk(prefDir_grating_vs_counterclk>=180) = prefDir_grating_vs_counterclk(prefDir_grating_vs_counterclk>=180)-180;
prefDir_ortho_vs_clockws(prefDir_ortho_vs_clockws>=180) = prefDir_ortho_vs_clockws(prefDir_ortho_vs_clockws>=180)-180;
prefDir_ortho_vs_counterclk(prefDir_ortho_vs_counterclk>=180) = prefDir_ortho_vs_counterclk(prefDir_ortho_vs_counterclk>=180)-180;
prefDir_clockws_vs_counterclk(prefDir_clockws_vs_counterclk>=180) = prefDir_clockws_vs_counterclk(prefDir_clockws_vs_counterclk>=180)-180;


% --------------------------------------------------
% Add up distributions across pairwise comparisons
% (1) Omit grating, just bars:
fitAng_all_bars = cat(1,fitAng_ortho_vs_clockws,fitAng_ortho_vs_counterclk,fitAng_clockws_vs_counterclk);
prefDir_all_bars = cat(1,prefDir_ortho_vs_clockws,prefDir_ortho_vs_counterclk,prefDir_clockws_vs_counterclk);

% (2) Include grating:
fitAng_bars_grating = cat(1,fitAng_grating_vs_ortho,fitAng_grating_vs_clockws,fitAng_grating_vs_counterclk,fitAng_ortho_vs_clockws,fitAng_ortho_vs_counterclk,fitAng_clockws_vs_counterclk);
prefDir_bars_grating = cat(1,prefDir_grating_vs_ortho,prefDir_grating_vs_clockws,prefDir_grating_vs_counterclk,prefDir_ortho_vs_clockws,prefDir_ortho_vs_counterclk,prefDir_clockws_vs_counterclk);

% --------------------------------------------------------------
% Calculate % of ROIs that are OS with grating and OS with bars
nROI_OS_grating = sum(grating_os); % # of ROIs that are OS with grating
nROI_OS_ortho = sum(ortho_os); % # of ROIs that are OS with ortho bars
nROI_OS_clockws = sum(clockws_os); % # of ROIs that are OS with clockwise bars
nROI_OS_counterclk = sum(counterclk_os); % # of ROIs that are OS with counterclockwise bars
nROI_OS_grating_ortho = sum(grating_ortho_os); % # of ROIs that are OS with both grating and ortho bars
nROI_OS_grating_clockws = sum(grating_clockws_os); % # of ROIs that are OS with both grating and clockwise bars
nROI_OS_grating_counterclk = sum(grating_counterclk_os); % # of ROIs that are OS with both grating and counterclockwise bars
nROI_OS_ortho_clockws = sum(ortho_clockws_os); % # of ROIs that are OS with both ortho and clockwise bars
nROI_OS_ortho_counterclk = sum(ortho_counterclk_os); % # of ROIs that are OS with both ortho and counterclockwise bars
nROI_OS_clockws_counterclk = sum(clockws_counterclk_os); % # of ROIs that are OS with both clockwise and counterclockwise bars
fr_grating = [nROI_OS_grating_ortho nROI_OS_grating_clockws nROI_OS_grating_counterclk]/nROI_OS_grating; %fraction of ROIs tuned with grating that are tuned in other conditions
fr_ortho = [nROI_OS_grating_ortho nROI_OS_ortho_clockws nROI_OS_ortho_counterclk]/nROI_OS_ortho; %fraction of ROIs tuned with ortho bars that are tuned in other conditions
fr_clockws = [nROI_OS_grating_clockws nROI_OS_ortho_clockws nROI_OS_clockws_counterclk]/nROI_OS_clockws; %fraction of ROIs tuned with clockwise bars that are tuned in other conditions
fr_counterclk = [nROI_OS_grating_counterclk nROI_OS_ortho_counterclk nROI_OS_clockws_counterclk]/nROI_OS_counterclk; %fraction of ROIs tuned with counterclockwise bars that are tuned in other conditions


% --------------------------------------------------
% Save the data
save(fullfile(filePath, ['data_' identifier '.mat']),'fitAng_grating_vs_ortho','fitAng_grating_vs_clockws','fitAng_grating_vs_counterclk',...
    'fitAng_ortho_vs_clockws','fitAng_ortho_vs_counterclk','fitAng_clockws_vs_counterclk','prefDir_grating_vs_ortho',...
    'prefDir_grating_vs_clockws','prefDir_grating_vs_counterclk','prefDir_ortho_vs_clockws','prefDir_ortho_vs_counterclk',...
    'prefDir_clockws_vs_counterclk','fitAng_all_bars','prefDir_all_bars','fitAng_bars_grating','prefDir_bars_grating',...
    'fr_grating','fr_ortho','fr_clockws','fr_counterclk')




%%
% Make figures
% -------------

% Distribution of change in fitted angle, grating vs orthogonal bars
fig1 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_grating_vs_ortho,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle grating vs ortho bars';
save_figures(filePath,result,identifier,fig1,figtitle)

% Distribution of change in preferred direction, grating vs orthogonal bars
fig2 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_grating_vs_ortho,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir grating vs ortho bars';
save_figures(filePath,result,identifier,fig2,figtitle)

% Distribution of change in fitted angle, grating vs clockwise bars
fig3 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_grating_vs_clockws,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle grating vs clockwise bars';
save_figures(filePath,result,identifier,fig3,figtitle)

% Distribution of change in preferred direction, grating vs clockwise bars
fig4 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_grating_vs_clockws,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir grating vs clockwise bars';
save_figures(filePath,result,identifier,fig4,figtitle)

% Distribution of change in fitted angle, grating vs counterclockwise bars
fig5 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_grating_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle grating vs counterclk bars';
save_figures(filePath,result,identifier,fig5,figtitle)

% Distribution of change in preferred direction, grating vs counterclockwise bars
fig6 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_grating_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir grating vs counterclk bars';
save_figures(filePath,result,identifier,fig6,figtitle)

% Distribution of change in fitted angle, ortho vs clockwise bars
fig5 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_ortho_vs_clockws,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle ortho vs clockwise bars';
save_figures(filePath,result,identifier,fig5,figtitle)

% Distribution of change in preferred direction, ortho vs clockwise bars
fig6 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_ortho_vs_clockws,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir ortho vs clockwise bars';
save_figures(filePath,result,identifier,fig6,figtitle)

% Distribution of change in fitted angle, ortho vs counterclockwise bars
fig7 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_ortho_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle ortho vs counterclk bars';
save_figures(filePath,result,identifier,fig7,figtitle)

% Distribution of change in preferred direction, ortho vs counterclockwise bars
fig8 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_ortho_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir ortho vs counterclk bars';
save_figures(filePath,result,identifier,fig8,figtitle)

% Distribution of change in fitted angle, clockwise vs counterclockwise bars
fig9 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_clockws_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle clockws vs counterclk bars';
save_figures(filePath,result,identifier,fig9,figtitle)

% Distribution of change in preferred direction, clockwise vs counterclockwise bars
fig10 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_clockws_vs_counterclk,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir clockws vs counterclk bars';
save_figures(filePath,result,identifier,fig10,figtitle)

% Distribution of change in fitted angle, all bars together
fig11 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_all_bars,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle all bars';
save_figures(filePath,result,identifier,fig11,figtitle)

% Distribution of change in preferred direction, all bars together
fig12 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_all_bars,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir all bars';
save_figures(filePath,result,identifier,fig12,figtitle)

% Distribution of change in fitted angle, all bars and grating together
fig13 = figure;
hold on
edges = 0:10:90;
histogram(fitAng_bars_grating,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 90]);
xticks([0 30 60 90]);
hold off;
figtitle = 'fitted angle all bars grating';
save_figures(filePath,result,identifier,fig13,figtitle)

% Distribution of change in preferred direction, all bars and grating together
fig14 = figure;
hold on
edges = 0:10:180;
histogram(prefDir_bars_grating,edges);
set(gca,'GridAlpha',0.8,'FontSize',18,'TickDir','out');
xlim([0 180]);
xticks([0 30 60 90 120 150 180]);
hold off;
figtitle = 'pref dir all bars grating';
save_figures(filePath,result,identifier,fig14,figtitle)


% Fractions of ROIs that are tuned under each condition
% ROIs tuned for grating and other conditions
fig15 = figure;
hold on
bar(fr_grating) 
ylim([0 1]);
set(gca,'XTick',[])
hold off
figtitle = 'fraction grating tuned';
save_figures(filePath,result,identifier,fig15,figtitle)

% ROIs tuned for ortho bars and other conditions
fig16 = figure;
hold on
bar(fr_ortho) 
ylim([0 1])
set(gca,'XTick',[])
hold off
figtitle = 'fraction ortho tuned';
save_figures(filePath,result,identifier,fig16,figtitle)

% ROIs tuned for clockwise bars and other conditions
fig17 = figure;
hold on
bar(fr_clockws) 
ylim([0 1])
set(gca,'XTick',[])
hold off
figtitle = 'fraction clockwise tuned';
save_figures(filePath,result,identifier,fig17,figtitle)

% ROIs tuned for counterclockwise bars and other conditions
fig18 = figure;
hold on
bar(fr_counterclk) 
ylim([0 1])
set(gca,'XTick',[])
hold off
figtitle = 'fraction counterclk tuned';
save_figures(filePath,result,identifier,fig18,figtitle)




%%
function save_figures(filePath,result,identifier,fig,figtitle)
savefig(fig,fullfile(filePath,result,[figtitle ' ' identifier]))
saveas(fig,fullfile(filePath,result,[figtitle ' ' identifier]),'png')
print(fig,fullfile(filePath,result,[figtitle ' ' identifier]),'-depsc','-tiff')
end

