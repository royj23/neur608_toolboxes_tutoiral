%% Getting started

% scrub scrub scrub
clear all
close all
clc

% If you are looking at this script, I presume you have already cloned the
% github repo for this tutorial. 

% We will need to retrieve the toolbox repositories from github onto your 
% local machine. To do this, open a terminal and change the current 
% directory to the location where you want the toolboxes to be cloned.
% For this tutorial, this directory will be defined as:
toolbox_dir = "/Users/Jess/Documents"; % !!! CHANGE ME

% Here is where I cloned the tutorials repo:
tutorial_dir = strcat(toolbox_dir, "/teaching/2024FALL_NEUR608/tutorial/neur608_toolboxes_tutorial"); % !!! CHANGE ME
addpath(genpath(strcat(tutorial_dir))); 

% Cloning BrainSpace and ENIGMA: run this in terminal
% toolbox_dir="/Users/Jess/Documents"
% cd $toolbox_dir
% git clone https://github.com/MICA-MNI/BrainSpace.git
% git clone https://github.com/MICA-MNI/ENIGMA.git

% Let's set up the directories for this tutorial
brainspace_dir = strcat(toolbox_dir, "/BrainSpace"); % where you cloned brainspace
enigma_dir = strcat(toolbox_dir, "/ENIGMA"); % where you cloned enigma

% We will add the repos to our path to have access to their functions
addpath(genpath(strcat(brainspace_dir,"/matlab"))); % brainspace
addpath(genpath(strcat(enigma_dir,"/matlab"))); % enigma

% Bonus: colormaps
% https://github.com/DrosteEffect/BrewerMap
addpath(genpath('/Users/Jess/Documents/BrewerMap'));


%% Let's load some data and look at it!
% As in Hettwer et al.

% Disorders we want to analyze
disorder_names = {'bipolar'; 'adhd'; 'asd'; 'depression'; 'ocd'; 'schizophrenia'};
CT_d = {size(disorder_names,1),1};

% Load summary statistics and cohen's D maps
for load_disorder = 1:size(disorder_names,1)

    % Summary statistics for this disorder
    sum_stats = load_summary_stats(disorder_names{load_disorder});
    
    % Depending on the disorder, some summary tables are stratified by age
    % group. If adult data is available, we load that. Otherwise we load
    % all the available data. Naming slightly differs for some disorders
    % We only look at cortical thickness here
    if isfield(sum_stats,"CortThick_case_vs_controls_adult")
        CT = sum_stats.CortThick_case_vs_controls_adult;
    elseif isfield(sum_stats,"CortThick_case_vs_controls")
        CT = sum_stats.CortThick_case_vs_controls;
    elseif isfield(sum_stats,"CortThick_case_vs_controls_meta_analysis")
        CT = sum_stats.CortThick_case_vs_controls_meta_analysis;
    end

    % Extract cohen's D case-control difference maps, which correct for
    % age, sex, ICV.
    CT_d{load_disorder} = CT.d_icv;

    % Plotting time
    % Prepare the data: go from parcellation to vertexwise
    CT_d_fsa5 = parcel_to_surface(CT_d{load_disorder}, 'aparc_fsa5');
    cbar_label = ['Cohen D for ', disorder_names{load_disorder}];

    % Project the results on the cortical surface
    f = figure,
    plot_cortical(CT_d_fsa5, 'surface_name', 'fsa5', ...
        'color_range', [-0.35 0.35], 'cmap', 'RdBu_r', ...
        'label_text', cbar_label)
end

close all


%% How patterns of cortical thickness change co-vary across disorders ?

% Restructure the cohen D data in as a big matrix: 
% rows are regions and columns are disorders
CT_d_all = horzcat(CT_d{:});

% Correlate all pairs of regions across disorders 
r_disorders = corr(CT_d_all');

% Plot the result
f = figure, 
    imagesc(r_disorders);
    set(gca,'XTick',[],'YTick',[]);
    axis square;
    colormap(flipud(brewermap([],"RdBu"))); caxis([-1 1]); colorbar;

% Thresholded
thresh = 80;
r_disorders_bin_thresh = bsxfun(@gt, r_disorders, prctile(r_disorders,thresh))'; 

% Plot the result
f = figure, 
    imagesc(r_disorders_bin_thresh);
    set(gca,'XTick',[],'YTick',[]);
    axis square;
    colormap(brewermap([],"Oranges")); colorbar;

% plot degree centrality map / cross-disorder covariance hub map
dc_disorders = sum(r_disorders_bin_thresh);

% Prepare
dc_disorders_fsa5 = parcel_to_surface(dc_disorders, 'aparc_fsa5');

f = figure,
    plot_cortical(dc_disorders_fsa5, 'surface_name', 'fsa5', ...
    'color_range', [0 40], 'cmap', 'viridis', ...
    'label_text', 'Degree centrality or hubs of strongest cross-disorder correlations')

close all


%% Macroscale organization of transdiagnostic covariance in cortical thickness alterations

% Building gradients from correlation matrix
gm = GradientMaps('kernel', 'na', 'approach', 'dm', 'n_components', 5); 
gm = gm.fit(r_disorders, 'sparsity', thresh);
gradients = gm.gradients{1};

% Scree plot
f = figure, 
    handles = scree_plot(gm.lambda{1});
    set(handles.axes, 'FontName', 'Gill Sans MT', 'FontSize', 12, ...
        'XLim', [0.9 6], 'XTick', [2 4 6], 'XTickLabel', {'2','4','6'}, ...
        'YLim', [0 0.45], 'YTick', [0 0.2 0.4], 'YTickLabel', {'0','0.2','0.4'});

% Plotting time!
for this_gradient = 1:2

    % Prepare the data: go from parcellation to vertexwise
    CT_psych_gradient = gradients(:,this_gradient);
    climits = [prctile(CT_psych_gradient,2.5), prctile(CT_psych_gradient,97.5)];
    CT_grad_fsa5 = parcel_to_surface(CT_psych_gradient, 'aparc_fsa5');
    cbar_label = char(strcat("Gradient ", string(this_gradient)));

    % Project the results on the cortical surface
    f = figure,
    plot_cortical(CT_grad_fsa5, 'surface_name', 'fsa5', ...
        'color_range', climits, 'cmap', 'viridis', ...
        'label_text', cbar_label)

end

% Contextualize with cytoarchitecture
f = figure,
class_mean_G1 = economo_koskinas_spider(gradients(:,1), 'axis_range', [-0.3 0.3]);

f = figure,
class_mean_G2 = economo_koskinas_spider(gradients(:,2), 'axis_range', [-0.3 0.3]);

close all


