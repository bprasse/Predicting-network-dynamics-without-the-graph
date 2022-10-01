function Fig_4( seed, N, num_observations, network_type )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%N is the number of nodes for random graph models
if nargin < 2
    N = 300;
end

%set the number of observations n
if nargin < 3
    num_observations = 1e2;
end

%the network type:
%'real_world' for the empirical networks in the main text of the manuscript
%'ER' are Erdos-Renyi random graphs
%'BA' are Barabasi-Albert random graphs
if nargin < 4
    network_type = 'ER'; %'real_world';'BA';
end

%the filename to save/read the results
filename = strcat( './results/Fig_3_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', network_type );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) ~= 2
    Fig_3_S5_S6_S7_S8
end

%load the data
load( filename )

%the model abbreviations
model_names = { 'LV', 'MP', 'MM', 'SIS', 'KUR', 'CW' };

%open a new large figure 
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over every model
for model_counter = 1:6
    
    subplot( 3, 2, model_counter )
    
    %get the model-specific true network B and surrogate network B_est
    switch model_counter
        case 1
            B = results_LV.B;
            B_est = results_LV.B_hat;
        case 2
            B = results_MP.B;
            B_est = results_MP.B_hat;
        case 3
            B = results_MM.B;
            B_est = results_MM.B_hat;
        case 4
            B = results_SIS.B;
            B_est = results_SIS.B_hat;
        case 5
            B = results_kuramoto.B;
            B_est = results_kuramoto.B_hat;
        case 6
            B = results_cw.B;
            B_est = results_cw.B_hat_cw_hat;
    end
    
    %the number of nodes
    N = size( B, 1 );
    
    %compute the degree distribution of the true networks B and B_est 
    in_degrees = transpose( sum( B > 0, 2 ));
    in_degrees_est = transpose( sum( B_est > 0, 2 ));
    
    in_degree_distribution = zeros( 1, N ); %degree=0, 1, ..., N-1
    in_degree_distribution_est = zeros( 1, N );
    
    for degree = 1:N
        in_degree_distribution( degree ) = sum( in_degrees < degree );
        in_degree_distribution_est( degree ) = sum( in_degrees_est < degree );
    end
    
    in_degree_distribution = in_degree_distribution/N;    
    in_degree_distribution_est = in_degree_distribution_est/N;
    
    %plot the true and estimated degree distribution    
    loglog( 1- in_degree_distribution, '-xb' )
    hold on
    loglog( 1- in_degree_distribution_est, '-or' )
    xlabel( 'Degree d' )
    ylabel( 'Distribution Pr[D>=d]')
    legend( 'In-degree true', 'In-degree est', 'Location', 'best')
    title( model_names{ model_counter } )
    
end
end


