function Fig_S13( seed, N, num_observations, random_network_model )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%N is the number of nodes for random graph models
if nargin < 2
    N=100;
end

%set the number of observations n
if nargin < 3
    num_observations = 2e2;
end

%the network type:
if nargin < 4
    random_network_model = 'BA';
end

%the filename to save/read the results
filename = strcat( './results/Fig_S13_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', random_network_model );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    %set the simulation parameters
    parameters = set_parameters;
    
    %set the maximum prediction time
    T_max = 2;
    
    %set the maximum initial nodal state
    x_init_max = 1;
    
    %the number of coefficients of the polynomial
    g_order = 10; % polynomial 1 + x + y + xy + x^2 + y^2 + x^3 +  x^2y +  y^2x + y^3
    
    %the number of networks
    num_networks = 1e2;
    
    %the observation time fraction t_obs/T_max
    T_obs_fraction = 0.25;
    
    %the number of values for the standard deviation sigma_g, that controls
    %the heterogeneity of the coupling functions g_i
    num_g_var = 6;
    
    %the minimum and maximum value for sigma_g
    g_std_min = 1e-3;
    g_std_max = 0.5;
    
    %generate the range of values for sigma_g
    g_std = [0 logspace( log10( g_std_min ), log10( g_std_max ), num_g_var-1 )];
    
    %pre-allocate the prediction errors
    error_pred = nan( num_g_var, num_networks  );
    
    %loop over values for the standard deviation sigma_g
    for g_std_count = 1:num_g_var
        g_std_i = g_std( g_std_count );
        
        %loop over all networks
        for network_i = 1:num_networks
            
            %the zero-one (unweighted) network
            A = create_network( N, random_network_model, [], [], [],  parameters.m0_BA,  parameters.m_BA, true );
            
            %the link weights
            link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
            
            %the weighted network
            B = zeros( N );
            B( A > 0 ) = link_weights;
            
            %all time samples
            t_all = linspace( 0, T_max, num_observations + 1 );
            
            %the sampling time
            delta_T = t_all( 2 ) - t_all( 1 );
            
            %the observation time samples
            t_obs = t_all( 1:ceil( T_obs_fraction*( num_observations + 1) ));
            
            %the prediction time samples
            t_prediction = t_all( ceil( T_obs_fraction*( num_observations + 1) )+1:end );
            
            
            %the coefficients for the self-interaction function f_i
            f_coeff = transpose( [ 0 -1 0 0]);
            
            %the mean values of the coefficients of the coupling funtion g_i,
            %across all nodes i
            g_coeff_mean = transpose( randn( 1, g_order ));
            
            %the coefficients of the coupling funtion g_i, for all nodes i
            g_coeff = g_coeff_mean*ones( 1, N ) +  g_std_i*transpose( randn( N, g_order ));
            g_coeff = g_coeff/100;
            
            %the initial nodal state
            x_init = x_init_max*rand( N, 1 );
            
            %generate the past nodal state sequence
            [ x, ~ ] = compute_nodal_states_g_het( x_init, f_coeff, g_coeff, B, t_obs );
            
            %generate the future nodal state sequence
            [ x_future, ~ ] = compute_nodal_states_g_het( x( :, end ), f_coeff, g_coeff, B, t_prediction);
            
            %the input parameters for the network
            %reconstruction algorithm
            input_param = [];
            input_param.f_coeff = f_coeff;
            input_param.g_coeff = g_coeff;
            
            %numerical differentiation
            dx = transpose( diff( transpose( x ) ) )/delta_T;
            
            %delete last nodal state such that x and dx have the same size
            x = x( :, 1:end-1 );
            
            %start the network reconstruction
            [ B_hat, ~ ] = network_reconstruction( x, dx, 'g_het', input_param, parameters );
            
            %predict the nodal state with the surrogate network B_hat
            x_future_hat = compute_nodal_states_g_het( x_future( :, 1 ), f_coeff, g_coeff, B_hat, t_prediction  );
            
            %save the prediction error
            error_pred( g_std_count, network_i ) = mean( mean( abs( x_future_hat - x_future ) ) );
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% save files  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save( strcat( filename, '.mat' ) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% plots  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fig S13: Prediction error versus heterogeneity of coupling functions g_i
figure
boxplot( transpose( error_pred ),'Labels', cellfun( @(x) num2str(x,'%.2e'), num2cell( g_std ) , 'UniformOutput', false ))
xlabel( 'Standard Deviation $$\sigma_g$$', 'Interpreter', 'Latex')
ylabel( 'Prediction Error $$\bar{\epsilon}$$', 'Interpreter', 'Latex')
end

