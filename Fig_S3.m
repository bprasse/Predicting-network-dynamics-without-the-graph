function Fig_S3( seed )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%the filename to save/read the results
filename = strcat( './results/Fig_S3_seed_',  num2str( seed ) );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    %the model abbreviations
    model_name = { 'LV', 'MP', 'MM', 'SIS', 'KUR', 'WC' };
    
    %the network sizes
    N_all = [100:100:1000];
    
    %the number of different network sizes
    num_N = length( N_all );
    
    %set the number of observations n
    num_observations = 1e3;
    
    %the number of different observation times t_obs
    num_t_samples = 20;
    
    %the different observation times t_obs, evenly spaced from 50 to num_observations
    t_samples = ceil( linspace( 50, num_observations, num_t_samples ) );
    
    %the network type:
    random_network_model = 'BA';
    
    %the number of networks (for each network size in N_all)
    num_networks = 100;
    
    %pre-allocate the rank of the nodal state
    rank_vs_t = nan( num_networks, num_t_samples, num_N, length( model_name ) );
    rank_vs_t_relative_to_N = nan( num_networks, num_t_samples, num_N, length( model_name ) );
    
    %set the parameters
    parameters = set_parameters;
    
    %set the maximum observation time
    parameters.T_max = 10*[ 8  ...% LV
        2  ...% MP
        4  ...% MM
        2 ...% SIS
        2 ...%KUR
        4 ]; %CW
    
    %set the network parameters
    m0_BA = parameters.m0_BA;
    m_BA = parameters.m_BA;
    
    %loop over all network sizes
    for N_i = 1:num_N
        N = N_all( N_i );
        
        %loop over all models
        for model_i = 1:length( model_name )
            
            %loop over networks
            for run_i = 1:num_networks
                
                %the zero-one (unweighted) network
                A = create_network( N, random_network_model, [], [], [], m0_BA, m_BA, true );
                
                %the link weights
                link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
                
                %the weighted network
                B = zeros( N );
                B( A > 0 ) = link_weights;
                
                %the observation time samples
                t_obs = linspace( 0, parameters.T_max( model_i ), num_observations );
                
                switch model_i
                    case 1
                        %Lotka-Volterra
                        
                        %the initial nodal state
                        x_init = parameters.LV.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        alpha = 1 + parameters.LV.sigma_alpha*( 2*rand( N, 1 ) - 1);
                        theta = 1 + parameters.LV.sigma_theta*( 2*rand( N, 1 ) - 1);
                        
                        %generate the nodal state sequence
                        [ x, dx ] = compute_nodal_states_lotka_volterra( x_init, alpha, theta, B, t_obs );
                        
                    case 2
                        %Mutualistic population
                        
                        %the initial nodal state
                        x_init = parameters.MP.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        alpha = 1 + parameters.MP.sigma_alpha*( 2*rand( N, 1 ) - 1 );
                        theta = 1 + parameters.MP.sigma_theta*( 2*rand( N, 1 ) - 1 );
                        
                        %generate the nodal state sequence
                        [ x, dx ] = compute_nodal_states_mutualistic_pop( x_init, alpha, theta, B, t_obs );
                    case 3
                        %Michaelis-Menten
                        
                        %the initial nodal state
                        x_init = parameters.MM.x_init_max*rand( N, 1 );
                        
                        %generate the nodal state sequence
                        [ x, dx ] = compute_nodal_states_michaelis_menten( x_init, parameters.MM.hill_coeff, B, t_obs  );
                        
                    case 4
                        %SIS
                        
                        %the initial nodal state
                        x_init = parameters.SIS.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        delta_init = 1 + parameters.SIS.sigma_delta*( 2*rand( N, 1 ) - 1);
                        R_0_init = eigs( diag( 1./sqrt( delta_init ))*B*diag( 1./sqrt( delta_init )), 1 );
                        delta = R_0_init./parameters.SIS.R_0_SIS*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B
                        
                        %generate the nodal state sequence
                        [ x, dx] = compute_nodal_states_SIS( x_init, delta, B, t_obs );
                    case 5
                        %Kuramoto
                        
                        %the initial nodal state
                        x_init = ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )*rand( N, 1 ) - ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )/2;
                        
                        %the parameters of the model
                        omega = parameters.kuramoto.sigma_omega*randn( N, 1 );
                        
                        %generate the nodal state sequence
                        [ x, dx] = compute_nodal_states_kuramoto( x_init, omega, B, t_obs);
                        
                    case 6
                        %Cowan-Wilson
                        
                        %the initial nodal state
                        x_init = parameters.cw.x_init_max*rand( N, 1 );
                        
                        %generate the nodal state sequence
                        [ x, dx] = compute_nodal_states_cowan_wilson( x_init, parameters.cw.tau, parameters.cw.mu, B, t_obs  );
                end
                
                %compute the absolute and relative rank of the nodal state
                %sequence x for all observation times t_obs
                for t_i = 1:num_t_samples
                    rank_vs_t( run_i, t_i, N_i, model_i ) = rank( x( :, 1:t_samples( t_i )) );
                    rank_vs_t_relative_to_N( run_i, t_i, N_i, model_i ) = rank_vs_t( run_i, t_i, N_i, model_i )/N;
                end
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% save file  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save( strcat( filename, '.mat' ) )
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot results, corresponding to Fig_S3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open new figure
figure

%Only plot the curves for a subset of t_obs (otherwise the plot is too
%busy)
t_obs_ind = [2 4 8 12 16 20 ];

%loop over all models
for model_i = 1:6    
    %compute the average
    data_plot_tmp = squeeze( mean( rank_vs_t_relative_to_N( :, t_obs_ind, :, model_i ), 1 ));
            
    %plot for model_i
    subplot( 3, 2, model_i )
    loglog( N_all, data_plot_tmp )
    hold
    loglog( N_all, 100./N_all, '--xk' )
    ylim([0 1])
    xlabel( 'Number of Nodes N' )
    ylabel( 'Relative Agitation m/N')
    title( model_name{ model_i })
    
    if model_i==1
        legend("t_{obs} = " + string( t_samples( t_obs_ind ) ), 'Location', 'best', 'NumColumns',ceil( length(t_samples( t_obs_ind ))/2))
    end    
end

end

