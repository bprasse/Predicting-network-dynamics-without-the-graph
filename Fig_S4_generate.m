function Fig_S4_generate( seed )
%generate and save the data that is required for Fig_S4
%this function is computationally heavy, ideally run on a cluster with
%different seeds on different nodes

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%the filename to save/read the results
filename = strcat( './results/Fig_S4_seed_', num2str( seed ));

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) ~= 2
    
    %set the number of observations n
    num_observations = 1e3;
    
    %the network sizes
    N_all = fliplr( [5:5:15 20:5:50 75 100:100:500 750 1000] );
    
    %the number of different networks
    num_networks = 10;
    
    %set the parameters
    parameters = set_parameters;
    
    %the maximum observation time for all models
    T_obs = parameters.T_max;
    
    %pre-allocate the rank of the linear system
    rank_vs_N_LV = nan( num_networks, length( N_all ));
    rank_vs_N_MP = nan( num_networks, length( N_all ));
    rank_vs_N_MM = nan( num_networks, length( N_all ));
    rank_vs_N_SIS = nan( num_networks, length( N_all ));
    rank_vs_N_kuramoto = nan( num_networks, length( N_all ));
    rank_vs_N_cw = nan( num_networks, length( N_all ));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %redefine variables for parfor loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma_alpha_LV = parameters.LV.sigma_alpha;
    sigma_theta_LV = parameters.LV.sigma_theta;
    
    sigma_alpha_MP = parameters.MP.sigma_alpha;
    sigma_theta_MP = parameters.MP.sigma_theta;
    
    sigma_delta_SIS = parameters.SIS.sigma_delta;
    R_0_SIS = parameters.SIS.R_0_SIS;
    
    sigma_omega_kuramoto = parameters.kuramoto.sigma_omega;
    
    hill_coeff = parameters.MM.hill_coeff;
    
    tau_cw = parameters.cw.tau;
    mu_cw = parameters.cw.mu;
    
    random_network_model = parameters.random_network_model;
    m0_BA = parameters.m0_BA;
    m_BA = parameters.m_BA;
    
    min_link_weight = parameters.min_link_weight;
    max_link_weight = parameters.max_link_weight;
    
    x_init_max_LV = parameters.LV.x_init_max;
    x_init_max_MP = parameters.MP.x_init_max;
    x_init_max_MM = parameters.MM.x_init_max;
    x_init_max_SIS = parameters.SIS.x_init_max;
    x_init_kuramoto_max = parameters.kuramoto.x_init_max; x_init_kuramoto_min = parameters.kuramoto.x_init_min;
    x_init_max_CW = parameters.cw.x_init_max;
    
    T_obs_LV = T_obs( 1 );
    T_obs_MP = T_obs( 2 );
    T_obs_MM = T_obs( 3 );
    T_obs_SIS = T_obs( 4 );
    T_obs_kuramoto = T_obs( 5 );
    T_obs_CW = T_obs( 6 );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop over all networks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
    parfor curr_network = 1:num_networks
        
        %redefine some variables for parfor loop
        rank_vs_N_LV_this_parfor = nan( 1, length( N_all ) );
        rank_vs_N_MP_this_parfor = nan( 1, length( N_all ) );
        rank_vs_N_MM_this_parfor = nan( 1, length( N_all ) );
        rank_vs_N_SIS_this_parfor = nan( 1, length( N_all ) );
        rank_vs_N_kuramoto_this_parfor = nan( 1, length( N_all ) );
        rank_vs_N_cw_this_parfor = nan( 1, length( N_all ) );
        
        %loop over all network sizes
        for N = N_all
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Generate the network (the same for all models)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %the unweighted network
            A = create_network( N, random_network_model, [], [], [], m0_BA, m_BA, true );
            
            %the link weights
            link_weights = min_link_weight + ( max_link_weight - min_link_weight )*rand( nnz( A ), 1);
            
            %the weighted network
            B = zeros( N );
            B( A>0 ) = link_weights;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Model parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Lotka-Volterra
            alpha_LV = 1 + sigma_alpha_LV*( 2*rand( N, 1 ) - 1);
            theta_LV = 1 + sigma_theta_LV*( 2*rand( N, 1 ) - 1);
            
            %Mutualistic population
            alpha_MP = 1 + sigma_alpha_MP*( 2*rand( N, 1 ) - 1 );
            theta_MP = 1 + sigma_theta_MP*( 2*rand( N, 1 ) - 1 );
            
            %Michaelis-Menten: no parameters
            
            %SIS
            delta_init = 1 + sigma_delta_SIS*( 2*rand( N, 1 ) - 1);
            R_0_init = eigs( diag( 1./sqrt( delta_init ))*B*diag( 1./sqrt( delta_init )), 1 );
            multiplicity_delta = R_0_init./R_0_SIS;
            delta_SIS = multiplicity_delta*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B
            
            %Kuramoto
            omega_kuramoto = sigma_omega_kuramoto*randn( N, 1 );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% the initial nodal state
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            x_init_LV = x_init_max_LV*rand( N, 1 );
            x_init_MP = x_init_max_MP*rand( N, 1 );
            x_init_MM = x_init_max_MM*rand( N, 1 );
            x_init_SIS = x_init_max_SIS*rand( N, 1 );
            x_init_kuramoto = ( x_init_kuramoto_max - x_init_kuramoto_min )*rand( N, 1 ) - ( x_init_kuramoto_max - x_init_kuramoto_min )/2;
            x_init_cw = x_init_max_CW*rand( N, 1 );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% generate the nodal state sequences
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Lotka-Volterra
            t_obs_LV = linspace( 0, T_obs_LV, num_observations );
            [ x_LV, dx_LV ] = compute_nodal_states_lotka_volterra( x_init_LV, alpha_LV, theta_LV, B, t_obs_LV );
            
            %Mutualistic population
            t_obs_MP = linspace( 0, T_obs_MP, num_observations  );
            [ x_MP, dx_MP ] = compute_nodal_states_mutualistic_pop( x_init_MP, alpha_MP, theta_MP, B, t_obs_MP );
            
            %Michaelis-Menten
            t_obs_MM = linspace( 0, T_obs_MM, num_observations );
            [ x_MM, dx_MM ] = compute_nodal_states_michaelis_menten( x_init_MM, hill_coeff, B, t_obs_MM  );
            
            %SIS
            t_obs_SIS = linspace( 0, T_obs_SIS, num_observations );
            [ x_SIS, dx_SIS ] = compute_nodal_states_SIS( x_init_SIS, delta_SIS, B, t_obs_SIS );
            
            %Kuramoto
            t_obs_kuramoto = linspace( 0, T_obs_kuramoto, num_observations );
            [ x_kuramoto, dx_kuramoto ] = compute_nodal_states_kuramoto( x_init_kuramoto, omega_kuramoto, B, t_obs_kuramoto );
            
            %Cowan-Wilson
            t_obs_cw = linspace( 0, T_obs_CW, num_observations );
            [ x_cw, dx_cw  ] = compute_nodal_states_cowan_wilson( x_init_cw, tau_cw, mu_cw, B, t_obs_cw  );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute the rank of the linear systems
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Lotka-Volterra
            rank_vs_N_LV_tmp = 0;
            
            for node = 1:N
                [ ~, F_LV_i ] = linear_system_lotka_volterra( x_LV, dx_LV, alpha_LV( node ), theta_LV( node ), node );
                
                rank_vs_N_LV_tmp =  rank_vs_N_LV_tmp +  rank( F_LV_i );
            end
            rank_vs_N_LV_this_parfor( N_all == N ) =  rank_vs_N_LV_tmp/N;
            
            %Mutualistic population
            rank_vs_N_MP_tmp = 0;
            
            for node = 1:N
                
                [ ~, F_MP_i ] = linear_system_mutualistic_pop( x_MP, dx_MP, alpha_LV( node ), theta_LV( node ), node );
                
                rank_vs_N_MP_tmp =  rank_vs_N_MP_tmp +  rank( F_MP_i );
            end
            rank_vs_N_MP_this_parfor( N_all == N ) =  rank_vs_N_MP_tmp/N;
            
            %Michaelis-Menten
            rank_vs_N_MM_tmp = 0;
            
            for node = 1:N
                
                [ ~, F_MM_i ] = linear_system_michaelis_menten( x_MM, dx_MM, hill_coeff, node );
                
                rank_vs_N_MM_tmp =  rank_vs_N_MM_tmp +  rank( F_MM_i( :, 2:end) );
            end
            rank_vs_N_MM_this_parfor( N_all == N ) =  rank_vs_N_MM_tmp/N;
            
            %SIS
            rank_vs_N_SIS_tmp = 0;
            
            for node = 1:N
                
                [ ~, F_SIS_i ] = linear_system_SIS( x_SIS, dx_SIS, delta_SIS( node ), node );
                
                rank_vs_N_SIS_tmp =  rank_vs_N_SIS_tmp +  rank( F_SIS_i );
            end
            rank_vs_N_SIS_this_parfor( N_all == N ) =  rank_vs_N_SIS_tmp/N;
            
            %Kuramoto
            rank_vs_N_kuramoto_tmp = 0;
            for node = 1:N
                [ ~, F_kuramoto_i ] = linear_system_kuramoto( x_kuramoto, dx_kuramoto, omega_kuramoto( node ), node );
                
                rank_vs_N_kuramoto_tmp =  rank_vs_N_kuramoto_tmp +  rank( F_kuramoto_i );
            end
            rank_vs_N_kuramoto_this_parfor( N_all == N ) =  rank_vs_N_kuramoto_tmp/N;
            
            %Cowan-Wilson
            rank_vs_N_cw_tmp = 0;
            for node = 1:N
                
                [ ~, F_cw_i ] = linear_system_cowan_wilson(  x_cw, dx_cw, tau_cw, mu_cw, node );
                
                rank_vs_N_cw_tmp =  rank_vs_N_cw_tmp +  rank( F_cw_i );
            end
            rank_vs_N_cw_this_parfor( N_all == N ) =  rank_vs_N_cw_tmp/N;
        end
        
        %save the results
        rank_vs_N_LV( curr_network, : ) = rank_vs_N_LV_this_parfor;
        rank_vs_N_MP( curr_network, : ) = rank_vs_N_MP_this_parfor;
        rank_vs_N_MM( curr_network, : ) = rank_vs_N_MM_this_parfor;
        rank_vs_N_SIS( curr_network, : ) = rank_vs_N_SIS_this_parfor;
        rank_vs_N_kuramoto( curr_network, : ) = rank_vs_N_kuramoto_this_parfor;
        rank_vs_N_cw( curr_network, : ) = rank_vs_N_cw_this_parfor;
    end
    
    %obtain the computation time
    computation_time = toc;
    
    %display the computation time
    fprintf( 'time: less than %i Minutes \n', ceil( computation_time/60 ) )
    
    %save the results
    save( filename )
end

end

