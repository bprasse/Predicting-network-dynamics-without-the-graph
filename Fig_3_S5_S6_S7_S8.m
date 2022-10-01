function Fig_3_S5_S6_S7_S8( seed, N, num_observations, network_type )

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

%only for plotting: num_N_plot is the number of nodal traces that are plotted
num_N_plot = 6;

%the filename to save/read the results
filename = strcat( './results/Fig_3_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', network_type );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% set the simulation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    parameters = set_parameters;
    
    %set the ratio t_obs/T_max:
    T_obs_fraction = 1/5;
    
    %the link probability for ER networks
    p_ER = 0.05;
    
    %the parameters for the BA networks
    m0_BA = 3;
    m_BA = 3;
    
    %save the results of the different models here:
    results_LV = [];
    results_SIS = [];
    results_MP = [];
    results_MM = [];
    results_kuramoto = [];
    results_cw = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 1: Lotka-Volterra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_LV = create_network( [], 'little_rock_lake',[], [],[], [], [], [] );
        case 'ER'
            A_LV = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
            parameters.LV.x_init_max = 0.1*parameters.LV.x_init_max;
        case 'BA'
            A_LV = create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
    end
    
    %the network size/number of nodes
    N_LV = size( A_LV, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_LV ), 1);
    
    %the weighted network
    results_LV.B = zeros( N_LV );
    results_LV.B( A_LV > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_LV.spectrum = zeros( 1, N_LV );
    results_LV.spectrum_hat = zeros( 1, N_LV );
    
    %the parameters of the model
    results_LV.alpha = 1 + parameters.LV.sigma_alpha*( 2*rand( N_LV, 1 ) - 1);
    results_LV.theta = 1 + parameters.LV.sigma_theta*( 2*rand( N_LV, 1 ) - 1);
    
    %the initial nodal state
    results_LV.x_init = parameters.LV.x_init_max*rand( N_LV, 1 );
    
    %the observation time samples
    results_LV.t_obs = linspace( 0, parameters.T_max( 1 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_LV = results_LV.t_obs( 2 ) - results_LV.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_LV = results_LV.t_obs( end ):delta_T_LV:parameters.T_max( 1 );
    
    %generate the past nodal state sequence
    [ x_LV, ~ ] = compute_nodal_states_lotka_volterra( results_LV.x_init, results_LV.alpha, results_LV.theta, results_LV.B, results_LV.t_obs );
    
    %generate the future nodal state sequence
    [ x_LV_future, ~ ] = compute_nodal_states_lotka_volterra( x_LV( :, end ), results_LV.alpha, results_LV.theta, results_LV.B, t_prediction_LV );
    
    %numerical differentiation
    dx_LV = transpose( diff( transpose( x_LV ) ) )/delta_T_LV;
    
    %delete last nodal state such that x and dx have the same size
    x_LV = x_LV( :, 1:end-1 );
    results_LV.t_obs( end ) = [];
    
    %the input parameters to the network reconstruction
    input_param = [];
    input_param.alpha = results_LV.alpha;
    input_param.theta = results_LV.theta;
    
    %start the network reconstruction
    fprintf('Running Lotka-Volterra network reconstruction (N = %i)...\n', N_LV)
    [ B_hat_LV, rho_LV ] = network_reconstruction( x_LV, dx_LV, 'LV', input_param, parameters );
    fprintf( 'Finished Lotka-Volterra network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_LV_hat = compute_nodal_states_lotka_volterra( results_LV.x_init, results_LV.alpha, results_LV.theta, B_hat_LV, results_LV.t_obs  );
    x_LV_future_hat = compute_nodal_states_lotka_volterra( x_LV( :, end ), results_LV.alpha, results_LV.theta, B_hat_LV, t_prediction_LV  );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_LV.B( : )>0, B_hat_LV( : ),1);
    results_LV.AUC = AUC;
    
    %save the surrogate matrix
    results_LV.B_hat = B_hat_LV;
    
    %save the regularisation parameter
    results_LV.rho = rho_LV;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_LV.spectrum = sort( abs( eig( results_LV.B ) ));
    results_LV.spectrum_hat = sort( abs( eig( B_hat_LV ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_LV.e_fit =  mean(mean( abs( x_LV - x_LV_hat )));
    results_LV.e_pred = mean(mean( abs( x_LV_future - x_LV_future_hat )));
    
    %display results
    fprintf( '-> AUC = %.2g \n', results_LV.AUC );
    fprintf( '-> e_fit = %.2d \n', results_LV.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_LV.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 2: Mutualistic population
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_MP = create_network( [], 'mutualistic',[], [],[], [], [], [] );
        case 'ER'
            A_MP = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
            parameters.T_max( 2 ) = 10*parameters.T_max( 2 );
            parameters.MP.x_init_max = 0.2*parameters.MP.x_init_max;
        case 'BA'
            A_MP =  create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
            parameters.T_max( 2 ) = 10*parameters.T_max( 2 );
            parameters.MP.x_init_max = 0.1*parameters.MP.x_init_max;
    end
    
    %the network size/number of nodes
    N_MP = size( A_MP, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_MP ), 1);
    
    %the weighted network
    results_MP.B = zeros( N_MP );
    results_MP.B ( A_MP > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_MP.spectrum = zeros( 1, N_MP );
    results_MP.spectrum_hat = zeros( 1, N_MP );
    
    %the parameters of the model
    results_MP.alpha = 1 + parameters.MP.sigma_alpha*( 2*rand( N_MP, 1 ) - 1 );
    results_MP.theta = 1 + parameters.MP.sigma_theta*( 2*rand( N_MP, 1 ) - 1 );
    
    %the initial nodal state
    results_MP.x_init = parameters.MP.x_init_max*rand( N_MP, 1 );
    
    %the observation time samples
    results_MP.t_obs = linspace( 0, parameters.T_max( 2 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_MP = results_MP.t_obs( 2 ) - results_MP.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_MP = results_MP.t_obs( end ):delta_T_MP:parameters.T_max( 2 );
    
    %generate the past nodal state sequence
    [ x_MP, ~ ] = compute_nodal_states_mutualistic_pop( results_MP.x_init, results_MP.alpha, results_MP.theta, results_MP.B, results_MP.t_obs );
    
    %generate the future nodal state sequence
    [ x_MP_future, ~ ] = compute_nodal_states_mutualistic_pop( x_MP( :, end ), results_MP.alpha, results_MP.theta, results_MP.B, t_prediction_MP  );
    
    %numerical differentiation
    dx_MP = transpose( diff( transpose( x_MP ) ) )/delta_T_MP;
    
    %delete last nodal state such that x and dx have the same size
    x_MP = x_MP( :, 1:end-1 );
    results_MP.t_obs( end ) = [];
    
    %the input parameters to the network reconstruction
    input_param = [];
    input_param.alpha = results_MP.alpha;
    input_param.theta = results_MP.theta;
    
    %start the network reconstruction
    fprintf('Running mutualistic population dynamics network reconstruction (N = %i)...\n', N_MP)
    [ B_hat_MP, rho_MP ] = network_reconstruction( x_MP, dx_MP, 'MP', input_param, parameters );
    fprintf( 'Finished Lotka-Volterra network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_MP_hat = compute_nodal_states_mutualistic_pop( results_MP.x_init, results_MP.alpha, results_MP.theta, B_hat_MP, results_MP.t_obs );
    x_MP_future_hat = compute_nodal_states_mutualistic_pop( x_MP( :, end ), results_MP.alpha, results_MP.theta, B_hat_MP, t_prediction_MP  );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_MP.B( : )>0, B_hat_MP( : ),1);
    results_MP.AUC = AUC;
    
    %save the surrogate matrix
    results_MP.B_hat = B_hat_MP;
    
    %save the regularisation parameter
    results_MP.rho = rho_MP;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_MP.spectrum = sort( abs( eig( results_MP.B ) ));
    results_MP.spectrum_hat = sort( abs( eig( B_hat_MP ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_MP.e_fit =  mean(mean( abs( x_MP - x_MP_hat )));
    results_MP.e_pred =  mean(mean( abs( x_MP_future - x_MP_future_hat )));
    
    %display results
    fprintf( 'Finished mutualistic population dynamics network reconstruction: \n' );
    fprintf( '-> AUC = %.2g \n', results_MP.AUC );
    fprintf( '-> e_fit = %.2d \n', results_MP.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_MP.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 3: Michaelis-Menten
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_MM = create_network( [], 'yeast_regulatory',[], [],[], [], [], [] );
        case 'ER'
            A_MM = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
        case 'BA'
            A_MM = create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
    end
    
    %the network size/number of nodes
    N_MM = size( A_MM, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_MM ), 1);
    
    %the weighted network
    results_MM.B = zeros( N_MM );
    results_MM.B( A_MM > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_MM.spectrum = zeros( 1, N_MM );
    results_MM.spectrum_hat = zeros( 1, N_MM );
    
    %the initial nodal state
    results_MM.x_init = parameters.MM.x_init_max*rand( N_MM, 1 );
    
    %the observation time samples
    results_MM.t_obs = linspace( 0, parameters.T_max( 3 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_MM = results_MM.t_obs( 2 ) - results_MM.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_MM = results_MM.t_obs( end ):delta_T_MM:parameters.T_max( 3 );
    
    %generate the past nodal state sequence
    [ x_MM, ~ ] = compute_nodal_states_michaelis_menten( results_MM.x_init, parameters.MM.hill_coeff, results_MM.B, results_MM.t_obs  );
    
    %generate the future nodal state sequence
    [ x_MM_future, ~ ] = compute_nodal_states_michaelis_menten( x_MM( :, end ), parameters.MM.hill_coeff, results_MM.B, t_prediction_MM );
    
    %numerical differentiation
    dx_MM = transpose( diff( transpose( x_MM ) ) )/delta_T_MM;
    
    %delete last nodal state such that x and dx have the same size
    x_MM = x_MM( :, 1:end-1 );
    results_MM.t_obs( end ) = [];
    
    %start the network reconstruction
    fprintf('Running Michaelis-Menten network reconstruction (N = %i)...\n', N_MM)
    [ B_hat_MM, rho_MM ] = network_reconstruction( x_MM, dx_MM, 'MM', [], parameters );
    fprintf( 'Finished Michaelis-Menten network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_MM_hat = compute_nodal_states_michaelis_menten( results_MM.x_init, parameters.MM.hill_coeff, B_hat_MM, results_MM.t_obs  );
    x_MM_future_hat = compute_nodal_states_michaelis_menten( x_MM( :, end ), parameters.MM.hill_coeff, B_hat_MM, t_prediction_MM  );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_MM.B( : )>0, B_hat_MM( : ),1);
    results_MM.AUC = AUC;
    
    %save the surrogate matrix
    results_MM.B_hat = B_hat_MM;
    
    %save the regularisation parameter
    results_MM.rho = rho_MM;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_MM.spectrum = sort( abs( eig( results_MM.B ) ));
    results_MM.spectrum_hat = sort( abs( eig( B_hat_MM ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_MM.e_fit =  mean(mean( abs( x_MM - x_MM_hat )));
    results_MM.e_pred = mean(mean( abs( x_MM_future - x_MM_future_hat )));
    
    %display results
    fprintf( '-> AUC = %.2g \n', results_MM.AUC );
    fprintf( '-> e_fit = %.2d \n', results_MM.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_MM.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 4: SIS epidemics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_SIS = create_network( [], 'infectious',[], [],[], [], [], [] );
        case 'ER'
            A_SIS = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
        case 'BA'
            A_SIS = create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
    end
    
    %the network size/number of nodes
    N_SIS = size( A_SIS, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_SIS ), 1);
    
    %the weighted network
    results_SIS.B = zeros( N_SIS );
    results_SIS.B( A_SIS > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_SIS.spectrum = zeros( 1, N_SIS );
    results_SIS.spectrum_hat = zeros( 1, N_SIS );
    
    %the parameters of the model
    delta_init = 1 + parameters.SIS.sigma_delta*( 2*rand( N_SIS, 1 ) - 1);
    R_0_init = eigs( diag( 1./sqrt( delta_init ))*results_SIS.B*diag( 1./sqrt( delta_init )), 1 );
    multiplicity_delta = R_0_init./parameters.SIS.R_0_SIS;
    results_SIS.delta = multiplicity_delta*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B
    
    %the initial nodal state
    results_SIS.x_init = parameters.SIS.x_init_max*rand( N_SIS, 1 );
    
    %the observation time samples
    results_SIS.t_obs = linspace( 0, parameters.T_max( 4 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_SIS = results_SIS.t_obs( 2 ) - results_SIS.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_SIS = results_SIS.t_obs( end ):delta_T_SIS:parameters.T_max( 4 );
    
    %generate the past nodal state sequence
    [ x_SIS, ~ ] = compute_nodal_states_SIS( results_SIS.x_init, results_SIS.delta, results_SIS.B, results_SIS.t_obs );
    
    %generate the future nodal state sequence
    [ x_SIS_future, ~ ] = compute_nodal_states_SIS( x_SIS( :, end ), results_SIS.delta, results_SIS.B, t_prediction_SIS );
    
    %numerical differentiation
    dx_SIS = transpose( diff( transpose( x_SIS ) ) )/delta_T_SIS;
    
    %delete last nodal state such that x and dx have the same size
    x_SIS = x_SIS( :, 1:end-1 );
    results_SIS.t_obs( end ) = [];
    
    %the input parameters to the network reconstruction
    input_param = [];
    input_param.delta = results_SIS.delta;
    
    %start the network reconstruction
    fprintf('Running SIS epidemics network reconstruction (N = %i)...\n', N_SIS)
    [ B_hat_SIS, rho_SIS ] = network_reconstruction( x_SIS, dx_SIS, 'SIS', input_param, parameters );
    fprintf( 'Finished SIS epidemics network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_SIS_hat = compute_nodal_states_SIS( results_SIS.x_init, results_SIS.delta, B_hat_SIS, results_SIS.t_obs );
    x_SIS_future_hat = compute_nodal_states_SIS( x_SIS( :, end ), results_SIS.delta, B_hat_SIS, t_prediction_SIS );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_SIS.B( : )>0, B_hat_SIS( : ),1);
    results_SIS.AUC = AUC;
    
    %save the surrogate matrix
    results_SIS.B_hat = B_hat_SIS;
    
    %save the regularisation parameter
    results_SIS.rho = rho_SIS;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_SIS.spectrum = sort( abs( eig( results_SIS.B ) ));
    results_SIS.spectrum_hat = sort( abs( eig( B_hat_SIS ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_SIS.e_fit =  mean(mean( abs( x_SIS - x_SIS_hat )));
    results_SIS.e_pred =  mean(mean( abs( x_SIS_future - x_SIS_future_hat )));
    
    %display results
    fprintf( '-> AUC = %.2g \n', results_SIS.AUC );
    fprintf( '-> e_fit = %.2d \n', results_SIS.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_SIS.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 5: Kuramoto oscillators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_kuramoto = create_network( [], 'human_connectome',[], [],[], [], [], [] );
        case 'ER'
            A_kuramoto = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
        case 'BA'
            A_kuramoto = create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
    end
    
    %the network size/number of nodes
    N_kuramoto = size( A_kuramoto, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_kuramoto ), 1);
    
    %the weighted network
    results_kuramoto.B = zeros( N_kuramoto );
    results_kuramoto.B( A_kuramoto > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_kuramoto.spectrum = zeros( 1, N_kuramoto );
    results_kuramoto.spectrum_hat = zeros( 1, N_kuramoto );
    
    %the parameters of the model
    results_kuramoto.omega = parameters.kuramoto.sigma_omega*randn( N_kuramoto, 1 );
    
    %the initial nodal state
    results_kuramoto.x_init = ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )*rand( N_kuramoto, 1 ) - ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )/2;
    
    %the observation time samples
    results_kuramoto.t_obs = linspace( 0, parameters.T_max( 5 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_kuramoto = results_kuramoto.t_obs( 2 ) - results_kuramoto.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_kuramoto = results_kuramoto.t_obs( end ):delta_T_kuramoto:parameters.T_max( 5 );
    
    %generate the past nodal state sequence
    [ x_kuramoto, ~ ] = compute_nodal_states_kuramoto( results_kuramoto.x_init, results_kuramoto.omega, results_kuramoto.B, results_kuramoto.t_obs  );
    
    %generate the future nodal state sequence
    [ x_kuramoto_future, ~ ] = compute_nodal_states_kuramoto( x_kuramoto( :, end ), results_kuramoto.omega, results_kuramoto.B, t_prediction_kuramoto );
    
    %numerical differentiation
    dx_kuramoto = transpose( diff( transpose( x_kuramoto ) ) )/delta_T_kuramoto;
    
    %delete last nodal state such that x and dx have the same size
    x_kuramoto = x_kuramoto( :, 1:end-1 );
    results_kuramoto.t_obs( end ) = [];
    
    %the input parameters to the network reconstruction
    input_param = [];
    input_param.omega = results_kuramoto.omega;
    
    %start the network reconstruction
    fprintf('Running Kuramoto network reconstruction (N = %i)...\n', N_kuramoto)
    [ B_hat_kuramoto, rho_kuramoto ] = network_reconstruction( x_kuramoto, dx_kuramoto, 'kuramoto', input_param, parameters );
    fprintf( 'Finished Kuramoto network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_kuramoto_hat = compute_nodal_states_kuramoto( results_kuramoto.x_init, results_kuramoto.omega, B_hat_kuramoto, results_kuramoto.t_obs  );
    x_kuramoto_future_hat = compute_nodal_states_kuramoto( x_kuramoto( :, end ), results_kuramoto.omega, B_hat_kuramoto, t_prediction_kuramoto  );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_kuramoto.B( : )>0, B_hat_kuramoto( : ),1);
    results_kuramoto.AUC = AUC;
    
    %save the surrogate matrix
    results_kuramoto.B_hat = B_hat_kuramoto;
    
    %save the regularisation parameter
    results_kuramoto.rho = rho_kuramoto;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_kuramoto.spectrum = sort( abs( eig( results_kuramoto.B ) ));
    results_kuramoto.spectrum_hat = sort( abs( eig( B_hat_kuramoto ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_kuramoto.e_fit =  mean(mean( abs( x_kuramoto - x_kuramoto_hat )));
    results_kuramoto.e_pred = mean(mean( abs( x_kuramoto_future - x_kuramoto_future_hat )));
    
    %display results
    fprintf( '-> AUC = %.2g \n', results_kuramoto.AUC );
    fprintf( '-> e_fit = %.2d \n', results_kuramoto.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_kuramoto.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Model 6: Cowan-Wilson
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the zero-one (unweighted) network
    switch network_type
        case 'real_world'
            A_cw = create_network( [], 'c_elegans_neuronal',[], [],[], [], [], [] );
        case 'ER'
            A_cw = create_network( N, 'ER_undirected', p_ER, [],[], [], [], [] );
        case 'BA'
            A_cw = create_network( N, 'BA',[], [],[], m0_BA, mm_BA, [] );
    end
        
    %the network size/number of nodes
    N_cw = size( A_cw, 1 );
    
    %the link weights
    link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_cw ), 1);
    
    %the weighted network
    results_cw.B = zeros( N_cw );
    results_cw.B( A_cw > 0 ) = link_weights;
    
    %pre-allocate the spectrum of the true and the surrogate matrix
    results_cw.spectrum = zeros( 1, N_cw );
    results_cw.spectrum_hat = zeros( 1, N_cw );
    
    %the initial nodal state
    results_cw.x_init = parameters.cw.x_init_max*rand( N_cw, 1 );
    
    %the observation time samples
    results_cw.t_obs = linspace( 0, parameters.T_max( 6 )*T_obs_fraction, num_observations + 1 );
    
    %the sampling time
    delta_T_cw = results_cw.t_obs( 2 ) - results_cw.t_obs( 1 );
    
    %the prediction time samples
    t_prediction_cw = results_cw.t_obs( end ):delta_T_cw:parameters.T_max( 6 );
    
    %generate the past nodal state sequence
    [ x_cw, ~  ] = compute_nodal_states_cowan_wilson( results_cw.x_init, parameters.cw.tau, parameters.cw.mu, results_cw.B, results_cw.t_obs  );
    
    %generate the future nodal state sequence
    [ x_cw_future, ~ ] = compute_nodal_states_cowan_wilson( x_cw( :, end ), parameters.cw.tau, parameters.cw.mu, results_cw.B, t_prediction_cw  );
    
    %numerical differentiation
    dx_cw = transpose( diff( transpose( x_cw ) ) )/delta_T_cw;
    
    %delete last nodal state such that x and dx have the same size
    x_cw = x_cw( :, 1:end-1 );
    results_cw.t_obs( end ) = [];
    
    %start the network reconstruction
    fprintf('Running Cowan-Wilson network reconstruction (N = %i)...\n', N_cw)
    [ B_hat_cw, rho_cw ] = network_reconstruction( x_cw, dx_cw, 'cowan_wilson', [], parameters );
    fprintf( 'Finished Cowan-Wilson network reconstruction: \n' );
    
    %predict the nodal state with the surrogate network B_hat
    x_cw_hat = compute_nodal_states_cowan_wilson( results_cw.x_init, parameters.cw.tau, parameters.cw.mu,  B_hat_cw, results_cw.t_obs  );
    x_cw_future_hat = compute_nodal_states_cowan_wilson( x_cw( :, end ), parameters.cw.tau, parameters.cw.mu, B_hat_cw, t_prediction_cw  );
    
    %obtain and save the AUC of the surrogate network and the true network
    [~,~,~,AUC]=perfcurve( results_cw.B( : )>0, B_hat_cw( : ),1);
    results_cw.AUC = AUC;
    
    %save the surrogate matrix
    results_cw.B_hat_cw_hat = B_hat_cw;
    
    %save the regularisation parameter
    results_cw.rho = rho_cw;
    
    %save the eigenvalue spectra of the true and surrogate network
    results_cw.spectrum = sort( abs( eig( results_cw.B ) ));
    results_cw.spectrum_hat = sort( abs( eig( B_hat_cw ) ));
    
    %the error of the fit and the prediction based on the surrogate network
    results_cw.e_fit =  mean(mean( abs( x_cw - x_cw_hat )));
    results_cw.e_pred = mean(mean( abs( x_cw_future - x_cw_future_hat )));
    
    %display results
    fprintf( '-> AUC = %.2g \n', results_cw.AUC );
    fprintf( '-> e_fit = %.2d \n', results_cw.e_fit  );
    fprintf( '-> e_pred = %.2d \n', results_cw.e_pred  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% save results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save( strcat( filename, '.mat' ) )
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot results, corresponding to:
%Fig 3 (network_type == 'real_world' ); or
%Fig S5 (network_type == 'BA' ); or
%Fig S6 (network_type == 'ER' )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a new large figure 
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over all models 
for model_count = 1:6
    
    switch model_count
        case 1
            model_name = 'LV';
            x = x_LV;
            x_future = x_LV_future;
            x_hat = x_LV_hat;
            x_future_hat = x_LV_future_hat;
            t_plot = [results_LV.t_obs t_prediction_LV];
            AUC = results_LV.AUC;
            e_pred = results_LV.e_pred;
            t_obs_end = results_LV.t_obs( end );
        case 2
            model_name = 'MP';
            x = x_MP;
            x_future = x_MP_future;
            x_hat = x_MP_hat;
            x_future_hat = x_MP_future_hat;
            t_plot = [results_MP.t_obs t_prediction_MP];
            AUC = results_MP.AUC;
            e_pred = results_MP.e_pred;
            t_obs_end = results_MP.t_obs( end );
        case 3
            model_name = 'MM';
            x = x_MM;
            x_future = x_MM_future;
            x_hat = x_MM_hat;
            x_future_hat = x_MM_future_hat;
            t_plot = [results_MM.t_obs t_prediction_MM];
            AUC = results_MM.AUC;
            e_pred = results_MM.e_pred;
            t_obs_end = results_MM.t_obs( end );
        case 4
            model_name = 'SIS';
            x = x_SIS;
            x_future = x_SIS_future;
            x_hat = x_SIS_hat;
            x_future_hat = x_SIS_future_hat;
            t_plot = [results_SIS.t_obs t_prediction_SIS];
            AUC = results_SIS.AUC;
            e_pred = results_SIS.e_pred;
            t_obs_end = results_SIS.t_obs( end );
        case 5
            model_name = 'KUR';
            x = x_kuramoto;
            x_future = x_kuramoto_future;
            x_hat = x_kuramoto_hat;
            x_future_hat = x_kuramoto_future_hat;
            t_plot = [results_kuramoto.t_obs t_prediction_kuramoto];
            AUC = results_kuramoto.AUC;
            e_pred = results_kuramoto.e_pred;
            t_obs_end = results_kuramoto.t_obs( end );
        case 6
            model_name = 'WC';
            x = x_cw;
            x_future = x_cw_future;
            x_hat = x_cw_hat;
            x_future_hat = x_cw_future_hat;
            t_plot = [results_cw.t_obs t_prediction_cw];
            AUC = results_cw.AUC;
            e_pred = results_cw.e_pred;
            t_obs_end = results_cw.t_obs( end );
    end
    
    
    % obtain N_plot, the nodes whose state is plotted 
    N_plot = nan( 1, num_N_plot+1 );
    
    %x_end_quantised has N_plot+1 values, ranging from the minimum to the 
    %maximum nodal state at the last time point 
    x_end_quantised = linspace( min( x_future( :, end ) ), max( x_future( :, end ) ), num_N_plot+1 );
    
    %the nodes in N_plot are obtain by finding the nodes in x( :, end ) 
    %that are closes to the values in x_end_quantised
    for curr_N_plot = 1:(num_N_plot+1)
        [~, indd] = min( abs( x_future( :, end ) - x_end_quantised( curr_N_plot ) ));
        N_plot( curr_N_plot ) = indd;
    end
    
    
    subplot( 2, 3, model_count )
    plot( t_plot, [x( N_plot, : ) x_future( N_plot, : )]', 'b' )
    hold
    y1=get(gca,'ylim');
    plot( [ t_obs_end t_obs_end ], y1, '-om' );
    plot( t_plot, [x_hat( N_plot, : ) x_future_hat( N_plot, : )]', '--r' )
    title( [model_name, '; AUC = ', num2str( AUC,'%.2f' ), '; e_{pred}=', num2str( e_pred,'%.2f') ])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot results, corresponding to:
%Fig S7 (network_type == 'BA' ); or
%Fig S8 (network_type == 'ER' ):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a new large figure 
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over all models 
for model_count = 1:6
    
    switch model_count
        case 1
            model_name = 'LV';
            results = results_LV;
        case 2
            model_name = 'MP';
            results = results_MP;
        case 3
            model_name = 'MM';
            results = results_MM;
        case 4
            model_name = 'SIS';
            results = results_SIS;
        case 5
            model_name = 'KUR';
            results = results_kuramoto;
        case 6
            model_name = 'WC';
            results = results_cw;
    end
    
    subplot( 2, 3, model_count )
    stem( results.spectrum )
    hold
    stem( results.spectrum_hat, 'xr' )
    xlabel( 'Index i' )
    ylabel( 'Magnitude eigenvalue |\lambda_i|' )
    title( model_name )
end
end





