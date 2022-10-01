function Fig_S10_S14( seed, N, num_observations, random_network_model )

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
filename = strcat( './results/Fig_S10_S14_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', random_network_model );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    %set the simulation parameters
    parameters = set_parameters;
    
    %set the maximum prediction time
    parameters.T_max = [ 8  ...% LV
        2  ...% MP
        4  ...% MM
        2 ...% SIS
        2 ...%KUR
        4 ]; %CW
    
    %set the maximum prediction time T_max by a factor T_max_mult higher
    %(Figs. S10 and S14 explore various large ranges of the prediction time horizon,
    %which requires T_max to be large.)
    T_max_mult = 10;
    
    %set the maximum initial nodal state
    parameters.LV.x_init_max = 1;
    parameters.MP.x_init_max = 1;
    parameters.MM.x_init_max = 1;
    parameters.SIS.x_init_max = 0.1;
    parameters.kuramoto.x_init_max = pi/4;
    parameters.cw.x_init_max = 1;
    
    %the number of networks
    num_networks = 1e2;
    
    %the number of different observation times t_obs
    num_T_obs = 10;
    
    %generate the set of different observation times t_obs
    T_obs_fraction_all = linspace( 0.05, 3,  num_T_obs )/T_max_mult;
    
    %set the number of observations n
    num_observations = ceil( num_observations*T_max_mult );
    
    %pre-allocate the network reconstruction accuracy AUC
    AUC = nan( num_T_obs, num_networks, 6  );
    
    %pre-allocate the number of connected components
    connected_comp_all = nan( num_T_obs, num_networks, 6  );
    
    %pre-allocate the run time
    time_all = nan( num_T_obs, num_networks, 6 );
    
    %pre-allocate the prediction errors
    error_pred = nan( num_T_obs, num_observations, num_networks, 6  );
        
    %loop over all networks
    for network_i = 1:num_networks
        
        %the zero-one (unweighted) network
        A = create_network( N, random_network_model, [], [], [],  parameters.m0_BA,  parameters.m_BA, true );
        
        %the link weights
        link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
        
        %the weighted network
        B = zeros( N );
        B( A > 0 ) = link_weights;
        
        %loop over all observation time fractions t_obs/T_max
        for T_obs_i = 1:num_T_obs
            
            T_obs_fraction = T_obs_fraction_all( T_obs_i );
            
            %loop over all models
            for model_count = 1:6
                
                %all time samples
                t_all = linspace( 0, parameters.T_max( model_count ), num_observations + 1 );
                
                %the sampling time
                delta_T = t_all( 2 ) - t_all( 1 );
                
                %the observation time samples
                t_obs = t_all( 1:ceil( T_obs_fraction*( num_observations + 1) ));
                
                %the prediction time samples
                t_prediction = t_all( ceil( T_obs_fraction*( num_observations + 1) )+1:end );
                
                
                switch model_count
                    case 1
                        model_name = 'LV';
                        
                        %the initial nodal state
                        x_init = parameters.LV.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        alpha = 1 + parameters.LV.sigma_alpha*( 2*rand( N, 1 ) - 1);
                        theta = 1 + parameters.LV.sigma_theta*( 2*rand( N, 1 ) - 1);
                        
                        %generate the past nodal state sequence
                        [ x, ~ ] = compute_nodal_states_lotka_volterra( x_init, alpha, theta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] = compute_nodal_states_lotka_volterra( x( :, end ), alpha, theta, B, t_prediction);
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                        input_param.alpha = alpha;
                        input_param.theta = theta;
                    case 2
                        model_name = 'MP';
                        
                        %the initial nodal state
                        x_init = parameters.MP.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        alpha = 1 + parameters.MP.sigma_alpha*( 2*rand( N, 1 ) - 1 );
                        theta = 1 + parameters.MP.sigma_theta*( 2*rand( N, 1 ) - 1 );
                        
                        %generate the past nodal state sequence
                        [ x, ~ ] = compute_nodal_states_mutualistic_pop( x_init, alpha, theta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] = compute_nodal_states_mutualistic_pop( x( :, end ), alpha, theta, B, t_prediction );
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                        input_param.alpha = alpha;
                        input_param.theta = theta;
                    case 3
                        model_name = 'MM';
                        
                        %the initial nodal state
                        x_init = parameters.MM.x_init_max*rand( N, 1 );
                        
                        %generate the past nodal state sequence
                        [ x, ~ ] = compute_nodal_states_michaelis_menten( x_init, parameters.MM.hill_coeff, B, t_obs  );
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] =compute_nodal_states_michaelis_menten( x( :, end ), parameters.MM.hill_coeff, B, t_prediction  );
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                    case 4
                        model_name = 'SIS';
                        
                        %the initial nodal state
                        x_init = parameters.SIS.x_init_max*rand( N, 1 );
                        
                        %the parameters of the model
                        delta_init = 1 + parameters.SIS.sigma_delta*( 2*rand( N, 1 ) - 1);
                        R_0_init = eigs( diag( 1./sqrt( delta_init ))*B*diag( 1./sqrt( delta_init )), 1 );
                        delta = R_0_init./parameters.SIS.R_0_SIS*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B
                        
                        %generate the past nodal state sequence
                        [ x, ~] = compute_nodal_states_SIS( x_init, delta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] = compute_nodal_states_SIS( x( :, end ), delta, B, t_prediction );
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                        input_param.delta = delta;
                    case 5
                        model_name = 'kuramoto';
                        
                        %the initial nodal state
                        x_init = ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )*rand( N, 1 ) - ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )/2;
                        
                        %the parameters of the model
                        omega = parameters.kuramoto.sigma_omega*randn( N, 1 );
                        
                        %generate the past nodal state sequence
                        [ x, ~] = compute_nodal_states_kuramoto( x_init, omega, B, t_obs);
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] = compute_nodal_states_kuramoto( x( :, end ), omega, B, t_prediction);
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                        input_param.omega = omega;
                    case 6
                        model_name = 'cw';
                        
                        %the initial nodal state
                        x_init = parameters.cw.x_init_max*rand( N, 1 );
                        
                        %generate the past nodal state sequence
                        [ x, ~] = compute_nodal_states_cowan_wilson( x_init, parameters.cw.tau, parameters.cw.mu, B, t_obs  );
                        
                        %generate the future nodal state sequence
                        [ x_future, ~ ] = compute_nodal_states_cowan_wilson( x( :, end ), parameters.cw.tau, parameters.cw.mu, B, t_prediction  );
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                    otherwise
                        error('unknown model')
                end
                
                %numerical differentiation
                dx = transpose( diff( transpose( x ) ) )/delta_T;
                
                %delete last nodal state such that x and dx have the same size
                x = x( :, 1:end-1 );
                
                %start the network reconstruction
                tic
                if model_count <6
                    [ B_hat, ~ ] = network_reconstruction( x, dx, model_name, input_param, parameters );
                else
                    [ B_hat, ~ ] = network_reconstruction( x, dx, 'cowan_wilson', input_param, parameters );
                end
                
                %save the runtime
                time_all(  T_obs_i, network_i, model_count ) = toc;
                
                %compute and save the AUC
                [~,~,~,AUC_i]=perfcurve( B( : )>0, B_hat( : ),1);
                AUC(  T_obs_i, network_i, model_count ) = AUC_i;
                
                %compute and save the number of connected components
                [connectedComponents, ~] = graphconncomp(sparse(B_hat), 'Directed', false);
                connected_comp_all(  T_obs_i, network_i, model_count ) = connectedComponents;
                
                %predict the nodal state with the surrogate network B_hat
                switch model_count
                    case 1
                        x_future_hat = compute_nodal_states_lotka_volterra( x_future( :, 1 ), alpha, theta, B_hat, t_prediction  );
                    case 2
                        x_future_hat = compute_nodal_states_mutualistic_pop( x_future( :, 1 ), alpha, theta, B_hat, t_prediction  );
                    case 3
                        x_future_hat = compute_nodal_states_michaelis_menten( x_future( :, 1 ), parameters.MM.hill_coeff, B_hat, t_prediction  );
                    case 4
                        x_future_hat = compute_nodal_states_SIS( x_future( :, 1 ), delta, B_hat, t_prediction);
                    case 5
                        x_future_hat = compute_nodal_states_kuramoto( x_future( :, 1 ), omega, B_hat, t_prediction);
                    case 6
                        x_future_hat = compute_nodal_states_cowan_wilson( x_future( :, 1 ), parameters.cw.tau, parameters.cw.mu, B_hat, t_prediction);
                end
                
                %save the prediction error
                error_pred(  T_obs_i, 1:length( t_prediction ), network_i, model_count ) = mean( abs( x_future_hat - x_future ) );
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% save file  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save( strcat( filename, '.mat' ) )
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% plots  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fig S10: Prediction error versus n_future
figure
for model_count = 1:6
    switch model_count
        case 1
            model_name = 'LV';
        case 2
            model_name = 'MP';
        case 3
            model_name = 'MM';
        case 4
            model_name = 'SIS';
        case 5
            model_name = 'kuramoto';
        case 6
            model_name = 'cw';
    end
    
    subplot( 3, 2, model_count )
    
    T_obs_ind = 1:length( T_obs_fraction_all );
    
    error_pred_mean = mean( error_pred( T_obs_ind, :, :, model_count ), 3 );
    loglog( error_pred_mean' )
    xlabel('Prediction time horizon $$\tilde{t}/\Delta t $$', 'Interpreter', 'LaTeX')
    ylabel( 'Prediction error $$\epsilon(\tilde{t})$$' , 'Interpreter', 'LaTeX')
    
    if model_count==1
        legend("t_{obs}/T_{max} = " + string(round( T_obs_fraction_all( T_obs_ind ), 2 ) ), 'Location', 'best', 'NumColumns',ceil( length(T_obs_ind)/2))
    end
    
    title( model_name )
end

%Fig S14: Histogramm of number of connected components
figure
for model_count = 1:6
    switch model_count
        case 1
            model_name = 'LV';
        case 2
            model_name = 'MP';
        case 3
            model_name = 'MM';
        case 4
            model_name = 'SIS';
        case 5
            model_name = 'kuramoto';
        case 6
            model_name = 'cw';
    end
    subplot( 3, 2, model_count )
    tmp_data = connected_comp_all( :, :, model_count );
    
    frac_connected = sum( tmp_data( : ) == 1 )/length( tmp_data(:) );
    
    histogram( tmp_data( : ),'Normalization','pdf')
    title(  model_name )
    title( strcat( model_name, '; connected: ', num2str( frac_connected )) )
    ylabel( 'Frequency' )
    xlabel( 'Connected components in $$\hat{A}$$' , 'Interpreter', 'LaTeX')
    
end


end

