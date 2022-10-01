function Fig_S11_S18( seed, num_observations, random_network_model )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%set the number of observations n
if nargin < 2
    num_observations = 2e2;
end

%the network type:
if nargin < 3
    random_network_model = 'BA';
end

%the filename to save/read the results
filename = strcat( './results/Fig_S11_S18_seed_',  num2str( seed ), '_n_', num2str( num_observations ), '_', random_network_model );

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
    
    %set the maximum initial nodal state
    parameters.LV.x_init_max = 1;
    parameters.MP.x_init_max = 1;
    parameters.MM.x_init_max = 1;
    parameters.SIS.x_init_max = 0.1;
    parameters.kuramoto.x_init_max = pi/4;
    parameters.cw.x_init_max = 1;
    
    %the number of networks
    num_networks = 1e2;
    
    %the range of values for the network size N
    N_all = [500, 400, 300, 250, 200, 150, 100, 75, 50, 25];
    
    %the number of different values for N
    num_N = length( N_all );
    
    %the observation time fractions t_obs/T_max
    T_obs_fraction = 0.5;
    
    %pre-allocate the network reconstruction accuracy AUC
    AUC = nan( num_N, num_networks, 6 );
    
    %pre-allocate the run time
    time_all = nan( num_N, num_networks, 6 );
    
    %pre-allocate the prediction errors
    error_pred = nan( num_N, num_observations, num_networks, 6 );
    
    %obtain average dominant eigenvalue lambda_1 across the random networks
    %with N=100 nodes
    lambda_1_N_100 = 0;
    num_run_lambda_1 = 1e3;
    for random_network_i = 1:num_run_lambda_1
        
        %the zero-one (unweighted) network
        A = create_network( 100, random_network_model, [], [], [],  parameters.m0_BA,  parameters.m_BA, true );
        
        %the link weights
        link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
        
        %the weighted network
        B = zeros( 100 );
        B( A > 0 ) = link_weights;
        
        %add the current dominant eigenvalue of B to lambda_1_N_100
        lambda_1_N_100 = lambda_1_N_100 + eigs( B, 1 );
    end
    
    %obtain the mean of the dominant eigenvalue by dividing by num_run_lambda_1
    lambda_1_N_100 = lambda_1_N_100/num_run_lambda_1;
    
    %loop over num_networks realisations of the random network model
    for network_i = 1:num_networks
        
        %loop over all network sizes in N_all
        for N = N_all
            N_i = find( N == N_all );
            
            %the zero-one (unweighted) network
            A = create_network( N, random_network_model, [], [], [],  parameters.m0_BA,  parameters.m_BA, true );
            
            %the link weights
            link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
            
            %the weighted network
            B = zeros( N );
            B( A > 0 ) = link_weights;
            
            %scale the network B such that its dominant eigenvalue is in
            %[0.9, 1.1]*lambda_1_N_100
            B = lambda_1_N_100/eigs( B, 1 )*B*( 0.9 + 0.2*rand );
            
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
                        [ x, dx ] = compute_nodal_states_lotka_volterra( x_init, alpha, theta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] = compute_nodal_states_lotka_volterra( x( :, end ), alpha, theta, B, t_prediction);
                        
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
                        [ x, dx ] = compute_nodal_states_mutualistic_pop( x_init, alpha, theta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] = compute_nodal_states_mutualistic_pop( x( :, end ), alpha, theta, B, t_prediction );
                        
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
                        [ x, dx ] = compute_nodal_states_michaelis_menten( x_init, parameters.MM.hill_coeff, B, t_obs  );
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] =compute_nodal_states_michaelis_menten( x( :, end ), parameters.MM.hill_coeff, B, t_prediction  );
                        
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
                        [ x, dx] = compute_nodal_states_SIS( x_init, delta, B, t_obs );
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] = compute_nodal_states_SIS( x( :, end ), delta, B, t_prediction );
                        
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
                        [ x, dx] = compute_nodal_states_kuramoto( x_init, omega, B, t_obs);
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] = compute_nodal_states_kuramoto( x( :, end ), omega, B, t_prediction);
                        
                        %the input parameters for the network
                        %reconstruction algorithm
                        input_param = [];
                        input_param.omega = omega;
                    case 6
                        model_name = 'cw';
                        
                        %the initial nodal state
                        x_init = parameters.cw.x_init_max*rand( N, 1 );
                        
                        %generate the past nodal state sequence
                        [ x, dx] = compute_nodal_states_cowan_wilson( x_init, parameters.cw.tau, parameters.cw.mu, B, t_obs  );
                        
                        %generate the future nodal state sequence
                        [ x_future, dx_future ] = compute_nodal_states_cowan_wilson( x( :, end ), parameters.cw.tau, parameters.cw.mu, B, t_prediction  );
                        
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
                time_all( N_i, network_i, model_count ) = toc;
                
                %compute and save the AUC
                [~,~,~,AUC_i]=perfcurve( B( : )>0, B_hat( : ),1);
                AUC( N_i, network_i, model_count ) = AUC_i;
                
                % reset the lastwarn message and id
                lastwarn('', '');
                
                %predict the nodal state with the surrogate network B_hat
                switch model_count
                    case 1
                        x_future_pred = compute_nodal_states_lotka_volterra( x_future( :, 1 ), alpha, theta, B_hat, t_prediction  );
                    case 2
                        x_future_pred = compute_nodal_states_mutualistic_pop( x_future( :, 1 ), alpha, theta, B_hat, t_prediction  );
                    case 3
                        x_future_pred = compute_nodal_states_michaelis_menten( x_future( :, 1 ), parameters.MM.hill_coeff, B_hat, t_prediction  );
                    case 4
                        x_future_pred = compute_nodal_states_SIS( x_future( :, 1 ), delta, B_hat, t_prediction);
                    case 5
                        x_future_pred = compute_nodal_states_kuramoto( x_future( :, 1 ), omega, B_hat, t_prediction);
                    case 6
                        x_future_pred = compute_nodal_states_cowan_wilson( x_future( :, 1 ), parameters.cw.tau, parameters.cw.mu, B_hat, t_prediction);
                end
                
                %save the prediction error
                error_pred(  N_i, 1:length( t_prediction ), network_i, model_count ) = mean( abs( x_future_pred - x_future ) );
                
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

%Fig S11: Prediction error versus N
figure

%plot the number of observations n for which a curve is plotted
n_plot = 25:25:100;

[~, sort_ind ] = sort( N_all );
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
    
    data_plot_tmp = mean(  error_pred( :, n_plot, :, model_count ), 3 );  
    data_plot_tmp = squeeze( data_plot_tmp' );
    data_plot_tmp = data_plot_tmp( :, sort_ind );
    
    subplot( 3, 2, model_count )
    
    plot( N_all( sort_ind ), data_plot_tmp', '-x' )
    
    xlabel('Number of Nodes N')
    ylabel( 'Prediction Error' )
    
    title( model_name )
    
end


%Fig S18: Computation time versus N
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
    
    data_plot_tmp = mean(  time_all( :, :, model_count ), 2 );
    data_plot_tmp = data_plot_tmp( N_all>=N_min );
    data_plot_tmp = data_plot_tmp( sort_ind )';
        
    subplot( 3, 2, model_count )
    p = polyfit( log( N_all ),...
        log( data_plot_tmp ),2);
    
    f1 = polyval(p,log( N_all));
    
    loglog( N_all, data_plot_tmp, '-x' )
    hold
    loglog( N_all, exp( f1 ), '--ro' )
    legend( 'data' ,'fit')
    
    title( strcat( model_name, '; c1 = ', num2str( p( 1 )), '; c2 = ', num2str( p( 2 )) , '; c3 = ', num2str( p( 3 )) ))
    
    xlabel('Number of Nodes N')
    ylabel( 'Runtime [s]' )    
end

end

