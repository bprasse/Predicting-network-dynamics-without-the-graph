function Fig_S15_S16_S17( seed, N, num_observations, random_network_model )

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
filename = strcat( './results/Fig_S15_S16_S17_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', random_network_model );

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
    
    %the observation time fraction t_obs/T_max
    T_obs_fraction = 0.25;
    
    %the maximum number of traces S
    num_traces = 20;
    
    %pre-allocate the network reconstruction accuracy AUC
    AUC = nan( num_traces, num_networks, 6 );
    
    %pre-allocate the prediction errors
    error_pred = nan( num_traces, num_networks, 6 );
    
    %loop over all networks
    for network_i = 1:num_networks
        
        %the zero-one (unweighted) network
        A = create_network( N, random_network_model, [], [], [],  parameters.m0_BA,  parameters.m_BA, true );
        
        %the link weights
        link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
        
        %the weighted network
        B = zeros( N );
        B( A > 0 ) = link_weights;
        
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
            
            %x_all is the concatenation of the nodal state of all traces
            x_all = [];
            
            %x_all is the concatenation of the (approximate) derivative of the nodal state of all traces
            dx_all = [];
            
            %loop over num_traces traces
            for trace_i = 1:num_traces
                
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
                t_obs( end ) = [];
                
                %concatenate the traces
                x_all = [ x_all x ];
                dx_all = [ dx_all dx ];
                
                %start the network reconstruction, based on the
                %*concatenated* traces x_all
                if model_count <6
                    [ B_hat, ~ ] = network_reconstruction( x_all, dx_all, model_name, input_param, parameters );
                else
                    [ B_hat, ~ ] = network_reconstruction( x_all, dx_all, 'cowan_wilson', input_param, parameters );
                end
                
                %compute and save the AUC
                [~,~,~,AUC_i]=perfcurve( B( : )>0, B_hat( : ),1);
                AUC(  trace_i, network_i, model_count ) = AUC_i;
                
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
                error_pred(  trace_i, network_i, model_count ) = mean( mean( abs( x_future_hat - x_future ) ) );
                
            end
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

%Fig S15: Prediction error versus the number of traces
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
    boxplot( transpose( reshape( error_pred( :, :, model_count ), num_traces, num_networks ) ))
    xlabel('Number of Trajectories')
    ylabel( 'Prediction Error $$\bar{\epsilon}$$', 'Interpreter', 'Latex' )
    
    title( model_name )
end

%Fig S16: Similarity of surrogate and true network versus the number of traces
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
    boxplot( transpose( reshape( AUC( :, :, model_count ), num_traces, num_networks ) ))
    xlabel('Number of Trajectories')
    ylabel( 'AUC' )
    title( model_name )
end

%Fig S17: Prediction error versus network similarity, across all number of
%trajectories
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
    
    error_all = error_pred( :, :, model_count );
    AUC_all =  AUC( :, :, model_count );
    
    subplot( 3, 2, model_count )
    scatter( AUC_all(:), error_all(:))
    xlabel('AUC')
    ylabel( 'Prediction Error $$\bar{\epsilon}$$', 'Interpreter', 'Latex' )
    title( model_name )
    
end

end


