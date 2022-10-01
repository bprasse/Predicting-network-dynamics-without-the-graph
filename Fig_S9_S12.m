function Fig_S9_S12( seed, N, num_observations, random_network_model )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%the number of nodes
if nargin < 2
    N=100;
end

%the number of observations
if nargin < 3
    num_observations = 2e2;
end

%the random graph model
if nargin < 4
    random_network_model = 'BA';
end

%the filename to save/read the results
filename = strcat( './results/Fig_S9_S12_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', random_network_model );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    %the number of networks
    num_networks = 1e2;
    
    %set the observation time fraction t_obs/T_max
    T_obs_fraction_all = [ 0.25 0.5 ];
    
    %set the parameters
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
    
    %only for plotting: num_N_plot is the number of nodal traces that are plotted
    num_N_plot = 6;
    
    %the number of different noise levels (standard deviations)
    num_noise = 10;
    
    %the minimum standard deviation
    std_noise_min = 1e-4;
    
    %the maximum standard deviation
    std_noise_max = 1e-2;
    
    %all values for the standard deviation, logarithmically evenly spaced
    std_noise = [0 logspace( log10( std_noise_min ), log10( std_noise_max ), num_noise-1 )];
    
    %pre-allocate the prediction errors
    error_pred = nan( length( T_obs_fraction_all ), num_noise, num_networks, 6  );
        
    %loop over all standard deviations
    for noise_count = 1:num_noise
        std_i = std_noise( noise_count );
        
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
            for T_obs_count = 1:length( T_obs_fraction_all )
                
                T_obs_fraction = T_obs_fraction_all( T_obs_count );
                
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
                    
                    %obtain the maximum past nodal state for each node
                    x_max = max( x, [], 2 );
                    
                    %obtain the maximum future nodal state for each node
                    x_max_future = max( x_future, [], 2 );
                    
                    %add noise to the past nodal state sequence
                    x = x + std_i*( x_max*ones( 1, length( t_obs )) ).*randn( N, length( t_obs ));
                    
                    %add noise to the future nodal state sequence
                    x_future = x_future + std_i*( x_max_future*ones( 1, length( t_prediction )) ).*randn( N, length( t_prediction ));
                    
                    %numerical differentiation
                    dx = transpose( diff( transpose( x ) ) )/delta_T;
                    
                    %delete last nodal state such that x and dx have the same size
                    x = x( :, 1:end-1 );
                    t_obs( end ) = [];
                    
                    %start the network reconstruction
                    if model_count <6
                        [ B_hat, ~ ] = network_reconstruction( x, dx, model_name, input_param, parameters );
                    else
                        [ B_hat, ~ ] = network_reconstruction( x, dx, 'cowan_wilson', input_param, parameters );
                    end
                    
                    % reset the lastwarn message and id, to avoid the potential error
                    % MATLAB:ode45:IntegrationTolNotMet further below
                    lastwarn('', '');
                    
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
                    
                    %read last warning
                    [warnMsg, warnId] = lastwarn();
                    
                    if isempty( warnId )
                        %no error, save the prediction error
                        error_pred(  T_obs_count, noise_count, network_i, model_count ) = mean( mean( abs( x_future_hat - x_future ) ) );
                    else
                        if  strcmp( warnId, 'MATLAB:ode45:IntegrationTolNotMet' )
                            %an error occurred, save 'Inf' as prediction error
                            error_pred(  T_obs_count, noise_count, network_i, model_count ) = Inf;                            
                        else
                            error( strcat( 'Unknown routine for: ', warnId ) )
                        end
                    end
                                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% plot nodal trajectories, corresponding to Fig_S9
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if network_i==1 && noise_count == 1
                        
                        % obtain N_plot, the nodes whose state is plotted
                        N_plot = nan( 1, num_N_plot );
                        x_vector_quantisation = x_future( :, end);
                        x_end_quantised = linspace( min( x_vector_quantisation( :, end ) ), max( x_vector_quantisation( :, end ) ), num_N_plot );
                        
                        
                        %the nodes in N_plot are obtain by finding the nodes in x( :, end )
                        %that are closes to the values in x_end_quantised
                        for curr_N_plot = 1:(num_N_plot+1)
                            [~, indd] = min( abs( x_future( :, end ) - x_end_quantised( curr_N_plot ) ));
                            N_plot( curr_N_plot ) = indd;
                        end
                        
                        %plot the trajectories for the current model
                        subplot( 3, 2, model_count )
                        plot( [t_obs t_prediction], [x( N_plot, : ) x_future( N_plot, : )]', 'b' )
                        title( model_name )
                        
                    end
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
%%% plot results, corresponding to Fig_S12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for T_obs_count = 1:length( T_obs_fraction_all )
    figure
    sgtitle( strcat( 't_{obs}/T_{max} = ', num2str( T_obs_fraction_all( T_obs_count ) )))
    
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
        error_pred_i = squeeze( error_pred( T_obs_count, :, :, model_count ));
        inf_ind = ( error_pred_i == Inf );
        error_pred_i( inf_ind ) = nan;
        num_inf_error_pred = sum( inf_ind( : ) );
        
        subplot( 3, 2, model_count )
        boxplot( transpose( error_pred_i ))
        ylabel( 'Prediction error epsilon; ' )
        if num_inf_error_pred>0
            title( strcat( model_name, '; Num. infinite error: ', num2str( num_inf_error_pred ) ))
        else
            title(  model_name )
        end
    end
    
end

end


