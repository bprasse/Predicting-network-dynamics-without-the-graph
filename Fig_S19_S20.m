function Fig_S19_S20( seed, N, num_observations, random_network_model )

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
filename = strcat( './results/Fig_S19_S20_seed_',  num2str( seed ), '_N_', num2str( N ), '_n_', num2str( num_observations ), '_', random_network_model );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    %set the simulation parameters
    parameters = set_parameters;
    
    %set the maximum prediction time
    parameters.T_max = [ 8  ...% LV
        2 ...1  ...% MP
        4  ...% MM
        2 ... 4 ...% SIS
        2 ...%KUR
        4 ]; %CW
    
    %the standard deviation of the initial perturbtation xi_i
    std_init_error = 0.001;
    
    %the number of initial perturbations (for which xi_i is drawn)
    num_error_runs =  1e2;
    
    %the number of networks
    num_networks = 1e2;
    
    %pre-allocate the prediction errors
    error_pred = nan( num_observations, num_networks, 6 );
    
    %pre-allocate the prediction parameter
    Lambda = nan( num_observations, num_networks, 6 );
    
    %loop over all models
    for model_count = 1:6
        
        %loop over all networks
        for network_i = 1:num_networks
            
            %the observation time fraction t_obs/T_max
            T_obs_fraction = 0.05 + 0.45*rand;
            
            %the zero-one (unweighted) network
            A = create_network( N, random_network_model, [], [], [], parameters.m0_BA, parameters.m_BA, true );
            
            %the link weights
            link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
            
            %the weighted network
            B = zeros( N );
            B( A > 0 ) = link_weights;
            
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
                    [ x_future, dx_future ] = compute_nodal_states_lotka_volterra( x( :, end ), alpha, theta, B, t_prediction);
                    
                    %the input parameters for the network
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
                    [ x, ~ ] = compute_nodal_states_michaelis_menten( x_init, parameters.MM.hill_coeff, B, t_obs  );
                    
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
                    [ x, ~] = compute_nodal_states_SIS( x_init, delta, B, t_obs );
                    
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
                    [ x, ~] = compute_nodal_states_kuramoto( x_init, omega, B, t_obs);
                    
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
                    [ x, ~] = compute_nodal_states_cowan_wilson( x_init, parameters.cw.tau, parameters.cw.mu, B, t_obs  );
                    
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
            if model_count <6
                [ B_hat, ~ ] = network_reconstruction( x, dx, model_name, input_param, parameters );
            else
                [ B_hat, ~ ] = network_reconstruction( x, dx, 'cowan_wilson', input_param, parameters );
            end
            
            %predict the nodal state with the surrogate network B_hat
            switch model_count
                case 1
                    x_future_hat = compute_nodal_states_lotka_volterra( x( :, end ), alpha, theta, B_hat, t_prediction  );
                case 2
                    x_future_hat = compute_nodal_states_mutualistic_pop( x( :, end ), alpha, theta, B_hat, t_prediction  );
                case 3
                    x_future_hat = compute_nodal_states_michaelis_menten( x( :, end ), parameters.MM.hill_coeff, B_hat, t_prediction  );
                case 4
                    x_future_hat = compute_nodal_states_SIS( x( :, end ), delta, B_hat, t_prediction);
                case 5
                    x_future_hat = compute_nodal_states_kuramoto( x( :, end ), omega, B_hat, t_prediction);
                case 6
                    x_future_hat = compute_nodal_states_cowan_wilson( x( :, end ), parameters.cw.tau, parameters.cw.mu, B_hat, t_prediction);
            end
            
            %the deviation of the predicted nodal state from the true nodal
            %state versus time
            delta_pred = zeros( size( x_future ));
            delta_pred( :, 1 ) = x_future_hat( :, 1 ) - x_future( :, 1 );
            
            %the prediction error versus time
            error_pred( 1, network_i, model_count ) = mean( abs( delta_pred( :, 1 ) ) );
            
            %compute the deviation and prediction error at every time
            for t=1:size( x_future,2 )-1
                delta_pred( :, t+1 ) = x_future_hat( :, t+1 ) - x_future( :, t+1 );
                error_pred( t+1, network_i, model_count ) = mean( abs( delta_pred( :, t+1 ) ) );
            end
                        
            %the number of agitation modes
            rank_POD = rank( x );
            
            %obtain the agitation modes as the singular vectors of the
            %nodal state trajectory x
            [ singular_vectors, ~, ~] = svds( x, rank_POD  );
            
            %the agitation modes
            Y = singular_vectors( :, 1:rank_POD );
            
            %Y_orth contains as columns all vectors in R^N that are
            %orthogonal to all the every agiation mode (columns in Y)
            Y_orth = null( transpose( Y ) );
            
            %the projection matrices on the column space of Y and Y_orth,
            %respectively
            P_X = Y * transpose( Y );
            P_X_orth = Y_orth * transpose( Y_orth );
            
            %average error of linearisation around trajectory
            rel_error_L2_traj = zeros( 1, size( x_future,2 ));
            
            %loop over num_error_runs perturbations
            for error_j=1:num_error_runs
                %x_traj denotes the linearisation of the nodal state around its trajectory
                x_traj = zeros( size( x_future ));
                
                %the deviation of the linearisation x_traj to the true
                %nodal state x
                delta_traj = zeros( size( x_future ));
                
                %the error of the linearisation around trajectory, for the
                %current run error_j (which determines the initial
                %perturbation xi)
                rel_error_j_L2_traj = zeros( 1, size( x_future,2 ));
                
                %initialise the linearisation around trajectory at time 
                %t=t_obs as the true nodal state x_future( :, 1 ) times a
                %perturbation with stand. deviation std_init_error
                x_traj( :, 1 ) = x_future( :, 1 ).*( 1 + std_init_error*randn( N, 1 ) ) ;
                
                %initialise the deviation and the error of the
                %linearisation
                delta_traj( :, 1 ) =  x_traj( :, 1 ) - x_future( :, 1 ) ;
                rel_error_j_L2_traj( 1 ) = norm( delta_traj( :, 1 ) );
                
                %iterate the linearisation of the trajectory x_traj over time
                for t=1:size( x_future,2 )-1
                    %obtain Jacobian at the point x_future( :, t )
                    F_LTI = linearise( input_param, B, model_name, x_future( :, t ), parameters );
                    
                    %restricted the Jacobian to the subspace spanned by the
                    %agitation modes 
                    F_LTI_restricted = P_X*F_LTI;
                    
                    %update the linearisation, and the deviation and the error of the
                    %linearisation
                    delta_traj( :, t+1 ) = delta_traj( :, t ) +  delta_T*( F_LTI_restricted* delta_traj( :, t ) - P_X_orth*dx_future( :, t ));
                    x_traj( :, t+1 ) = x_future( :, t ) +  delta_traj( :, t+1 );
                    rel_error_j_L2_traj( t+1 ) = norm( delta_traj( :, t+1 )  );
                end
                
                %normalise the error of the linearisation by error at the
                %initial time t_obs
                rel_error_j_L2_traj = rel_error_j_L2_traj/norm(  delta_traj( :, 1 ) );
                
                %add the error rel_error_j_L2_traj of the linearisation of the current
                %perturbation error_j to rel_error_L2_traj
                rel_error_L2_traj = rel_error_L2_traj + rel_error_j_L2_traj;
            end
            
            %obtain the average error of the linearisation (averaged over
            %all num_error_runs perturbations)
            rel_error_L2_traj = rel_error_L2_traj/num_error_runs;
            
            %compute and store the predictability parameter
            Lambda( 1:length( rel_error_L2_traj ), network_i, model_count) =  log( rel_error_L2_traj );
            Lambda( 1, network_i, model_count) = Lambda( 2, network_i, model_count);            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% save files  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save( strcat( filename, '.mat' ) )
    
end

%Figure S19: prediction error versus the sign of the predictability
%parameter Lambda
figure

%pre-allocate the correlation of the predictability coefficient Lambda and
%the prediction error error_pred, for all Lambda>0 and models
rho_diverging_delta_t =  nan( 6, num_networks );

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
    
    %pre-allocate the correlation of the predictability coefficient Lambda and
%the prediction error error_pred, for all Lambda and the current model_count
    rho = nan( 1, num_networks );
    
    %the maximum Lambda across all times t
    max_Lambda = nan( 1, num_networks );
    
    %loop over all networks
    for network_i = 1:num_networks
        
        %find the largest time index for which a prediction was made (since
        %T_obs is generated randomly, the number of predicted time steps
        %varies across different values for network_i)
        ind_max = find( isnan( Lambda( :, network_i, model_count ) ), 1)-1;
        
        %compute and save the correlation of Lambda and log( error_pred )
        cov_matrix = corrcoef(  Lambda( 2:ind_max, network_i, model_count ),log(error_pred(   2:ind_max, network_i, model_count )) );        
        rho( network_i ) = cov_matrix( 1,2 );
        
        %save the maximum value of Lambda across all times for the current network_i
        max_Lambda( network_i ) = max( Lambda( 2:end, network_i, model_count ) ) ;
    end
    
    %find all network_i for which Lambda>=0 at some time t
    ind_diverging_delta_t = find( max_Lambda>=0 );
    
    %find all network_i for which Lambda<0 at all times t
    ind_converging_delta_t = find( max_Lambda<0 );
    
    %save the correlation of Lambda and log( error_pred ) for all diverging
    %runs (given by ind_diverging_delta_t)
    rho_diverging_delta_t( model_count, 1:length( ind_diverging_delta_t ) ) = rho( ind_diverging_delta_t );
    
    %the prediction error at the maximum prediction time 
    error_pred_t_end = nan( 1, num_networks );
    
    %loop over all networks
    for network_i=1:num_networks
        
        %find the largest time index for which a prediction was made (since
        %T_obs is generated randomly, the number of predicted time steps
        %varies across different values for network_i)
        ind_max = find( isnan( error_pred( :, network_i, model_count ) ), 1)-1;
        
        %save the prediction error at the maximum prediction time 
        error_pred_t_end( network_i ) = error_pred(  ind_max, network_i, model_count );
    end
       
    %plot the prediction error at the maximum prediction time, constrained 
    %to of Lambda>=0 versus constrained to Lambda<0
    subplot( 3, 2, model_count )
    hold on
    boxplot( transpose( [ error_pred_t_end(  ind_converging_delta_t ) nan( 1, num_networks-length( ind_converging_delta_t ) );...
        error_pred_t_end( ind_diverging_delta_t )  nan( 1, num_networks-length( ind_diverging_delta_t ) )  ] ), 'Labels',{'Lambda_{max}<0', 'Lambda_{max}>=0'} )
    ylabel( 'Prediction Error $$\epsilon(T_{max})$$', 'Interpreter', 'Latex' )
    title( model_name )
end

%Figure S20: correlation of the prediction error and the predictability
%parameter Lambda, constrained to Lambda>=0
figure
boxplot( transpose( rho_diverging_delta_t ), 'Labels',{'LV','MP','MM','SIS','KUR','WC'} )
ylabel( 'Correlation $$\Lambda(t), \log(\epsilon(t))$$', 'Interpreter', 'Latex' )

end



function F_LTI = linearise( input_param, B, model, x_0, parameters )
%obtain the Jacobian at the point x_0

N = length( x_0 );

B_diag = diag( B( eye( N )>0 ) );
B_no_diag = B - B_diag;

switch model
    case 'LV'
        f_i_jacobian = diag( input_param.alpha - 2*input_param.theta.*x_0 );
        
        g_1st_arg_jacobian = -ones( N, 1 )*transpose( x_0 );
        g_1st_arg_jacobian( eye( N )>0 ) = 0;
        
        g_2nd_arg_jacobian = -x_0*ones( 1, N );
        g_2nd_arg_jacobian( eye( N )>0 )= -2*x_0;        
    case 'MP'
        
        f_i_jacobian = diag( input_param.alpha - 2*input_param.theta.*x_0 );
        
        g_1st_arg_jacobian = ones( N, 1 )*transpose( x_0.^2./( 1 + x_0.^2 ) );
        g_1st_arg_jacobian( eye( N )>0 ) = 0;
        
        g_2nd_arg_jacobian = x_0*transpose( 2*x_0./( 1 + x_0.^2 ).^2  );
        g_2nd_arg_jacobian( eye( N )>0 )= x_0.^2.*(x_0.^2+3)./( 1 + x_0.^2 ).^2;        
    case 'MM'
        if parameters.MM.hill_coeff ~= 2
            error( 'Derivative is only correct if the Hill coefficient equals h=2!' )
        end
        f_i_jacobian = -eye( N );
        
        g_1st_arg_jacobian = zeros( N );
        
        g_2nd_arg_jacobian = ones( N, 1 )*transpose( 2*x_0./( 1 + x_0.^2 ).^2  );
        
    case 'SIS'
        f_i_jacobian = diag( -input_param.delta  );
        
        g_1st_arg_jacobian = -ones( N, 1 )*transpose( x_0 );
        g_1st_arg_jacobian( eye( N )>0 ) = 0;
        
        g_2nd_arg_jacobian = (1-x_0)*ones( 1, N );
        g_2nd_arg_jacobian( eye( N )>0 )= 1-2*x_0;        
    case 'kuramoto'
        f_i_jacobian = zeros( N );
        
        g_1st_arg_jacobian =  cos( x_0*ones( 1, N ) -ones( N, 1 )*transpose( x_0 ));
        g_1st_arg_jacobian( eye( N )>0 ) = 0;
        
        g_2nd_arg_jacobian = -cos( x_0*ones( 1, N ) -ones( N, 1 )*transpose( x_0 ));
        g_2nd_arg_jacobian( eye( N )>0 )= 0;
        
    case 'cw'
        f_i_jacobian = -eye( N );
        
        g_1st_arg_jacobian = zeros( N );
        
        g_2nd_arg_jacobian = ones( N, 1 )*transpose( parameters.cw.tau*exp( -parameters.cw.tau*(x_0 - parameters.cw.mu ) )./( 1 + exp( -parameters.cw.tau*(x_0 - parameters.cw.mu ) )).^2 );
   otherwise
        error( 'Unknown model' )
end

F_LTI = f_i_jacobian + diag( sum( B_no_diag.*g_1st_arg_jacobian, 2 )) + B.*g_2nd_arg_jacobian;

end




