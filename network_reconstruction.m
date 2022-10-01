function [ B_hat, rho_best ] = network_reconstruction( x, dx, model, model_parameters, parameters )
%Obtain an estimate B_T_est of the infection rate matrix and an estimate delta_T_est of the curing rate vector by the constrained lasso method.
%The penalty parameter rho is set by hold-out cross validation.

N = size( x, 1 );
num_observations = size( x, 2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise parameter estimates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_hat = nan( N );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% partition the data for hold-out validation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_train_max = ceil( num_observations*( 1 - parameters.lasso_options.hold_out_ratio ) );
training_indices = 1:ind_train_max;
validation_indices = ind_train_max+1:num_observations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% estimate parameters for every rho value, choose final paramters by
%%% hold-out validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = nan( N, 1 );
theta = nan( N, 1 );
delta = nan( N, 1 );
omega= nan( N, 1 );
h = nan;
tau = nan;
mu = nan;
f_coeff = nan( N, 4 );
g_coeff = nan( N, 10 );

switch model
    case 'LV'
        alpha = model_parameters.alpha;
        theta = model_parameters.theta;
    case 'MP'
        alpha = model_parameters.alpha;
        theta = model_parameters.theta;
    case 'MM'
        h = parameters.MM.hill_coeff;
    case 'SIS'
        delta = model_parameters.delta;
    case 'kuramoto'
        omega= model_parameters.omega;
    case 'cowan_wilson'
        tau = parameters.cw.tau;
        mu = parameters.cw.mu;
    case 'g_het'
        f_coeff = model_parameters.f_coeff;
        g_coeff = model_parameters.g_coeff;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%redefine some variables for parfor loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_min_ratio = parameters.lasso_options.rho_min_ratio;
rho_max_ratio = parameters.lasso_options.rho_max_ratio;
num_rho = parameters.lasso_options.num_rho;
link_rounding_threshold = parameters.lasso_options.link_rounding_threshold;
quadprog_options = parameters.quadprog_options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimates the links to every node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_best = nan( 1, N );

parfor node = 1:N
    
    V_i = nan;
    F_i = nan;
    
    switch model
        case 'LV'
            [ V_i, F_i ] = linear_system_lotka_volterra( x, dx, alpha( node ), theta( node ), node );
        case 'MP'
            [ V_i, F_i ] = linear_system_mutualistic_pop( x, dx, alpha( node ), theta( node ), node );
        case 'MM'
            [ V_i, F_i ] = linear_system_michaelis_menten( x, dx, h, node );
        case 'SIS'
            [ V_i, F_i ] = linear_system_SIS( x, dx, delta( node ), node );
        case 'kuramoto'
            [ V_i, F_i ] = linear_system_kuramoto( x, dx,  omega( node ), node );
        case 'cowan_wilson'
            [ V_i, F_i ] = linear_system_cowan_wilson( x, dx, tau, mu, node );
        case 'g_het'
            [ V_i, F_i ] = linear_system_g_het( x, dx, f_coeff, g_coeff, node );
            
    end
    
    %training set
    V_i_train = V_i( training_indices );
    F_i_train = F_i( training_indices, : );
    
    %validation set
    V_i_valid = V_i( validation_indices );
    F_i_valid = F_i( validation_indices, : );
    
    %maximum value for the regularisation parameter rho
    rho_max = 2*max( abs( transpose( F_i_train )*V_i_train ));
    
    %logarithmically equally spaced values for rho
    rho_all = rho_max*logspace( log10( rho_min_ratio ), ...
        log10( rho_max_ratio ), num_rho );
    
    best_fit = inf;
    
    for rho = rho_all
        
        %solve the constrained linear least-squares problem for the current
        %value of the regularisation parameter rho on the training set F_i_train, V_i_train
        x_sol =  const_lin_lsq( F_i_train, V_i_train, rho, quadprog_options );
        
        %round small estimates to zero
        x_sol( x_sol < link_rounding_threshold ) = 0;
        
        %error fit_new with respect to the Euclidean norm for the current value of rho
        fit_new = norm( F_i_valid*x_sol - V_i_valid )^2;
        
        %if error fit_new is better than best error best_fit, save rho and
        %fit_new as best value
        if fit_new < best_fit
            best_fit = fit_new;
            rho_best( node ) = rho;
        end
        
        %check if the solution x_sol for current rho is all zero
        all_zero = all( x_sol == 0 );
        
        %if yes, then abort (x_sol will be all zeor for any rho larger than
        %the current rho)
        if all_zero && best_fit < inf
            break
        end
    end
    
    %solve the constrained linear least-squares problem for the best
    %value of the regularisation parameter rho_best on the whole set F_i, V_i
    x_sol =  const_lin_lsq( F_i, V_i, rho_best( node ), quadprog_options );
    
    %round small estimates to zero
    x_sol( x_sol < link_rounding_threshold ) = 0;
    
    %save x_sol as row of the surrogate matrix
    B_hat( node, : ) = transpose( x_sol );
end

end
