function parameters = set_parameters

%%%%%%%%%%%%%%%%%%%%%%
%%%Network parameters
%%%%%%%%%%%%%%%%%%%%%%
parameters.plot = true; 

%random graph model 
parameters.random_network_model = 'BA';  

%parameter for BA random graphs
parameters.m0_BA = 3; 
parameters.m_BA = 3; 

%minimum and maximum link weight
parameters.min_link_weight = 0.5;
parameters.max_link_weight = 1.5;

%%%%%%%%%%%%%%%%%%%%%%
%%% Lotka-Volterra 
%%%%%%%%%%%%%%%%%%%%%%
parameters.LV.sigma_alpha = 0.5;
parameters.LV.sigma_theta = 0.5;
parameters.LV.x_init_max = 1;

%%%%%%%%%%%%%%%%%%%%%%
%%%Mutualistic population
%%%%%%%%%%%%%%%%%%%%%%
parameters.MP.sigma_alpha = 0.5;
parameters.MP.sigma_theta = 0.5;
parameters.MP.x_init_max = 20;

%%%%%%%%%%%%%%%%%%%%%%
%%% Michaelis-Menten parameters
%%%%%%%%%%%%%%%%%%%%%%
parameters.MM.hill_coeff = 2;
parameters.MM.x_init_max = 2;

%%%%%%%%%%%%%%%%%%%%%%
%%% SIS parameters
%%%%%%%%%%%%%%%%%%%%%%
parameters.SIS.R_0_SIS = 1.5;
parameters.SIS.sigma_delta = 0.5;
parameters.SIS.x_init_max = 0.1;

%%%%%%%%%%%%%%%%%%%%%%
%%%Kuramoto parameters
%%%%%%%%%%%%%%%%%%%%%%
parameters.kuramoto.sigma_omega = 0.1*pi;
parameters.kuramoto.x_init_min = -pi/4;
parameters.kuramoto.x_init_max = pi/4;

%%%%%%%%%%%%%%%%%%%%%%
%%%Cowan-Wilson parameters
%%%%%%%%%%%%%%%%%%%%%%
parameters.cw.tau = 1;
parameters.cw.mu = 1;
parameters.cw.x_init_max = 10;

%%%%%%%%%%%%%%%%%%%%%% 
%%%Maximum time window
%%%%%%%%%%%%%%%%%%%%%%
parameters.T_max = [5 0.025 3 0.5 1 4]; 

%%%%%%%%%%%%%%%%%%%%%% 
%%%Solver parameters for obtaining the surrogate network
%%%%%%%%%%%%%%%%%%%%%%
parameters.quadprog_options = optimoptions('quadprog', 'Display', 'off' );
parameters.quadprog_options.OptimalityTolerance = 1e-8;    
parameters.quadprog_options.MaxIterations = 2e2;           

parameters.lasso_options.rho_min_ratio  = 1e-6;
parameters.lasso_options.rho_max_ratio  = 1e-2;
parameters.lasso_options.link_rounding_threshold = 1e-2;

parameters.lasso_options.num_rho = 20;
parameters.lasso_options.hold_out_ratio  = 0.2; 

end

