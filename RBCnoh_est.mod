% Model from Noh (2018)
% Baseline model
% WITH investment
% =========================================================================
% Estimation
% -------------------------------------------------------------------------

% DATASET
data(file='temp_stoch_simul_obs.mat');
%varobs C_obs I_obs H_obs;
%---------------------------------
% measurement errors
% 0 - Include,  1 - Exclude
@#ifndef NO_ME 
@#define NO_ME = 0 
@#endif

 % 0 - large, 1 - small
@#ifndef MAGNITUDE
@#define MAGNITUDE=0
@#endif
%-------------------------------------------------

@#ifndef LINEAR_KALMAN
@#define LINEAR_KALMAN = 0
@#endif

@#ifndef NON_LINEAR_KALMAN
@#define NON_LINEAR_KALMAN = 0
@#endif

@#ifndef ALGO_SIR
@#define ALGO_SIR = 1
@#endif

% Kollmann 2013
@#ifndef GAUSSIAN_APPROXIMATION 
@#define GAUSSIAN_APPROXIMATION = 0
@#endif

% Meyer-Gohde 2013
@#ifndef NON_CENTRAL_MEAN 
@#define NON_CENTRAL_MEAN = 0
@#endif

% Meyer-Gohde 2013
@#ifndef NON_CENTRAL_STEADY 
@#define NON_CENTRAL_STEADY = 0
@#endif    

% Holden 2018
@#ifndef SECOND_ORDER_EXTENDED_KALMAN 
@#define SECOND_ORDER_EXTENDED_KALMAN = 0
@#endif

% Holden 2018 - pruning
@#ifndef SECOND_ORDER_EXTENDED_KALMAN_P 
@#define SECOND_ORDER_EXTENDED_KALMAN_P = 0
@#endif

var C H I Z mu K Y C_obs H_obs I_obs ;
varobs C_obs I_obs H_obs;

varexo epsZ epsM; 

parameters sigma alpha Hss betta  delta rhoZ rhoM psi; 

alpha=0.67;
betta=0.99;
delta=0.025;
sigma=2.75;
rhoZ=0.95;
rhoM=0.72;

Hss=0.3;

model; 
    psi*(1-H)^(-sigma)=(1/C)*(alpha*Z*(K(-1)/H)^(1-alpha)); 
    C(+1)/(betta*C*mu)=(1-alpha)*Z(+1)*(K/H(+1))^(-alpha)+((1-delta)/mu); 
    K-((1-delta)*K(-1))=mu*I; 
    I=Y-C;
    Y=Z*(K(-1)^(1-alpha))*(H^alpha); 
    log(Z)=rhoZ*log(Z(-1))+epsZ; 
    log(mu)=rhoM*log(mu(-1))+epsM;
    C_obs=C;
    I_obs=I;
    H_obs=H;
end;

steady_state_model;
Z=1;
mu=1;
H=Hss;
K=H*(((1/(mu*betta))-((1-delta)/mu))*(1/((1-alpha)*Z)))^(-(1/alpha));
Y=Z*(K^(1-alpha))*(H^alpha);
I=(delta*K)/mu;
C=Y-I;
psi=(alpha/C)*(Y/H)*(1-H)^sigma;
C_obs=C;
I_obs=I;
H_obs=H;
end;

% =========================================================================
% Steady-State, Checks and Diagnostics
% =========================================================================
steady;                 % compute steady state given the starting values
resid;                  % check the residuals of model equations evaluated at steady state
check;                  % check Blanchard-Kahn conditions
model_diagnostics;      % check obvious model errors
options_.debug=1;


shocks;
var epsZ; stderr 0.007;
var epsM; stderr 0.06;
end;

% =========================================================================
% Specify Economic Priors
% =========================================================================
estimated_params;
% --------------------------------------------------------------------------------------------------------------------------------------------------
%PARAMETER_NAME, INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND, PRIOR_SHAPE,   PRIOR_MEAN, PRIOR_STANDARD_ERROR, PRIOR_3RD_PARAMETER, PRIOR_4TH_PARAMETER;
% --------------------------------------------------------------------------------------------------------------------------------------------------
%Hss,                       Hss,            ,            , uniform_pdf,             ,                     ,                  0.1,                  0.6;
alpha,                   alpha,            ,            , uniform_pdf,             ,                     ,              0.3,                  1;
betta,                   betta,            ,            , uniform_pdf,             ,                     ,                 0.75,                  0.995;
delta,                   delta,            ,            , uniform_pdf,             ,                     ,              0.01,               0.05;
sigma,                   sigma,            ,            , uniform_pdf,             ,                     ,              0.00001,                100;
rhoZ,                     rhoZ,            ,            , uniform_pdf,             ,                     ,              0.00001,                  1;
rhoM,                     rhoM,            ,            , uniform_pdf,             ,                     ,              0.00001,                  1;
stderr epsM,              0.06,            ,            , uniform_pdf,             ,                     ,              0.00001,                100;
stderr epsZ,              0.007,            ,            , uniform_pdf,             ,                     ,              0.00001,                100;

@#if NO_ME
% no measurement errors
@#else
    @#if MAGNITUDE
          %small
    stderr C_obs, 0.004,,, uniform_pdf,,, 0.00001, 100;
    stderr I_obs, 0.019,,, uniform_pdf,,, 0.00001, 100;
    stderr H_obs, 0.009,,, uniform_pdf,,, 0.00001, 100;
    @#endif
    %large
     stderr C_obs, 0.004,,, uniform_pdf,,, 0.00001, 100;
    stderr I_obs, 0.019,,, uniform_pdf,,, 0.00001, 100;
    stderr H_obs, 0.009,,, uniform_pdf,,, 0.00001, 100; 
  
    
 @#endif
    end;



@#if LINEAR_KALMAN
  estimation( order=1
      ,mode_compute=9, mh_replic=0, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 ) 
      %, prior_trunc=0
      , nograph, mode_check, graph_format=fig) C_obs I_obs H_obs;
@#endif

@#if NON_LINEAR_KALMAN
  	estimation( order=2,
        filter_algorithm=nlkf, proposal_approximation=cubature
        , mode_compute = 9, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 1000000, 'MaxIter', 1000000, 'TolCon', 1.49011611938477e-08 ),mh_replic=0
        , nograph, mode_check, graph_format=fig) C_obs I_obs H_obs;
@#endif

@#if ALGO_SIR
    estimation( order=2
        , number_of_particles=300
        , resampling=systematic
        , mode_compute = 9, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 1000000, 'MaxIter', 1000000, 'TolCon', 1.49011611938477e-08 ), mh_replic=0
        %, prior_trunc=0,
        , nograph
        , mode_check, graph_format=fig) C_obs I_obs H_obs;
    @#endif

% Kollmann 2015 
@#if GAUSSIAN_APPROXIMATION
    
        %stoch_simul( order = 3, pruning, irf = 0, periods = 0 , nograph);
    
    options_.gaussian_approximation = 1;
    options_.underlying_order = 2; 			// make sure estimation uses a 2rd order underlying approximation
    
    options_.cova_compute = 1;              // disable computing the hessian

estimation(lik_init = 1
    , mode_compute = 9, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 ),mh_replic=0
    %, prior_trunc=0
    , nograph, mode_check, graph_format=fig) C_obs I_obs H_obs;
@#endif

% Meyer-Gohde 2013
@#if NON_CENTRAL_MEAN 
    options_.non_central_approximation = 2; // 2 approximates around the mean, 1 approximates around the RSS
    options_.underlying_order = 2; 			// make sure estimation uses a 3rd order underlying approximation
    // disable computing the hessian
    options_.cova_compute = 1;

stoch_simul( order = 3, pruning, irf = 0, periods = 0 , nograph);

estimation( lik_init = 1
    , mode_compute = 9, mh_replic=0, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 )
    %, prior_trunc=0
    , nograph, mode_check, graph_format=fig) C_obs I_obs H_obs;
@#endif

% Meyer-Gohde 2013
@#if NON_CENTRAL_STEADY 
    options_.non_central_approximation = 1; // 2 approximates around the mean, 1 approximates around the RSS
    options_.underlying_order = 2; 			// make sure estimation uses a 3rd order underlying approximation
// disable computing the hessian
    options_.cova_compute = 1;

stoch_simul( order = 3, pruning, irf = 0, periods = 0, nograph );

estimation( lik_init = 1
    , mode_compute = 9, optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 ), mh_replic=0
    %, prior_trunc=0
    , mode_check, nograph, graph_format=fig) C_obs I_obs H_obs;

@#endif


% Holden 2018 - no pruning 
@#if SECOND_ORDER_EXTENDED_KALMAN 

    options_.extended_kalman_filter = 1; % - 0  off - 1 on    
    options_.add_empty_presamples=0;
    // disable computing the hessian
    options_.cova_compute = 1;
    
stoch_simul( order = 2, irf = 0, periods = 0, nograph, pruning );

options_.pruning=0;

estimation(lik_init = 1
    , order=2, mode_compute = 9,  optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 ), mh_replic=0
    , mode_check, nograph, graph_format=fig) C_obs I_obs H_obs;

@#endif


% Holden 2018 - pruning 
@#if SECOND_ORDER_EXTENDED_KALMAN_P 

    options_.extended_kalman_filter = 1; % - 0  off - 1 on    
    options_.add_empty_presamples=0;
    // disable computing the hessian
    options_.cova_compute = 1;
    options_.pruning=1;

stoch_simul( order = 2, pruning, irf = 0, periods = 0, nograph );

estimation(lik_init = 1
, order=2, mode_compute = 9,  optim = ( 'UseParallel', 'always', 'Display', 'iter', 'MaxFunEvals', 100000, 'MaxIter', 100000, 'TolCon', 1.49011611938477e-08 ), mh_replic=0
, mode_check, nograph, graph_format=fig) C_obs I_obs H_obs;
@#endif

     oo_.rodolfo_M=M_.params;
     oo_.rodolfo_S=M_.Sigma_e;
% 
         set_param_value('alpha',oo_.posterior_mode.parameters.alpha);
         set_param_value('betta',oo_.posterior_mode.parameters.betta);
         set_param_value('delta',oo_.posterior_mode.parameters.delta);
         set_param_value('sigma',oo_.posterior_mode.parameters.sigma);
         set_param_value('rhoZ',oo_.posterior_mode.parameters.rhoZ);
         set_param_value('rhoM',oo_.posterior_mode.parameters.rhoM);
%     
        M_.Sigma_e(1,1)=oo_.posterior_mode.shocks_std.epsM.^2;
        M_.Sigma_e(2,2)=oo_.posterior_mode.shocks_std.epsZ.^2;

%        set_Sigma_e_value('Sigma_e',((oo_.posterior_mode.shocks_std.epsM).^2))
% 
%     oo_.rodolfo_set_matrix=get_posterior_parameters('mode',M_,estim_params_,oo_,options_);

stoch_simul(order=2, periods=600, nograph);
