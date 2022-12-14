%%
%Medium Scale New Keynesian Model with Heterogeneous Expectations
%Main file

clear
clc

%% Set seed state
seed = 567;
rng(seed);

%% Characteristics of model

number_endogenous_variables = 24; %Number of endogenous variables
number_exogenous_variables = 7; %Number of exogenous variables
number_aux_variables = 4; %Number of auxiliary variables
number_observed_variables = 7; %Number of observable variables

%Number of variables in state transition equation
number_state_variables = number_endogenous_variables + number_exogenous_variables + number_aux_variables; 

%% Load data set

%Data set used for estimation
%1 = pre great moderation
%2 = great moderation
%3 = pre great recession
%4 = full data set
data_set_identifier = 4;

[data,Sigma_hat,c,first_observation,Sigma_u_sd,theta,prior_information] = msnk_alh_load_data_set(data_set_identifier);

%% Metropolis-Hastings (MH) characteristics

number_draws = 100000; %Number of MH draws
burn_proportion = 0;  %Proportion of draws to discard in MH chain
percentile = 0.1; %Percentile for interval estimates
log_marginal_likelihood_tau = 0.95; %Tuning parameter used in the calculation of the modified harmonic mean

%% Get value of log posterior kernel with the initial parameter values
[log_prior,log_likelihood,log_posterior,hee,A,B,s,error_count] = msnk_alh_log_posterior_calculate(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,data,theta,prior_information,Sigma_u_sd,first_observation);
                                                                                                          
%% Run Metropolis-Hastings algorithm
[mh_theta,mh_theta_log_prior,mh_theta_log_likelihood,mh_theta_log_posterior_kernel,acceptance_rate] = msnk_alh_random_walk_metropolis_hastings_algorithm(theta,Sigma_hat,c,number_draws,number_endogenous_variables,data,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,prior_information,Sigma_u_sd,burn_proportion,first_observation);
[log_posterior_max,index_number] = max(mh_theta_log_posterior_kernel);
theta_max = mh_theta(index_number,:)';

%% Save the Metropolis-Hastings draws
mh_theta_5 = mh_theta;
mh_theta_log_prior_5 = mh_theta_log_prior;
mh_theta_log_likelihood_5 = mh_theta_log_likelihood;
mh_theta_log_posterior_kernel_5 = mh_theta_log_posterior_kernel;

%% Build Metropolis-Hastings chain
% mh_theta = [mh_theta_1;mh_theta_2;mh_theta_3;mh_theta_4;mh_theta_5;mh_theta_6;mh_theta_7;mh_theta_8;mh_theta_9;mh_theta_10];
% mh_theta_log_prior = [mh_theta_log_prior_1;mh_theta_log_prior_2;mh_theta_log_prior_3;mh_theta_log_prior_4;mh_theta_log_prior_5;mh_theta_log_prior_6;mh_theta_log_prior_7;mh_theta_log_prior_8;mh_theta_log_prior_9;mh_theta_log_prior_10];
% mh_theta_log_likelihood = [mh_theta_log_likelihood_1;mh_theta_log_likelihood_2;mh_theta_log_likelihood_3;mh_theta_log_likelihood_4;mh_theta_log_likelihood_5;mh_theta_log_likelihood_6;mh_theta_log_likelihood_7;mh_theta_log_likelihood_8;mh_theta_log_likelihood_9;mh_theta_log_likelihood_10];
% mh_theta_log_posterior_kernel = [mh_theta_log_posterior_kernel_1;mh_theta_log_posterior_kernel_2;mh_theta_log_posterior_kernel_3;mh_theta_log_posterior_kernel_4;mh_theta_log_posterior_kernel_5;mh_theta_log_posterior_kernel_6;mh_theta_log_posterior_kernel_7;mh_theta_log_posterior_kernel_8;mh_theta_log_posterior_kernel_9;mh_theta_log_posterior_kernel_10];

mh_theta = [mh_theta_1;mh_theta_2;mh_theta_3;mh_theta_4;mh_theta_5];
mh_theta_log_prior = [mh_theta_log_prior_1;mh_theta_log_prior_2;mh_theta_log_prior_3;mh_theta_log_prior_4;mh_theta_log_prior_5];
mh_theta_log_likelihood = [mh_theta_log_likelihood_1;mh_theta_log_likelihood_2;mh_theta_log_likelihood_3;mh_theta_log_likelihood_4;mh_theta_log_likelihood_5];
mh_theta_log_posterior_kernel = [mh_theta_log_posterior_kernel_1;mh_theta_log_posterior_kernel_2;mh_theta_log_posterior_kernel_3;mh_theta_log_posterior_kernel_4;mh_theta_log_posterior_kernel_5];

% mh_theta = [mh_theta_1;mh_theta_2;mh_theta_3;mh_theta_4];
% mh_theta_log_prior = [mh_theta_log_prior_1;mh_theta_log_prior_2;mh_theta_log_prior_3;mh_theta_log_prior_4];
% mh_theta_log_likelihood = [mh_theta_log_likelihood_1;mh_theta_log_likelihood_2;mh_theta_log_likelihood_3;mh_theta_log_likelihood_4];
% mh_theta_log_posterior_kernel = [mh_theta_log_posterior_kernel_1;mh_theta_log_posterior_kernel_2;mh_theta_log_posterior_kernel_3;mh_theta_log_posterior_kernel_4];

%% Remove the burn-in draws
theta_post_burn = mh_theta((burn_proportion*number_draws)+1:end,:);
log_prior_post_burn = mh_theta_log_prior((burn_proportion*number_draws)+1:end,:);
log_likelihood_post_burn = mh_theta_log_likelihood((burn_proportion*number_draws)+1:end,:);
log_posterior_post_burn = mh_theta_log_posterior_kernel((burn_proportion*number_draws)+1:end,:);
                                                                                                                                                         
%% Get estimates
[estimates_point,estimates_interval] = bayesian_estimates(theta_post_burn,percentile);

%% Get value of log marginal likelihood using the modified harmonic mean
log_marginal_likelihood = modified_harmonic_mean(theta_post_burn,log_prior_post_burn,log_likelihood_post_burn,log_marginal_likelihood_tau);

%% Create posterior distribution plots
% figure(1)
% msnk_alh_posterior_distributions_learning(theta_post_burn);
% msnk_alh_posterior_distributions(theta_post_burn);

%% Create trace plots
% figure(2)
% msnk_alh_trace_plots(theta_post_burn,log_posterior_post_burn);

%% Create recursive mean plots
% figure(3)
% msnk_alh_recursive_mean_plots(theta_post_burn);

%% Convergence diagnostics
% convergence_diagnostics = gewconv_diagnostic(theta_post_burn);