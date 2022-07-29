%%
%Medium-scale New Keynesian Model with Heterogeneous Expectations
%Simulates the model 

clear
clc

%% Characteristics of model

number_endogenous_variables = 24; %Number of endogenous variables
number_exogenous_variables = 7; %Number of exogenous variables
number_aux_variables = 4; %Number of auxiliary variables
number_observed_variables = 7; %Number of observable variables

%Number of variables in state transition equation
number_state_variables = number_endogenous_variables + number_exogenous_variables + number_aux_variables;

%Load the parameter estimates
parameter_info = load('msnk_al_parameter_estimates.csv');

%Number of periods to simulate the model
T = 343;

%Number of initial periods per simulation to discard
discard = 100;

%Number of simulations to run
number_simulations = 1000;

%Number of autocorrelation coefficient lags to calculate
acor_lags = 10;

%Number of cross-correlation coefficient lags to calculate
ccor_lags = 10;

% Set seed state
seed = 12;
rng(seed);

%% Create storage for summary statistics

%Real GDP growth
dy_mean = zeros(number_simulations,1);
dy_var = zeros(number_simulations,1);
dy_sd = zeros(number_simulations,1);
dy_acf = zeros(acor_lags,number_simulations);

%Consumption growth
dc_mean = zeros(number_simulations,1);
dc_var = zeros(number_simulations,1);
dc_sd = zeros(number_simulations,1);
dc_acf = zeros(acor_lags,number_simulations);

%Investment growth
dinve_mean = zeros(number_simulations,1);
dinve_var = zeros(number_simulations,1);
dinve_sd = zeros(number_simulations,1);
dinve_acf = zeros(acor_lags,number_simulations);

%Real wage growth
dw_mean = zeros(number_simulations,1);
dw_var = zeros(number_simulations,1);
dw_sd = zeros(number_simulations,1);
dw_acf = zeros(acor_lags,number_simulations);

% Hours worked
labobs_mean = zeros(number_simulations,1);
labobs_var = zeros(number_simulations,1);
labobs_sd = zeros(number_simulations,1);
labobs_acf = zeros(acor_lags,number_simulations);

%Inflation
pinfobs_mean = zeros(number_simulations,1);
pinfobs_var = zeros(number_simulations,1);
pinfobs_sd = zeros(number_simulations,1);
pinfobs_acf = zeros(acor_lags,number_simulations);

%Nominal interest rate
robs_mean = zeros(number_simulations,1);
robs_var = zeros(number_simulations,1);
robs_sd = zeros(number_simulations,1);
robs_acf = zeros(acor_lags,number_simulations);

%Relative standard deviations (to dy)
rel_dy = zeros(number_simulations,1);
rel_dc = zeros(number_simulations,1);
rel_dinve = zeros(number_simulations,1);
rel_dw = zeros(number_simulations,1);
rel_labobs = zeros(number_simulations,1);
rel_pinfobs = zeros(number_simulations,1);
rel_robs = zeros(number_simulations,1);

%Correlation matrix
correlation_matrix = zeros(number_observed_variables,number_observed_variables,number_simulations);

%Output growth cross-correlations
ccor_dy_dc = zeros(ccor_lags,number_simulations);
ccor_dy_dinve = zeros(ccor_lags,number_simulations);
ccor_dy_dw = zeros(ccor_lags,number_simulations);
ccor_dy_labobs = zeros(ccor_lags,number_simulations);
ccor_dy_pinfobs = zeros(ccor_lags,number_simulations);
ccor_dy_robs = zeros(ccor_lags,number_simulations);

%Consumption growth cross-correlations
ccor_dc_dy = zeros(ccor_lags,number_simulations);
ccor_dc_dinve = zeros(ccor_lags,number_simulations);
ccor_dc_dw = zeros(ccor_lags,number_simulations);
ccor_dc_labobs = zeros(ccor_lags,number_simulations);
ccor_dc_pinfobs = zeros(ccor_lags,number_simulations);
ccor_dc_robs = zeros(ccor_lags,number_simulations);

%Investment growth cross-correlations
ccor_dinve_dy = zeros(ccor_lags,number_simulations);
ccor_dinve_dc = zeros(ccor_lags,number_simulations);
ccor_dinve_dw = zeros(ccor_lags,number_simulations);
ccor_dinve_labobs = zeros(ccor_lags,number_simulations);
ccor_dinve_pinfobs = zeros(ccor_lags,number_simulations);
ccor_dinve_robs = zeros(ccor_lags,number_simulations);

%Real wage growth cross-correlations
ccor_dw_dy = zeros(ccor_lags,number_simulations);
ccor_dw_dc = zeros(ccor_lags,number_simulations);
ccor_dw_dinve = zeros(ccor_lags,number_simulations);
ccor_dw_labobs = zeros(ccor_lags,number_simulations);
ccor_dw_pinfobs = zeros(ccor_lags,number_simulations);
ccor_dw_robs = zeros(ccor_lags,number_simulations);

%Hours worked cross-correlations
ccor_labobs_dy = zeros(ccor_lags,number_simulations);
ccor_labobs_dc = zeros(ccor_lags,number_simulations);
ccor_labobs_dinve = zeros(ccor_lags,number_simulations);
ccor_labobs_dw = zeros(ccor_lags,number_simulations);
ccor_labobs_pinfobs = zeros(ccor_lags,number_simulations);
ccor_labobs_robs = zeros(ccor_lags,number_simulations);

%Inflation cross-correlations
ccor_pinfobs_dy = zeros(ccor_lags,number_simulations);
ccor_pinfobs_dc = zeros(ccor_lags,number_simulations);
ccor_pinfobs_dinve = zeros(ccor_lags,number_simulations);
ccor_pinfobs_dw = zeros(ccor_lags,number_simulations);
ccor_pinfobs_labobs = zeros(ccor_lags,number_simulations);
ccor_pinfobs_robs = zeros(ccor_lags,number_simulations);

%Nominal interest rate cross-correlations
ccor_robs_dy = zeros(ccor_lags,number_simulations);
ccor_robs_dc = zeros(ccor_lags,number_simulations);
ccor_robs_dinve = zeros(ccor_lags,number_simulations);
ccor_robs_dw = zeros(ccor_lags,number_simulations);
ccor_robs_labobs = zeros(ccor_lags,number_simulations);
ccor_robs_pinfobs = zeros(ccor_lags,number_simulations);

%Errors
errors = zeros(number_simulations,3);

%% Parameters

%Fixed Parameters
ctou = 0.025;
cg = 0.18;
clandaw = 1.5;
curvp = 10;
curvw = 10;

% Estimated parameters
sig_ea = parameter_info(1);
sig_eb = parameter_info(2);
sig_eg = parameter_info(3);
sig_eqs = parameter_info(4);
sig_em = parameter_info(5);
sig_epinf = parameter_info(6);
sig_ew = parameter_info(7);
csadjcost = parameter_info(8);
csigma = parameter_info(9);
chabb = parameter_info(10);
cprobw = parameter_info(11);
csigl = parameter_info(12);
cprobp = parameter_info(13);
cindw = parameter_info(14);
cindp = parameter_info(15);
czcap = parameter_info(16);
cfc = parameter_info(17);
crpi = parameter_info(18);
crr = parameter_info(19);
cry = parameter_info(20);
crdy = parameter_info(21);
constepinf = parameter_info(22);
constebeta = parameter_info(23);
constelab = parameter_info(24);
ctrend = parameter_info(25);
calfa = parameter_info(26);
crhoa = parameter_info(27);
crhob = parameter_info(28);
crhog = parameter_info(29);
crhoqs = parameter_info(30);
crhoms = parameter_info(31);
crhopinf = parameter_info(32);
crhow = parameter_info(33);
gna = parameter_info(34);
gnb = parameter_info(35);
omegaa = parameter_info(36);

theta = [sig_ea;sig_eb;sig_eg;sig_eqs;sig_em;sig_epinf;sig_ew;csadjcost;csigma;chabb;cprobw;csigl;cprobp;cindw;cindp;czcap;cfc;crpi;crr;cry;crdy;constepinf;constebeta;constelab;ctrend;calfa;crhoa;crhob;crhog;crhoqs;crhoms;crhopinf;crhow;gna;gnb;omegaa];

%Composite parameters
cpie = 1 + (constepinf/100);
cbeta = 1/(1+constebeta/100);
cgamma = 1 + ctrend/100;
clandap = cfc;
cbetabar = cbeta*cgamma^(-csigma);
cr = cpie/(cbeta*cgamma^(-csigma));
crk = (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar = (1-(1-ctou)/cgamma);
cik = (1-(1-ctou)/cgamma)*cgamma;
clk = ((1-calfa)/calfa)*(crk/cw);
cky = cfc*(clk)^(calfa-1);
ciy = cik*cky;
ccy = 1 - cg - cik*cky;
crkky = crk*cky;
cwhlc =  (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
conster = (cr-1)*100;
omegab = 1-omegaa;

%% Find the heterogeneous expectations equilibrium and build the state-space matrices that aren't a function of agent beliefs

[H_1,hee,K_1,K_2,K_3,K_4,Sigma_epsilon,t,Psi_0,Psi_1,Psi_2,S_c] = msnk_alh_build_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,theta);

%% Start the simulation replications

for i = 1:number_simulations

    % Initialize the adaptive learning algorithm
    [s,A,B,R_A,R_B,~,~,error_count] = msnk_alh_initialize_learning(number_endogenous_variables,number_exogenous_variables,number_state_variables,T+2,hee);
    
    % Initialize the State-space matrices
    S_1 = zeros(number_state_variables,number_state_variables,T+2);
    S_2 = zeros(number_state_variables,number_exogenous_variables,T+2);
    
    %Initialize storage for individual state variables
    y = zeros(1,T+2);
    c = zeros(1,T+2);
    inve = zeros(1,T+2);
    w = zeros(1,T+2);
    lab = zeros(1,T+2);
    pinf = zeros(1,T+2);
    r = zeros(1,T+2);

    %Initialize storage for individual observational equivalent variables
    dy = zeros(1,T+2);
    dc = zeros(1,T+2);
    dinve = zeros(1,T+2);
    dw = zeros(1,T+2);
    labobs = zeros(1,T+2);
    pinfobs = zeros(1,T+2);
    robs = zeros(1,T+2);
    
    %Generate the shocks
    epsilon = normrnd(0,1,number_exogenous_variables,T+2);

    %Store the initial values of the S_1 and S_2 matrices by updating the adaptive learning state-space matrices that are functions of agent beliefs with the given values of the parameters and the initial values of the agent beliefs
    [S_1(:,:,2),S_2(:,:,2),~,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2),hee,A(:,:,1),B(:,:,1),R_A(:,:,1),R_B(:,:,1),S_1(:,:,1),S_2(:,:,1));
    
    %Generate the initial values of the state variables
    s(:,2) = S_1(:,:,2)*s(:,1) + S_2(:,:,2)*epsilon(:,2);
    
    %Simulate the model
    for j = 3:T+2   
        
        %Get the last two periods of the state variables
        s_last_two_periods = [s(:,j-2),s(:,j-1)];

        %Now update the agent beliefs and the moment matrices (the dependent variables are last period's state variables)
        [A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_beliefs(theta,s_last_two_periods,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1));

        %Now use the updated beliefs to update the adaptive learning state-space matrices that are functions of agent beliefs
        [S_1(:,:,j),S_2(:,:,j),test,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j),hee,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1),S_1(:,:,j-1),S_2(:,:,j-1));
      
        %Now update the values of the state variables to reflect the updated S_1 and S_2 matrix
        s(:,j) = S_1(:,:,j)*s(:,j-1) + S_2(:,:,j)*epsilon(:,j);

        %Get relevant variables
        y(j) = s(19,j); %Output
        c(j) = s(18,j); %Consumption
        inve(j) = s(16,j); %Investment
        w(j) = s(22,j); %Real wage
        lab(j) = s(20,j); %Hours worked
        pinf(j) = s(21,j); %Inflation
        r(j) = s(23,j); %Nominal interest rate

        %Calculate observational equivalent variables
        dy(j) = y(j) - y(j-1) + ctrend;
        dc(j) = c(j) - c(j-1) + ctrend;
        dinve(j) = inve(j) - inve(j-1) + ctrend;
        dw(j) = w(j) - w(j-1) + ctrend;
        labobs(j) = lab(j) + constelab;
        pinfobs(j) = pinf(j) + constepinf;
        robs(j) = r(j) + conster;
        
        %Record the errors for this iteration
        error_count(j,:) = test;
        
    end

    %Drop the initial values
    dy = (dy(discard+1:end))'; %Output growth
    dc = (dc(discard+1:end))'; %Consumption growth
    dinve = (dinve(discard+1:end))'; %Investment growth
    dw = (dw(discard+1:end))'; %Real wage growth
    labobs = (labobs(discard+1:end))'; %Labor hours
    pinfobs = (pinfobs(discard+1:end))'; %Inflation
    robs = (robs(discard+1:end))'; %Nominal interest rate

    %Store the simulation data
    data = [dy,dc,dinve,dw,labobs,pinfobs,robs]; 

    %% Calculate statistics for this individual simultation run

    %Calculate statistics for individual variables
    
    %Real GDP growth
    dy_mean(i) = mean(dy);
    dy_var(i) = var(dy);
    dy_sd(i) = sqrt(dy_var(i));
    dy_acf(:,i) = acf(dy,10);

    %Consumption growth
    dc_mean(i) = mean(dc);
    dc_var(i) = var(dc);
    dc_sd(i) = sqrt(dc_var(i));
    dc_acf(:,i) = acf(dc,10);

    %Investment growth
    dinve_mean(i) = mean(dinve);
    dinve_var(i) = var(dinve);
    dinve_sd(i) = sqrt(dinve_var(i));
    dinve_acf(:,i) = acf(dinve,10);

    %Real wage growth
    dw_mean(i) = mean(dw);
    dw_var(i) = var(dw);
    dw_sd(i) = sqrt(dw_var(i));
    dw_acf(:,i) = acf(dw,10);

    %Hours worked
    labobs_mean(i) = mean(labobs);
    labobs_var(i) = var(labobs);
    labobs_sd(i) = sqrt(labobs_var(i));
    labobs_acf(:,i) = acf(labobs,10);

    %Inflation
    pinfobs_mean(i) = mean(pinfobs);
    pinfobs_var(i) = var(pinfobs);
    pinfobs_sd(i) = sqrt(pinfobs_var(i));
    pinfobs_acf(:,i) = acf(pinfobs,10);

    %Nominal interest rate
    robs_mean(i) = mean(robs);
    robs_var(i) = var(robs);
    robs_sd(i) = sqrt(robs_var(i));
    robs_acf(:,i) = acf(robs,10);

    %Calculate relative volatilities
    rel_dy(i) = dy_sd(i)/dy_sd(i);
    rel_dc(i) = dc_sd(i)/dy_sd(i);
    rel_dinve(i) = dinve_sd(i)/dy_sd(i);
    rel_dw(i) = dw_sd(i)/dy_sd(i);
    rel_labobs(i) = labobs_sd(i)/dy_sd(i);
    rel_pinfobs(i) = pinfobs_sd(i)/dy_sd(i);
    rel_robs(i) = robs_sd(i)/dy_sd(i);
    
    %Calculate correlation matrix
    correlation_matrix(:,:,i) = corrcoef(data);
    
    %Calculate cross-correlations
    for k = 1:11

        %Output growth
        ccor_dy_dc(k,i) = cross_correlation(dy,dc,k-6);
        ccor_dy_dinve(k,i) = cross_correlation(dy,dinve,k-6);
        ccor_dy_dw(k,i) = cross_correlation(dy,dw,k-6);
        ccor_dy_labobs(k,i) = cross_correlation(dy,labobs,k-6);
        ccor_dy_pinfobs(k,i) = cross_correlation(dy,pinfobs,k-6);
        ccor_dy_robs(k,i) = cross_correlation(dy,robs,k-6);

        %Consumption growth
        ccor_dc_dy(k,i) = cross_correlation(dc,dy,k-6);
        ccor_dc_dinve(k,i) = cross_correlation(dc,dinve,k-6);
        ccor_dc_dw(k,i) = cross_correlation(dc,dw,k-6);
        ccor_dc_labobs(k,i) = cross_correlation(dc,labobs,k-6);
        ccor_dc_pinfobs(k,i) = cross_correlation(dc,pinfobs,k-6);
        ccor_dc_robs(k,i) = cross_correlation(dc,robs,k-6);

        %Investment growth
        ccor_dinve_dy(k,i) = cross_correlation(dinve,dy,k-6);
        ccor_dinve_dc(k,i) = cross_correlation(dinve,dc,k-6);
        ccor_dinve_dw(k,i) = cross_correlation(dinve,dw,k-6);
        ccor_dinve_labobs(k,i) = cross_correlation(dinve,labobs,k-6);
        ccor_dinve_pinfobs(k,i) = cross_correlation(dinve,pinfobs,k-6);
        ccor_dinve_robs(k,i) = cross_correlation(dinve,robs,k-6);

        %Real wage growth
        ccor_dw_dy(k,i) = cross_correlation(dw,dy,k-6);
        ccor_dw_dc(k,i) = cross_correlation(dw,dc,k-6);
        ccor_dw_dinve(k,i) = cross_correlation(dw,dinve,k-6);
        ccor_dw_labobs(k,i) = cross_correlation(dw,labobs,k-6);
        ccor_dw_pinfobs(k,i) = cross_correlation(dw,pinfobs,k-6);
        ccor_dw_robs(k,i) = cross_correlation(dw,robs,k-6);

        %Hours worked
        ccor_labobs_dy(k,i) = cross_correlation(labobs,dy,k-6);
        ccor_labobs_dc(k,i) = cross_correlation(labobs,dc,k-6);
        ccor_labobs_dinve(k,i) = cross_correlation(labobs,dinve,k-6);
        ccor_labobs_dw(k,i) = cross_correlation(labobs,dw,k-6);
        ccor_labobs_pinfobs(k,i) = cross_correlation(labobs,pinfobs,k-6);
        ccor_labobs_robs(k,i) = cross_correlation(labobs,robs,k-6);

        %Inflation
        ccor_pinfobs_dy(k,i) = cross_correlation(pinfobs,dy,k-6);
        ccor_pinfobs_dc(k,i) = cross_correlation(pinfobs,dc,k-6);
        ccor_pinfobs_dinve(k,i) = cross_correlation(pinfobs,dinve,k-6);
        ccor_pinfobs_dw(k,i) = cross_correlation(pinfobs,dw,k-6);
        ccor_pinfobs_labobs(k,i) = cross_correlation(pinfobs,labobs,k-6);
        ccor_pinfobs_robs(k,i) = cross_correlation(pinfobs,robs,k-6);

        %Nominal interest rate
        ccor_robs_dy(k,i) = cross_correlation(robs,dy,k-6);
        ccor_robs_dc(k,i) = cross_correlation(robs,dc,k-6);
        ccor_robs_dinve(k,i) = cross_correlation(robs,dinve,k-6);
        ccor_robs_dw(k,i) = cross_correlation(robs,dw,k-6);
        ccor_robs_labobs(k,i) = cross_correlation(robs,labobs,k-6);
        ccor_robs_pinfobs(k,i) = cross_correlation(robs,pinfobs,k-6);

    end
   
    %Errors
    errors(i,:) = sum(error_count);

end

%% Mean and Standard Deviation

%Output growth
results_dy_mean = median(dy_mean);
results_dy_sd = median(dy_sd);
results_dy_acf = median(dy_acf,2);

%Consumption growth
results_dc_mean = median(dc_mean);
results_dc_sd = median(dc_sd);
results_dc_acf = median(dc_acf,2);

%Investment growth
results_dinve_mean = median(dinve_mean);
results_dinve_sd = median(dinve_sd);
results_dinve_acf = median(dinve_acf,2);

%Real wage growth
results_dw_mean = median(dw_mean);
results_dw_sd = median(dw_sd);
results_dw_acf = median(dw_acf,2);

%Hours worked
results_labobs_mean = median(labobs_mean);
results_labobs_sd = median(labobs_sd);
results_labobs_acf = median(labobs_acf,2);

%Inflation
results_pinfobs_mean = median(pinfobs_mean);
results_pinfobs_sd = median(pinfobs_sd);
results_pinfobs_acf = median(pinfobs_acf,2);

%Nominal interest rate
results_robs_mean = median(robs_mean);
results_robs_sd = median(robs_sd);
results_robs_acf = median(robs_acf,2);

%Relative standard deviations (to dy)
results_rel_dy = median(rel_dy);
results_rel_dc = median(rel_dc);
results_rel_dinve = median(rel_dinve);
results_rel_dw = median(rel_dw);
results_rel_labobs = median(rel_labobs);
results_rel_pinfobs = median(rel_pinfobs);
results_rel_robs = median(rel_robs);

%% Correlation matrix

results_correlation_matrix = median(correlation_matrix,3);

%% Cross-Correlations

%Output growth
results_ccor_dy_dc = median(ccor_dy_dc,2);
results_ccor_dy_dinve = median(ccor_dy_dinve,2);
results_ccor_dy_dw = median(ccor_dy_dw,2);
results_ccor_dy_labobs = median(ccor_dy_labobs,2);
results_ccor_dy_pinfobs = median(ccor_dy_pinfobs,2);
results_ccor_dy_robs = median(ccor_dy_robs,2);

%Consumption growth
results_ccor_dc_dy = median(ccor_dc_dy,2);
results_ccor_dc_dinve = median(ccor_dc_dinve,2);
results_ccor_dc_dw = median(ccor_dc_dw,2);
results_ccor_dc_labobs = median(ccor_dc_labobs,2);
results_ccor_dc_pinfobs = median(ccor_dc_pinfobs,2);
results_ccor_dc_robs = median(ccor_dc_robs,2);

%Investment growth
results_ccor_dinve_dy = median(ccor_dinve_dy,2);
results_ccor_dinve_dc = median(ccor_dinve_dc,2);
results_ccor_dinve_dw = median(ccor_dinve_dw,2);
results_ccor_dinve_labobs = median(ccor_dinve_labobs,2);
results_ccor_dinve_pinfobs = median(ccor_dinve_pinfobs,2);
results_ccor_dinve_robs = median(ccor_dinve_robs,2);

%Real wage growth
results_ccor_dw_dy = median(ccor_dw_dy,2);
results_ccor_dw_dc = median(ccor_dw_dc,2);
results_ccor_dw_dinve = median(ccor_dw_dinve,2);
results_ccor_dw_labobs = median(ccor_dw_labobs,2);
results_ccor_dw_pinfobs = median(ccor_dw_pinfobs,2);
results_ccor_dw_robs = median(ccor_dw_robs,2);

%Hours worked
results_ccor_labobs_dy = median(ccor_labobs_dy,2);
results_ccor_labobs_dc = median(ccor_labobs_dc,2);
results_ccor_labobs_dinve = median(ccor_labobs_dinve,2);
results_ccor_labobs_dw = median(ccor_labobs_dw,2);
results_ccor_labobs_pinfobs = median(ccor_labobs_pinfobs,2);
results_ccor_labobs_robs = median(ccor_labobs_robs,2);

%Inflation
results_ccor_pinfobs_dy = median(ccor_pinfobs_dy,2);
results_ccor_pinfobs_dc = median(ccor_pinfobs_dc,2);
results_ccor_pinfobs_dinve = median(ccor_pinfobs_dinve,2);
results_ccor_pinfobs_dw = median(ccor_pinfobs_dw,2);
results_ccor_pinfobs_labobs = median(ccor_pinfobs_labobs,2);
results_ccor_pinfobs_robs = median(ccor_pinfobs_robs,2);

%Nominal interest rate
results_ccor_robs_dy = median(ccor_robs_dy,2);
results_ccor_robs_dc = median(ccor_robs_dc,2);
results_ccor_robs_dinve = median(ccor_robs_dinve,2);
results_ccor_robs_dw = median(ccor_robs_dw,2);
results_ccor_robs_labobs = median(ccor_robs_labobs,2);
results_ccor_robs_pinfobs = median(ccor_robs_pinfobs,2);
