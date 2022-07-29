%%
%Medium-scale New Keynesian Model with Heterogeneous Expectations
%Forecasts with the model 

%This algorithm comes from the following 2 sources:
%1. "DSGE Model-Based Forecasting" by Del Negro and Schorfheide, pages 15-16
%2. "Forecasting with DSGE Models" by Christoffel, Coenen, and Warne, page 18

clear
clc

%% Characteristics of model

number_endogenous_variables = 24; %Number of endogenous variables
number_exogenous_variables = 7; %Number of exogenous variables
number_aux_variables = 4; %Number of auxiliary variables
number_observed_variables = 7; %Number of observable variables

%Number of variables in state transition equation
number_state_variables = number_endogenous_variables + number_exogenous_variables + number_aux_variables;

%Load the MH parameter draws
parameter_info = load('msnk_al_parameter_draws.csv');

%Get the number of observations in the estimation sample
T = 245;

%Load the data set
%load('msnk_data_pre_great_recession.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
load('msnk_data_full_data_set.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
data_estimation_sample = [dy(1:T),dc(1:T),dinve(1:T),dw(1:T),labobs(1:T),pinfobs(1:T),robs(1:T)];
data_forecast_sample = [dy(T+1:end),dc(T+1:end),dinve(T+1:end),dw(T+1:end),labobs(T+1:end),pinfobs(T+1:end),robs(T+1:end)];

%Number of periods-ahead to forecast
H = 43;

%Total number of observations in full data set
TpH = T + H;

%Number of simulations to run for each value of theta
m1 = 10000;

%Number of "thetas" to use
m2 = 1;

%Standard deviation of measurement error
Sigma_u_sd = 0.0;

%Set seed state
seed = 123;
rng(seed);

%% Create storage for summary statistics

%Storage for forecasts
forecasts = zeros(H,number_observed_variables,m1,m2);

%Storage for simulation errors
errors = zeros(m1,3,m2);

%% Start the "sampling the future" algorithm

for count_m2 = 1:m2
    
    %1. Draw theta
    
    %Fixed Parameters
    ctou = 0.025;
    cg = 0.18;
    clandaw = 1.5;
    curvp = 10;
    curvw = 10;
    
    % Estimated parameters
    sig_ea = parameter_info(1,count_m2);
    sig_eb = parameter_info(2,count_m2);
    sig_eg = parameter_info(3,count_m2);
    sig_eqs = parameter_info(4,count_m2);
    sig_em = parameter_info(5,count_m2);
    sig_epinf = parameter_info(6,count_m2);
    sig_ew = parameter_info(7,count_m2);
    csadjcost = parameter_info(8,count_m2);
    csigma = parameter_info(9,count_m2);
    chabb = parameter_info(10,count_m2);
    cprobw = parameter_info(11,count_m2);
    csigl = parameter_info(12,count_m2);
    cprobp = parameter_info(13,count_m2);
    cindw = parameter_info(14,count_m2);
    cindp = parameter_info(15,count_m2);
    czcap = parameter_info(16,count_m2);
    cfc = parameter_info(17,count_m2);
    crpi = parameter_info(18,count_m2);
    crr = parameter_info(19,count_m2);
    cry = parameter_info(20,count_m2);
    crdy = parameter_info(21,count_m2);
    constepinf = parameter_info(22,count_m2);
    constebeta = parameter_info(23,count_m2);
    constelab = parameter_info(24,count_m2);
    ctrend = parameter_info(25,count_m2);
    calfa = parameter_info(26,count_m2);
    crhoa = parameter_info(27,count_m2);
    crhob = parameter_info(28,count_m2);
    crhog = parameter_info(29,count_m2);
    crhoqs = parameter_info(30,count_m2);
    crhoms = parameter_info(31,count_m2);
    crhopinf = parameter_info(32,count_m2);
    crhow = parameter_info(33,count_m2);
    gna = parameter_info(34,count_m2);
    gnb = parameter_info(35,count_m2);
    omegaa = parameter_info(36,count_m2);

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

    
    %2. Get a draw of the state variables at time T using the Kalman filter
    
    %Find the heterogeneous expectations equilibrium and build the state-space matrices that aren't a function of agent beliefs
    [H_1,hee,K_1,K_2,K_3,K_4,Sigma_epsilon,t,Psi_0,Psi_1,Psi_2,S_c] = msnk_alh_build_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,theta);  

    %Initialize the adaptive learning algorithm
    [s,A,B,R_A,R_B,~,~,error_count] = msnk_alh_initialize_learning(number_endogenous_variables,number_exogenous_variables,number_state_variables,TpH,hee);

    % Initialize the State-space matrices
    S_1 = zeros(number_state_variables,number_state_variables,TpH);
    S_2 = zeros(number_state_variables,number_exogenous_variables,TpH);

    %Store the initial values of the S_1 and S_2 matrices by updating the adaptive learning state-space matrices that are functions of agent beliefs with the given values of the parameters and the initial values of the agent beliefs
    [S_1(:,:,2),S_2(:,:,2),~,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2),hee,A(:,:,1),B(:,:,1),R_A(:,:,1),R_B(:,:,1),S_1(:,:,1),S_2(:,:,1));

    %Initialize the Kalman filter
    [s(:,2),P_initial,~] = kalman_filter_initialize(number_state_variables,S_1(:,:,2),S_2(:,:,2),Sigma_epsilon,T);

    %Use the Kalman filter to get the values of the mean and covariance matrix of the state variables
    for period = 3:T
        
        %Get the new filter estimate of the mean and covariance matrix of the state variables for this particular data point
        [s(:,period),P,lik,~] = kalman_filter(data_estimation_sample(period,:)',s(:,period-1),P_initial,S_1(:,:,period-1),S_c,S_2(:,:,period-1),Psi_0,Psi_1,Psi_2,t,Sigma_epsilon,Sigma_u_sd);

        %Get the last two periods of the state variables
        s_last_two_periods = [s(:,period-2),s(:,period-1)];

        %Now update the agent beliefs and the moment matrices (the dependent variables are last period's state variables)
        [A(:,:,period),B(:,:,period),R_A(:,:,period),R_B(:,:,period)] = msnk_alh_update_beliefs(theta,s_last_two_periods,A(:,:,period-1),B(:,:,period-1),R_A(:,:,period-1),R_B(:,:,period-1));

        %Now use the updated beliefs to update the adaptive learning state-space matrices that are functions of agent beliefs
        [S_1(:,:,period),S_2(:,:,period),test,A(:,:,period),B(:,:,period),R_A(:,:,period),R_B(:,:,period)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,period),B(:,:,period),R_A(:,:,period),R_B(:,:,period),hee,A(:,:,period-1),B(:,:,period-1),R_A(:,:,period-1),R_B(:,:,period-1),S_1(:,:,period-1),S_2(:,:,period-1));

        %Update the MSE matrix for the next iteration
        P_initial = P;
        
        %Record the errors for this iteration
        error_count(period,:) = test;

    end
     
    %Get the filter estimate and the covariance matrix of the state variables
    s_filter_estimate = s(:,T);
    P_filter_estimate = nearestSPD(P);%The "nearestSPD" function ensures that the matrix P will be symmetric and postive semi-definite so that it can be used in the "mvnrnd" function
        
    for count_m1 = 1:m1

        %Generate a draw of the state variables at time T
        s(:,T-1) = mvnrnd(s_filter_estimate,P_filter_estimate)';
        s(:,T) = mvnrnd(s_filter_estimate,P_filter_estimate)';
        
        %3. Simulate a path of the state variables by using the draw of state variables as the initial values and generating a sequence of shocks

        %Initialize storage for individual state variables of interest
        y = zeros(1,TpH);
        c = zeros(1,TpH);
        inve = zeros(1,TpH);
        w = zeros(1,TpH);
        lab = zeros(1,TpH);
        pinf = zeros(1,TpH);
        r = zeros(1,TpH);

        %Initialize storage for individual observational equivalent variables
        model_dy = zeros(1,TpH);
        model_dc = zeros(1,TpH);
        model_dinve = zeros(1,TpH);
        model_dw = zeros(1,TpH);
        model_labobs = zeros(1,TpH);
        model_pinfobs = zeros(1,TpH);
        model_robs = zeros(1,TpH);

        %Generate the shocks
        epsilon = normrnd(0,1,number_exogenous_variables,TpH);

        %Simulate the model
        for j = (T+1):(TpH)   

            %Get the last two periods of the state variables
            s_last_two_periods = [s(:,j-2),s(:,j-1)];

            %Now update the agent beliefs (the dependent variables are last period's state variables)
            [A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_beliefs(theta,s_last_two_periods,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1));

            %Now use the updated beliefs to update the adaptive learning state-space matrices that are functions of agent beliefs
            [S_1(:,:,j),S_2(:,:,j),test,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j),hee,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1),S_1(:,:,j-1),S_2(:,:,j-1));

            %Now update the values of the state variables to reflect the updated S_1 and S_2 matrices
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
            model_dy(j) = y(j) - y(j-1) + ctrend;
            model_dc(j) = c(j) - c(j-1) + ctrend;
            model_dinve(j) = inve(j) - inve(j-1) + ctrend;
            model_dw(j) = w(j) - w(j-1) + ctrend;
            model_labobs(j) = lab(j) + constelab;
            model_pinfobs(j) = pinf(j) + constepinf;
            model_robs(j) = r(j) + conster;
            
            %Record the errors for this iteration
            error_count(j,:) = test;

        end

        %Drop the extraneous values
        model_dy = (model_dy(T+1:end))'; %Output growth
        model_dc = (model_dc(T+1:end))'; %Consumption growth
        model_dinve = (model_dinve(T+1:end))'; %Investment growth
        model_dw = (model_dw(T+1:end))'; %Real wage growth
        model_labobs = (model_labobs(T+1:end))'; %Labor hours
        model_pinfobs = (model_pinfobs(T+1:end))'; %Inflation
        model_robs = (model_robs(T+1:end))'; %Nominal interest rate

        forecasts(:,:,count_m1,count_m2) = [model_dy,model_dc,model_dinve,model_dw,model_labobs,model_pinfobs,model_robs];
        
        %Errors
        errors(count_m1,:,count_m2) = sum(error_count(T+1:TpH,:));

    end

end

forecasts_dy = reshape(forecasts(:,1,:,:),H,m1*m2);
forecasts_dc = reshape(forecasts(:,2,:,:),H,m1*m2);
forecasts_dinve = reshape(forecasts(:,3,:,:),H,m1*m2);
forecasts_dw = reshape(forecasts(:,4,:,:),H,m1*m2);
forecasts_labobs = reshape(forecasts(:,5,:,:),H,m1*m2);
forecasts_pinfobs = reshape(forecasts(:,6,:,:),H,m1*m2);
forecasts_robs = reshape(forecasts(:,7,:,:),H,m1*m2);

forecasts_point_dy = median(forecasts_dy,2);
forecasts_point_dc = median(forecasts_dc,2);
forecasts_point_dinve = median(forecasts_dinve,2);
forecasts_point_dw = median(forecasts_dw,2);
forecasts_point_labobs = median(forecasts_labobs,2);
forecasts_point_pinfobs = median(forecasts_pinfobs,2);
forecasts_point_robs = median(forecasts_robs,2);

forecasts_errors_squared_dy = (data_forecast_sample(:,1) - forecasts_point_dy).^2;
forecasts_errors_squared_dc = (data_forecast_sample(:,2) - forecasts_point_dc).^2;
forecasts_errors_squared_dinve = (data_forecast_sample(:,3) - forecasts_point_dinve).^2;
forecasts_errors_squared_dw = (data_forecast_sample(:,4) - forecasts_point_dw).^2;
forecasts_errors_squared_labobs = (data_forecast_sample(:,5) - forecasts_point_labobs).^2;
forecasts_errors_squared_pinfobs = (data_forecast_sample(:,6) - forecasts_point_pinfobs).^2;
forecasts_errors_squared_robs = (data_forecast_sample(:,7) - forecasts_point_robs).^2;

rmse.dy = sqrt(mean(forecasts_errors_squared_dy));
rmse.dc = sqrt(mean(forecasts_errors_squared_dc));
rmse.dinve = sqrt(mean(forecasts_errors_squared_dinve));
rmse.dw = sqrt(mean(forecasts_errors_squared_dw));
rmse.labobs = sqrt(mean(forecasts_errors_squared_labobs));
rmse.pinfobs = sqrt(mean(forecasts_errors_squared_pinfobs));
rmse.robs = sqrt(mean(forecasts_errors_squared_robs));