%%
%Medium-scale New Keynesian Model with Heterogeneous Expectations
%Generates IRFs

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
parameter_info = load('msnk_alh_parameter_estimates.csv');

%Number of periods to generate the IRFs
T = 10;

%% Create storage for the impulse responses and forecast error variance decompositions

%Impulse responses
resp = zeros(number_state_variables,number_exogenous_variables,T+2);

%Forecast error variance decompositions
%Gross price Markup
fevd_mc = zeros(T+1,number_exogenous_variables+1);
fevd_mc_prop = zeros(T+1,number_exogenous_variables+1);

%Capital utilization rate
fevd_zcap = zeros(T+1,number_exogenous_variables+1);
fevd_zcap_prop = zeros(T+1,number_exogenous_variables+1);

%Rental rate of capital
fevd_rk = zeros(T+1,number_exogenous_variables+1);
fevd_rk_prop = zeros(T+1,number_exogenous_variables+1);

%Capital services
fevd_k = zeros(T+1,number_exogenous_variables+1);
fevd_k_prop = zeros(T+1,number_exogenous_variables+1);

%Investment
fevd_inve = zeros(T+1,number_exogenous_variables+1);
fevd_inve_prop = zeros(T+1,number_exogenous_variables+1);

%Real value of existing capital stock
fevd_pk = zeros(T+1,number_exogenous_variables+1);
fevd_pk_prop = zeros(T+1,number_exogenous_variables+1);

%Consumption
fevd_c = zeros(T+1,number_exogenous_variables+1);
fevd_c_prop = zeros(T+1,number_exogenous_variables+1);

%Output
fevd_y = zeros(T+1,number_exogenous_variables+1);
fevd_y_prop = zeros(T+1,number_exogenous_variables+1);

%Hours worked 
fevd_lab = zeros(T+1,number_exogenous_variables+1);
fevd_lab_prop = zeros(T+1,number_exogenous_variables+1);

%Inflation 
fevd_pinf = zeros(T+1,number_exogenous_variables+1);
fevd_pinf_prop = zeros(T+1,number_exogenous_variables+1);

%Real wage
fevd_w = zeros(T+1,number_exogenous_variables+1);
fevd_w_prop = zeros(T+1,number_exogenous_variables+1);

%Nominal interest rate
fevd_r = zeros(T+1,number_exogenous_variables+1);
fevd_r_prop = zeros(T+1,number_exogenous_variables+1);

%Capital stock
fevd_kp = zeros(T+1,number_exogenous_variables+1);
fevd_kp_prop = zeros(T+1,number_exogenous_variables+1);

%% Parameters

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

%% Find the heterogeneous expectations equilibrium and build the state-space matrices that aren't a function of agent beliefs

[H_1,hee,K_1,K_2,K_3,K_4,Sigma_epsilon,t,Psi_0,Psi_1,Psi_2,S_c] = msnk_alh_build_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,theta);

%% Initialize all variables

% Initialize the adaptive learning algorithm
[s,A,B,R_A,R_B,~,~,error_count] = msnk_alh_initialize_learning(number_endogenous_variables,number_exogenous_variables,number_state_variables,T+2,hee);

% Initialize the state-space matrices
S_1 = zeros(number_state_variables,number_state_variables,T+2);
S_2 = zeros(number_state_variables,number_exogenous_variables,T+2);

%Generate the shocks (a one standard deviation shock in the first period (i.e., period 3), zero values for the shocks in all subsequent periods
epsilon = zeros(number_exogenous_variables,T+2);
epsilon(:,3) = ones(number_exogenous_variables,1); 

% Store the initial values of the S_1 and S_2 matrices by updating the adaptive learning state-space matrices that are functions of agent beliefs with the given values of the parameters and the initial values of the agent beliefs
[S_1(:,:,2),S_2(:,:,2),~,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,2),B(:,:,2),R_A(:,:,2),R_B(:,:,2),hee,A(:,:,1),B(:,:,1),R_A(:,:,1),R_B(:,:,1),S_1(:,:,1),S_2(:,:,1));

%% Generate the IRFs by simulating the model

for j = 3:T+2
    
    if j == 3
    
    %Get the last two periods of the state variables
    s_last_two_periods = [s(:,j-2),s(:,j-1)];
    
    %Now update the agent beliefs and moment matrices (the dependent variables are last period's state variables)
    [A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_beliefs(theta,s_last_two_periods,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1));
    
    %Now use the updated beliefs to update the adaptive learning state-space matrices that are functions of agent beliefs
    [S_1(:,:,j),S_2(:,:,j),test,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j),hee,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1),S_1(:,:,j-1),S_2(:,:,j-1));
    
    %Now update the values of the state variables to reflect the updated S_1 and S_2 matrices and the one standard deviation shock to all variables
    s(:,j) = S_1(:,:,j)*s(:,j-1) + S_2(:,:,j)*epsilon(:,j);
    
    %Store the response of the variables to the shocks
    resp(:,:,j) = S_2(:,:,j);
    
    %Record the errors for this iteration
    error_count(j,:) = test;
        
    else
        
    %Get the last two periods of the state variables
    s_last_two_periods = [s(:,j-2),s(:,j-1)];
    
    %Now update the agent beliefs and moment matrices (the dependent variables are last period's state variables)
    [A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_beliefs(theta,s_last_two_periods,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1));
    
    %Now use the updated beliefs to update the adaptive learning state-space matrices that are functions of agent beliefs
    [S_1(:,:,j),S_2(:,:,j),test,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j)] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A(:,:,j),B(:,:,j),R_A(:,:,j),R_B(:,:,j),hee,A(:,:,j-1),B(:,:,j-1),R_A(:,:,j-1),R_B(:,:,j-1),S_1(:,:,j-1),S_2(:,:,j-1));
    
    %Now update the values of the state variables to reflect the updated S_1 matrix and the fact that the shocks go to zero
    s(:,j) = S_1(:,:,j)*s(:,j-1) + S_2(:,:,j)*epsilon(:,j);
    
    %Store the response of the variables to the shocks
    resp(:,:,j) = S_1(:,:,j)*resp(:,:,j-1);
    
    %Record the errors for this iteration
    error_count(j,:) = test;
    
    end
end

%% IRFs

%Productivity shock
resp_a_mc(:,1) = (squeeze(resp(12,1,2:end))); %Gross price markup 
resp_a_zcap(:,1) = (squeeze(resp(13,1,2:end))); %Capital utilization rate 
resp_a_rk(:,1) = (squeeze(resp(14,1,2:end))); %Rental rate of capital 
resp_a_k(:,1) = (squeeze(resp(15,1,2:end))); %Capital services 
resp_a_inve(:,1) = (squeeze(resp(16,1,2:end))); %Investment 
resp_a_pk(:,1) = (squeeze(resp(17,1,2:end))); %Real value of existing capital stock 
resp_a_c(:,1) = (squeeze(resp(18,1,2:end))); %Consumption 
resp_a_y(:,1) = (squeeze(resp(19,1,2:end))); %Output 
resp_a_lab(:,1) = (squeeze(resp(20,1,2:end))); %Hours worked 
resp_a_pinf(:,1) = (squeeze(resp(21,1,2:end))); %Inflation 
resp_a_w(:,1) = (squeeze(resp(22,1,2:end))); %Real wage 
resp_a_r(:,1) = (squeeze(resp(23,1,2:end))); %Nominal interest rate 
resp_a_kp(:,1) = (squeeze(resp(24,1,2:end))); %Capital stock 

%Risk premium shock
resp_b_mc(:,1) = (squeeze(resp(12,2,2:end))); %Gross price markup 
resp_b_zcap(:,1) = (squeeze(resp(13,2,2:end))); %Capital utilization rate 
resp_b_rk(:,1) = (squeeze(resp(14,2,2:end))); %Rental rate of capital 
resp_b_k(:,1) = (squeeze(resp(15,2,2:end))); %Capital services 
resp_b_inve(:,1) = (squeeze(resp(16,2,2:end))); %Investment 
resp_b_pk(:,1) = (squeeze(resp(17,2,2:end))); %Real value of existing capital stock 
resp_b_c(:,1) = (squeeze(resp(18,2,2:end))); %Consumption 
resp_b_y(:,1) = (squeeze(resp(19,2,2:end))); %Output 
resp_b_lab(:,1) = (squeeze(resp(20,2,2:end))); %Hours worked 
resp_b_pinf(:,1) = (squeeze(resp(21,2,2:end))); %Inflation 
resp_b_w(:,1) = (squeeze(resp(22,2,2:end))); %Real wage 
resp_b_r(:,1) = (squeeze(resp(23,2,2:end))); %Nominal interest rate 
resp_b_kp(:,1) = (squeeze(resp(24,2,2:end))); %Capital stock 

%Exogenous spending shock
resp_g_mc(:,1) = (squeeze(resp(12,3,2:end))); %Gross price markup 
resp_g_zcap(:,1) = (squeeze(resp(13,3,2:end))); %Capital utilization rate 
resp_g_rk(:,1) = (squeeze(resp(14,3,2:end))); %Rental rate of capital 
resp_g_k(:,1) = (squeeze(resp(15,3,2:end))); %Capital services 
resp_g_inve(:,1) = (squeeze(resp(16,3,2:end))); %Investment 
resp_g_pk(:,1) = (squeeze(resp(17,3,2:end))); %Real value of existing capital stock 
resp_g_c(:,1) = (squeeze(resp(18,3,2:end))); %Consumption 
resp_g_y(:,1) = (squeeze(resp(19,3,2:end))); %Output 
resp_g_lab(:,1) = (squeeze(resp(20,3,2:end))); %Hours worked 
resp_g_pinf(:,1) = (squeeze(resp(21,3,2:end))); %Inflation 
resp_g_w(:,1) = (squeeze(resp(22,3,2:end))); %Real wage 
resp_g_r(:,1) = (squeeze(resp(23,3,2:end))); %Nominal interest rate 
resp_g_kp(:,1) = (squeeze(resp(24,3,2:end))); %Capital stock 

%Investment-specific technology shock
resp_qs_mc(:,1) = (squeeze(resp(12,4,2:end))); %Gross price markup 
resp_qs_zcap(:,1) = (squeeze(resp(13,4,2:end))); %Capital utilization rate 
resp_qs_rk(:,1) = (squeeze(resp(14,4,2:end))); %Rental rate of capital 
resp_qs_k(:,1) = (squeeze(resp(15,4,2:end))); %Capital services 
resp_qs_inve(:,1) = (squeeze(resp(16,4,2:end))); %Investment 
resp_qs_pk(:,1) = (squeeze(resp(17,4,2:end))); %Real value of existing capital stock 
resp_qs_c(:,1) = (squeeze(resp(18,4,2:end))); %Consumption 
resp_qs_y(:,1) = (squeeze(resp(19,4,2:end))); %Output 
resp_qs_lab(:,1) = (squeeze(resp(20,4,2:end))); %Hours worked 
resp_qs_pinf(:,1) = (squeeze(resp(21,4,2:end))); %Inflation 
resp_qs_w(:,1) = (squeeze(resp(22,4,2:end))); %Real wage 
resp_qs_r(:,1) = (squeeze(resp(23,4,2:end))); %Nominal interest rate 
resp_qs_kp(:,1) = (squeeze(resp(24,4,2:end))); %Capital stock 

%Monetary policy shock
resp_ms_mc(:,1) = (squeeze(resp(12,5,2:end))); %Gross price markup 
resp_ms_zcap(:,1) = (squeeze(resp(13,5,2:end))); %Capital utilization rate 
resp_ms_rk(:,1) = (squeeze(resp(14,5,2:end))); %Rental rate of capital 
resp_ms_k(:,1) = (squeeze(resp(15,5,2:end))); %Capital services 
resp_ms_inve(:,1) = (squeeze(resp(16,5,2:end))); %Investment 
resp_ms_pk(:,1) = (squeeze(resp(17,5,2:end))); %Real value of existing capital stock
resp_ms_c(:,1) = (squeeze(resp(18,5,2:end))); %Consumption 
resp_ms_y(:,1) = (squeeze(resp(19,5,2:end))); %Output 
resp_ms_lab(:,1) = (squeeze(resp(20,5,2:end))); %Hours worked 
resp_ms_pinf(:,1) = (squeeze(resp(21,5,2:end))); %Inflation 
resp_ms_w(:,1) = (squeeze(resp(22,5,2:end))); %Real wage 
resp_ms_r(:,1) = (squeeze(resp(23,5,2:end))); %Nominal interest rate 
resp_ms_kp(:,1) = (squeeze(resp(24,5,2:end))); %Capital stock

%Price markup shock
resp_spinf_mc(:,1) = (squeeze(resp(12,6,2:end))); %Gross price markup 
resp_spinf_zcap(:,1) = (squeeze(resp(13,6,2:end))); %Capital utilization rate 
resp_spinf_rk(:,1) = (squeeze(resp(14,6,2:end))); %Rental rate of capital 
resp_spinf_k(:,1) = (squeeze(resp(15,6,2:end))); %Capital services 
resp_spinf_inve(:,1) = (squeeze(resp(16,6,2:end))); %Investment 
resp_spinf_pk(:,1) = (squeeze(resp(17,6,2:end))); %Real value of existing capital stock 
resp_spinf_c(:,1) = (squeeze(resp(18,6,2:end))); %Consumption 
resp_spinf_y(:,1) = (squeeze(resp(19,6,2:end))); %Output 
resp_spinf_lab(:,1) = (squeeze(resp(20,6,2:end))); %Hours worked 
resp_spinf_pinf(:,1) = (squeeze(resp(21,6,2:end))); %Inflation 
resp_spinf_w(:,1) = (squeeze(resp(22,6,2:end))); %Real wage 
resp_spinf_r(:,1) = (squeeze(resp(23,6,2:end))); %Nominal interest rate 
resp_spinf_kp(:,1) = (squeeze(resp(24,6,2:end))); %Capital stock 

%Wage markup shock
resp_sw_mc(:,1) = (squeeze(resp(12,7,2:end))); %Gross price markup 
resp_sw_zcap(:,1) = (squeeze(resp(13,7,2:end))); %Capital utilization rate 
resp_sw_rk(:,1) = (squeeze(resp(14,7,2:end))); %Rental rate of capital 
resp_sw_k(:,1) = (squeeze(resp(15,7,2:end))); %Capital services 
resp_sw_inve(:,1) = (squeeze(resp(16,7,2:end))); %Investment 
resp_sw_pk(:,1) = (squeeze(resp(17,7,2:end))); %Real value of existing capital stock 
resp_sw_c(:,1) = (squeeze(resp(18,7,2:end))); %Consumption 
resp_sw_y(:,1) = (squeeze(resp(19,7,2:end))); %Output 
resp_sw_lab(:,1) = (squeeze(resp(20,7,2:end))); %Hours worked 
resp_sw_pinf(:,1) = (squeeze(resp(21,7,2:end))); %Inflation 
resp_sw_w(:,1) = (squeeze(resp(22,7,2:end))); %Real wage 
resp_sw_r(:,1) = (squeeze(resp(23,7,2:end))); %Nominal interest rate 
resp_sw_kp(:,1) = (squeeze(resp(24,7,2:end))); %Capital stock 

%% Generate figures

%Productivity Shock
figure(1)
subplot(3,3,1)
plot(0:T,resp_a_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_a_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_a_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_a_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_a_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_a_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_a_r(:,1));
title('r'); grid on

%Risk premium shock
figure(2)
subplot(3,3,1)
plot(0:T,resp_b_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_b_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_b_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_b_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_b_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_b_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_b_r(:,1));
title('r'); grid on

%Exogenous spending shock
figure(3)
subplot(3,3,1)
plot(0:T,resp_g_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_g_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_g_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_g_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_g_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_g_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_g_r(:,1));
title('r'); grid on

%Investment-specific technology shock
figure(4)
subplot(3,3,1)
plot(0:T,resp_qs_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_qs_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_qs_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_qs_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_qs_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_qs_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_qs_r(:,1));
title('r'); grid on

%Monetary policy shock
figure(5)
subplot(3,3,1)
plot(0:T,resp_ms_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_ms_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_ms_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_ms_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_ms_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_ms_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_ms_r(:,1));
title('r'); grid on

%Price markup shock
figure(6)
subplot(3,3,1)
plot(0:T,resp_spinf_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_spinf_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_spinf_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_spinf_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_spinf_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_spinf_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_spinf_r(:,1));
title('r'); grid on

%Wage markup shock
figure(7)
subplot(3,3,1)
plot(0:T,resp_sw_y(:,1))
title('y'); grid on
subplot(3,3,2)
plot(0:T,resp_sw_c(:,1))
title('c'); grid on
subplot(3,3,3)
plot(0:T,resp_sw_inve(:,1));
title('inve'); grid on
subplot(3,3,4)
plot(0:T,resp_sw_w(:,1));
title('w'); grid on
subplot(3,3,5)
plot(0:T,resp_sw_lab(:,1));
title('lab'); grid on
subplot(3,3,6)
plot(0:T,resp_sw_pinf(:,1));
title('pinf'); grid on
subplot(3,3,7)
plot(0:T,resp_sw_r(:,1));
title('r'); grid on

%% Forecast Error Variance Decomposition

%This comes from the powerpoint presentation "DSGE Simulation Techniques" by
%Martin Ellison; 

for i = 2:T+1
    
    %Gross price Markup
    fevd_mc(i,1) = fevd_mc(i-1,1) + (resp_a_mc(i)^2)*(sig_ea^2);
    fevd_mc(i,2) = fevd_mc(i-1,2) + (resp_b_mc(i)^2)*(sig_eb^2);
    fevd_mc(i,3) = fevd_mc(i-1,3) + (resp_g_mc(i)^2)*(sig_eg^2);
    fevd_mc(i,4) = fevd_mc(i-1,4) + (resp_qs_mc(i)^2)*(sig_eqs^2);
    fevd_mc(i,5) = fevd_mc(i-1,5) + (resp_ms_mc(i)^2)*(sig_em^2);
    fevd_mc(i,6) = fevd_mc(i-1,6) + (resp_spinf_mc(i)^2)*(sig_epinf^2);
    fevd_mc(i,7) = fevd_mc(i-1,7) + (resp_sw_mc(i)^2)*(sig_ew^2);
    fevd_mc(i,8) = fevd_mc(i,1) + fevd_mc(i,2) + fevd_mc(i,3) + fevd_mc(i,4) + fevd_mc(i,5) + fevd_mc(i,6) + fevd_mc(i,7);
    
    fevd_mc_prop(i,1) = fevd_mc(i,1)/fevd_mc(i,8);
    fevd_mc_prop(i,2) = fevd_mc(i,2)/fevd_mc(i,8);
    fevd_mc_prop(i,3) = fevd_mc(i,3)/fevd_mc(i,8);
    fevd_mc_prop(i,4) = fevd_mc(i,4)/fevd_mc(i,8);
    fevd_mc_prop(i,5) = fevd_mc(i,5)/fevd_mc(i,8);
    fevd_mc_prop(i,6) = fevd_mc(i,6)/fevd_mc(i,8);
    fevd_mc_prop(i,7) = fevd_mc(i,7)/fevd_mc(i,8);
    fevd_mc_prop(i,8) = sum(fevd_mc_prop(i,1:7));
    
    %Capital utilization rate
    fevd_zcap(i,1) = fevd_zcap(i-1,1) + (resp_a_zcap(i)^2)*(sig_ea^2);
    fevd_zcap(i,2) = fevd_zcap(i-1,2) + (resp_b_zcap(i)^2)*(sig_eb^2);
    fevd_zcap(i,3) = fevd_zcap(i-1,3) + (resp_g_zcap(i)^2)*(sig_eg^2);
    fevd_zcap(i,4) = fevd_zcap(i-1,4) + (resp_qs_zcap(i)^2)*(sig_eqs^2);
    fevd_zcap(i,5) = fevd_zcap(i-1,5) + (resp_ms_zcap(i)^2)*(sig_em^2);
    fevd_zcap(i,6) = fevd_zcap(i-1,6) + (resp_spinf_zcap(i)^2)*(sig_epinf^2);
    fevd_zcap(i,7) = fevd_zcap(i-1,7) + (resp_sw_zcap(i)^2)*(sig_ew^2);
    fevd_zcap(i,8) = fevd_zcap(i,1) + fevd_zcap(i,2) + fevd_zcap(i,3) + fevd_zcap(i,4) + fevd_zcap(i,5) + fevd_zcap(i,6) + fevd_zcap(i,7);
    
    fevd_zcap_prop(i,1) = fevd_zcap(i,1)/fevd_zcap(i,8);
    fevd_zcap_prop(i,2) = fevd_zcap(i,2)/fevd_zcap(i,8);
    fevd_zcap_prop(i,3) = fevd_zcap(i,3)/fevd_zcap(i,8);
    fevd_zcap_prop(i,4) = fevd_zcap(i,4)/fevd_zcap(i,8);
    fevd_zcap_prop(i,5) = fevd_zcap(i,5)/fevd_zcap(i,8);
    fevd_zcap_prop(i,6) = fevd_zcap(i,6)/fevd_zcap(i,8);
    fevd_zcap_prop(i,7) = fevd_zcap(i,7)/fevd_zcap(i,8);
    fevd_zcap_prop(i,8) = sum(fevd_zcap_prop(i,1:7));
    
    %Rental rate of capital  
    fevd_rk(i,1) = fevd_rk(i-1,1) + (resp_a_rk(i)^2)*(sig_ea^2);
    fevd_rk(i,2) = fevd_rk(i-1,2) + (resp_b_rk(i)^2)*(sig_eb^2);
    fevd_rk(i,3) = fevd_rk(i-1,3) + (resp_g_rk(i)^2)*(sig_eg^2);
    fevd_rk(i,4) = fevd_rk(i-1,4) + (resp_qs_rk(i)^2)*(sig_eqs^2);
    fevd_rk(i,5) = fevd_rk(i-1,5) + (resp_ms_rk(i)^2)*(sig_em^2);
    fevd_rk(i,6) = fevd_rk(i-1,6) + (resp_spinf_rk(i)^2)*(sig_epinf^2);
    fevd_rk(i,7) = fevd_rk(i-1,7) + (resp_sw_rk(i)^2)*(sig_ew^2);
    fevd_rk(i,8) = fevd_rk(i,1) + fevd_rk(i,2) + fevd_rk(i,3) + fevd_rk(i,4) + fevd_rk(i,5) + fevd_rk(i,6) + fevd_rk(i,7);
    
    fevd_rk_prop(i,1) = fevd_rk(i,1)/fevd_rk(i,8);
    fevd_rk_prop(i,2) = fevd_rk(i,2)/fevd_rk(i,8);
    fevd_rk_prop(i,3) = fevd_rk(i,3)/fevd_rk(i,8);
    fevd_rk_prop(i,4) = fevd_rk(i,4)/fevd_rk(i,8);
    fevd_rk_prop(i,5) = fevd_rk(i,5)/fevd_rk(i,8);
    fevd_rk_prop(i,6) = fevd_rk(i,6)/fevd_rk(i,8);
    fevd_rk_prop(i,7) = fevd_rk(i,7)/fevd_rk(i,8);
    fevd_rk_prop(i,8) = sum(fevd_rk_prop(i,1:7));
    
    %Capital services
    fevd_k(i,1) = fevd_k(i-1,1) + (resp_a_k(i)^2)*(sig_ea^2);
    fevd_k(i,2) = fevd_k(i-1,2) + (resp_b_k(i)^2)*(sig_eb^2);
    fevd_k(i,3) = fevd_k(i-1,3) + (resp_g_k(i)^2)*(sig_eg^2);
    fevd_k(i,4) = fevd_k(i-1,4) + (resp_qs_k(i)^2)*(sig_eqs^2);
    fevd_k(i,5) = fevd_k(i-1,5) + (resp_ms_k(i)^2)*(sig_em^2);
    fevd_k(i,6) = fevd_k(i-1,6) + (resp_spinf_k(i)^2)*(sig_epinf^2);
    fevd_k(i,7) = fevd_k(i-1,7) + (resp_sw_k(i)^2)*(sig_ew^2);
    fevd_k(i,8) = fevd_k(i,1) + fevd_k(i,2) + fevd_k(i,3) + fevd_k(i,4) + fevd_k(i,5) + fevd_k(i,6) + fevd_k(i,7);
    
    fevd_k_prop(i,1) = fevd_k(i,1)/fevd_k(i,8);
    fevd_k_prop(i,2) = fevd_k(i,2)/fevd_k(i,8);
    fevd_k_prop(i,3) = fevd_k(i,3)/fevd_k(i,8);
    fevd_k_prop(i,4) = fevd_k(i,4)/fevd_k(i,8);
    fevd_k_prop(i,5) = fevd_k(i,5)/fevd_k(i,8);
    fevd_k_prop(i,6) = fevd_k(i,6)/fevd_k(i,8);
    fevd_k_prop(i,7) = fevd_k(i,7)/fevd_k(i,8);
    fevd_k_prop(i,8) = sum(fevd_k_prop(i,1:7));
    
    %Investment
    fevd_inve(i,1) = fevd_inve(i-1,1) + (resp_a_inve(i)^2)*(sig_ea^2);
    fevd_inve(i,2) = fevd_inve(i-1,2) + (resp_b_inve(i)^2)*(sig_eb^2);
    fevd_inve(i,3) = fevd_inve(i-1,3) + (resp_g_inve(i)^2)*(sig_eg^2);
    fevd_inve(i,4) = fevd_inve(i-1,4) + (resp_qs_inve(i)^2)*(sig_eqs^2);
    fevd_inve(i,5) = fevd_inve(i-1,5) + (resp_ms_inve(i)^2)*(sig_em^2);
    fevd_inve(i,6) = fevd_inve(i-1,6) + (resp_spinf_inve(i)^2)*(sig_epinf^2);
    fevd_inve(i,7) = fevd_inve(i-1,7) + (resp_sw_inve(i)^2)*(sig_ew^2);
    fevd_inve(i,8) = fevd_inve(i,1) + fevd_inve(i,2) + fevd_inve(i,3) + fevd_inve(i,4) + fevd_inve(i,5) + fevd_inve(i,6) + fevd_inve(i,7);
    
    fevd_inve_prop(i,1) = fevd_inve(i,1)/fevd_inve(i,8);
    fevd_inve_prop(i,2) = fevd_inve(i,2)/fevd_inve(i,8);
    fevd_inve_prop(i,3) = fevd_inve(i,3)/fevd_inve(i,8);
    fevd_inve_prop(i,4) = fevd_inve(i,4)/fevd_inve(i,8);
    fevd_inve_prop(i,5) = fevd_inve(i,5)/fevd_inve(i,8);
    fevd_inve_prop(i,6) = fevd_inve(i,6)/fevd_inve(i,8);
    fevd_inve_prop(i,7) = fevd_inve(i,7)/fevd_inve(i,8);
    fevd_inve_prop(i,8) = sum(fevd_inve_prop(i,1:7));
    
    %Real value of existing capital stock
    fevd_pk(i,1) = fevd_pk(i-1,1) + (resp_a_pk(i)^2)*(sig_ea^2);
    fevd_pk(i,2) = fevd_pk(i-1,2) + (resp_b_pk(i)^2)*(sig_eb^2);
    fevd_pk(i,3) = fevd_pk(i-1,3) + (resp_g_pk(i)^2)*(sig_eg^2);
    fevd_pk(i,4) = fevd_pk(i-1,4) + (resp_qs_pk(i)^2)*(sig_eqs^2);
    fevd_pk(i,5) = fevd_pk(i-1,5) + (resp_ms_pk(i)^2)*(sig_em^2);
    fevd_pk(i,6) = fevd_pk(i-1,6) + (resp_spinf_pk(i)^2)*(sig_epinf^2);
    fevd_pk(i,7) = fevd_pk(i-1,7) + (resp_sw_pk(i)^2)*(sig_ew^2);
    fevd_pk(i,8) = fevd_pk(i,1) + fevd_pk(i,2) + fevd_pk(i,3) + fevd_pk(i,4) + fevd_pk(i,5) + fevd_pk(i,6) + fevd_pk(i,7);
    
    fevd_pk_prop(i,1) = fevd_pk(i,1)/fevd_pk(i,8);
    fevd_pk_prop(i,2) = fevd_pk(i,2)/fevd_pk(i,8);
    fevd_pk_prop(i,3) = fevd_pk(i,3)/fevd_pk(i,8);
    fevd_pk_prop(i,4) = fevd_pk(i,4)/fevd_pk(i,8);
    fevd_pk_prop(i,5) = fevd_pk(i,5)/fevd_pk(i,8);
    fevd_pk_prop(i,6) = fevd_pk(i,6)/fevd_pk(i,8);
    fevd_pk_prop(i,7) = fevd_pk(i,7)/fevd_pk(i,8);
    fevd_pk_prop(i,8) = sum(fevd_pk_prop(i,1:7));
    
    %Consumption
    fevd_c(i,1) = fevd_c(i-1,1) + (resp_a_c(i)^2)*(sig_ea^2);
    fevd_c(i,2) = fevd_c(i-1,2) + (resp_b_c(i)^2)*(sig_eb^2);
    fevd_c(i,3) = fevd_c(i-1,3) + (resp_g_c(i)^2)*(sig_eg^2);
    fevd_c(i,4) = fevd_c(i-1,4) + (resp_qs_c(i)^2)*(sig_eqs^2);
    fevd_c(i,5) = fevd_c(i-1,5) + (resp_ms_c(i)^2)*(sig_em^2);
    fevd_c(i,6) = fevd_c(i-1,6) + (resp_spinf_c(i)^2)*(sig_epinf^2);
    fevd_c(i,7) = fevd_c(i-1,7) + (resp_sw_c(i)^2)*(sig_ew^2);
    fevd_c(i,8) = fevd_c(i,1) + fevd_c(i,2) + fevd_c(i,3) + fevd_c(i,4) + fevd_c(i,5) + fevd_c(i,6) + fevd_c(i,7);
    
    fevd_c_prop(i,1) = fevd_c(i,1)/fevd_c(i,8);
    fevd_c_prop(i,2) = fevd_c(i,2)/fevd_c(i,8);
    fevd_c_prop(i,3) = fevd_c(i,3)/fevd_c(i,8);
    fevd_c_prop(i,4) = fevd_c(i,4)/fevd_c(i,8);
    fevd_c_prop(i,5) = fevd_c(i,5)/fevd_c(i,8);
    fevd_c_prop(i,6) = fevd_c(i,6)/fevd_c(i,8);
    fevd_c_prop(i,7) = fevd_c(i,7)/fevd_c(i,8);
    fevd_c_prop(i,8) = sum(fevd_c_prop(i,1:7));
    
    %Output 
    fevd_y(i,1) = fevd_y(i-1,1) + (resp_a_y(i)^2)*(sig_ea^2);
    fevd_y(i,2) = fevd_y(i-1,2) + (resp_b_y(i)^2)*(sig_eb^2);
    fevd_y(i,3) = fevd_y(i-1,3) + (resp_g_y(i)^2)*(sig_eg^2);
    fevd_y(i,4) = fevd_y(i-1,4) + (resp_qs_y(i)^2)*(sig_eqs^2);
    fevd_y(i,5) = fevd_y(i-1,5) + (resp_ms_y(i)^2)*(sig_em^2);
    fevd_y(i,6) = fevd_y(i-1,6) + (resp_spinf_y(i)^2)*(sig_epinf^2);
    fevd_y(i,7) = fevd_y(i-1,7) + (resp_sw_y(i)^2)*(sig_ew^2);
    fevd_y(i,8) = fevd_y(i,1) + fevd_y(i,2) + fevd_y(i,3) + fevd_y(i,4) + fevd_y(i,5) + fevd_y(i,6) + fevd_y(i,7);
    
    fevd_y_prop(i,1) = fevd_y(i,1)/fevd_y(i,8);
    fevd_y_prop(i,2) = fevd_y(i,2)/fevd_y(i,8);
    fevd_y_prop(i,3) = fevd_y(i,3)/fevd_y(i,8);
    fevd_y_prop(i,4) = fevd_y(i,4)/fevd_y(i,8);
    fevd_y_prop(i,5) = fevd_y(i,5)/fevd_y(i,8);
    fevd_y_prop(i,6) = fevd_y(i,6)/fevd_y(i,8);
    fevd_y_prop(i,7) = fevd_y(i,7)/fevd_y(i,8);
    fevd_y_prop(i,8) = sum(fevd_y_prop(i,1:7));
    
    %Hours worked
    fevd_lab(i,1) = fevd_lab(i-1,1) + (resp_a_lab(i)^2)*(sig_ea^2);
    fevd_lab(i,2) = fevd_lab(i-1,2) + (resp_b_lab(i)^2)*(sig_eb^2);
    fevd_lab(i,3) = fevd_lab(i-1,3) + (resp_g_lab(i)^2)*(sig_eg^2);
    fevd_lab(i,4) = fevd_lab(i-1,4) + (resp_qs_lab(i)^2)*(sig_eqs^2);
    fevd_lab(i,5) = fevd_lab(i-1,5) + (resp_ms_lab(i)^2)*(sig_em^2);
    fevd_lab(i,6) = fevd_lab(i-1,6) + (resp_spinf_lab(i)^2)*(sig_epinf^2);
    fevd_lab(i,7) = fevd_lab(i-1,7) + (resp_sw_lab(i)^2)*(sig_ew^2);
    fevd_lab(i,8) = fevd_lab(i,1) + fevd_lab(i,2) + fevd_lab(i,3) + fevd_lab(i,4) + fevd_lab(i,5) + fevd_lab(i,6) + fevd_lab(i,7);
    
    fevd_lab_prop(i,1) = fevd_lab(i,1)/fevd_lab(i,8);
    fevd_lab_prop(i,2) = fevd_lab(i,2)/fevd_lab(i,8);
    fevd_lab_prop(i,3) = fevd_lab(i,3)/fevd_lab(i,8);
    fevd_lab_prop(i,4) = fevd_lab(i,4)/fevd_lab(i,8);
    fevd_lab_prop(i,5) = fevd_lab(i,5)/fevd_lab(i,8);
    fevd_lab_prop(i,6) = fevd_lab(i,6)/fevd_lab(i,8);
    fevd_lab_prop(i,7) = fevd_lab(i,7)/fevd_lab(i,8);
    fevd_lab_prop(i,8) = sum(fevd_lab_prop(i,1:7));
    
    %Inflation
    fevd_pinf(i,1) = fevd_pinf(i-1,1) + (resp_a_pinf(i)^2)*(sig_ea^2);
    fevd_pinf(i,2) = fevd_pinf(i-1,2) + (resp_b_pinf(i)^2)*(sig_eb^2);
    fevd_pinf(i,3) = fevd_pinf(i-1,3) + (resp_g_pinf(i)^2)*(sig_eg^2);
    fevd_pinf(i,4) = fevd_pinf(i-1,4) + (resp_qs_pinf(i)^2)*(sig_eqs^2);
    fevd_pinf(i,5) = fevd_pinf(i-1,5) + (resp_ms_pinf(i)^2)*(sig_em^2);
    fevd_pinf(i,6) = fevd_pinf(i-1,6) + (resp_spinf_pinf(i)^2)*(sig_epinf^2);
    fevd_pinf(i,7) = fevd_pinf(i-1,7) + (resp_sw_pinf(i)^2)*(sig_ew^2);
    fevd_pinf(i,8) = fevd_pinf(i,1) + fevd_pinf(i,2) + fevd_pinf(i,3) + fevd_pinf(i,4) + fevd_pinf(i,5) + fevd_pinf(i,6) + fevd_pinf(i,7);
    
    fevd_pinf_prop(i,1) = fevd_pinf(i,1)/fevd_pinf(i,8);
    fevd_pinf_prop(i,2) = fevd_pinf(i,2)/fevd_pinf(i,8);
    fevd_pinf_prop(i,3) = fevd_pinf(i,3)/fevd_pinf(i,8);
    fevd_pinf_prop(i,4) = fevd_pinf(i,4)/fevd_pinf(i,8);
    fevd_pinf_prop(i,5) = fevd_pinf(i,5)/fevd_pinf(i,8);
    fevd_pinf_prop(i,6) = fevd_pinf(i,6)/fevd_pinf(i,8);
    fevd_pinf_prop(i,7) = fevd_pinf(i,7)/fevd_pinf(i,8);
    fevd_pinf_prop(i,8) = sum(fevd_pinf_prop(i,1:7));
    
    %Real wage
    fevd_w(i,1) = fevd_w(i-1,1) + (resp_a_w(i)^2)*(sig_ea^2);
    fevd_w(i,2) = fevd_w(i-1,2) + (resp_b_w(i)^2)*(sig_eb^2);
    fevd_w(i,3) = fevd_w(i-1,3) + (resp_g_w(i)^2)*(sig_eg^2);
    fevd_w(i,4) = fevd_w(i-1,4) + (resp_qs_w(i)^2)*(sig_eqs^2);
    fevd_w(i,5) = fevd_w(i-1,5) + (resp_ms_w(i)^2)*(sig_em^2);
    fevd_w(i,6) = fevd_w(i-1,6) + (resp_spinf_w(i)^2)*(sig_epinf^2);
    fevd_w(i,7) = fevd_w(i-1,7) + (resp_sw_w(i)^2)*(sig_ew^2);
    fevd_w(i,8) = fevd_w(i,1) + fevd_w(i,2) + fevd_w(i,3) + fevd_w(i,4) + fevd_w(i,5) + fevd_w(i,6) + fevd_w(i,7);
    
    fevd_w_prop(i,1) = fevd_w(i,1)/fevd_w(i,8);
    fevd_w_prop(i,2) = fevd_w(i,2)/fevd_w(i,8);
    fevd_w_prop(i,3) = fevd_w(i,3)/fevd_w(i,8);
    fevd_w_prop(i,4) = fevd_w(i,4)/fevd_w(i,8);
    fevd_w_prop(i,5) = fevd_w(i,5)/fevd_w(i,8);
    fevd_w_prop(i,6) = fevd_w(i,6)/fevd_w(i,8);
    fevd_w_prop(i,7) = fevd_w(i,7)/fevd_w(i,8);
    fevd_w_prop(i,8) = sum(fevd_w_prop(i,1:7));
    
    %Nominal interest rate
    fevd_r(i,1) = fevd_r(i-1,1) + (resp_a_r(i)^2)*(sig_ea^2);
    fevd_r(i,2) = fevd_r(i-1,2) + (resp_b_r(i)^2)*(sig_eb^2);
    fevd_r(i,3) = fevd_r(i-1,3) + (resp_g_r(i)^2)*(sig_eg^2);
    fevd_r(i,4) = fevd_r(i-1,4) + (resp_qs_r(i)^2)*(sig_eqs^2);
    fevd_r(i,5) = fevd_r(i-1,5) + (resp_ms_r(i)^2)*(sig_em^2);
    fevd_r(i,6) = fevd_r(i-1,6) + (resp_spinf_r(i)^2)*(sig_epinf^2);
    fevd_r(i,7) = fevd_r(i-1,7) + (resp_sw_r(i)^2)*(sig_ew^2);
    fevd_r(i,8) = fevd_r(i,1) + fevd_r(i,2) + fevd_r(i,3) + fevd_r(i,4) + fevd_r(i,5) + fevd_r(i,6) + fevd_r(i,7);
    
    fevd_r_prop(i,1) = fevd_r(i,1)/fevd_r(i,8);
    fevd_r_prop(i,2) = fevd_r(i,2)/fevd_r(i,8);
    fevd_r_prop(i,3) = fevd_r(i,3)/fevd_r(i,8);
    fevd_r_prop(i,4) = fevd_r(i,4)/fevd_r(i,8);
    fevd_r_prop(i,5) = fevd_r(i,5)/fevd_r(i,8);
    fevd_r_prop(i,6) = fevd_r(i,6)/fevd_r(i,8);
    fevd_r_prop(i,7) = fevd_r(i,7)/fevd_r(i,8);
    fevd_r_prop(i,8) = sum(fevd_r_prop(i,1:7));
    
    %Capital stock
    fevd_kp(i,1) = fevd_kp(i-1,1) + (resp_a_kp(i)^2)*(sig_ea^2);
    fevd_kp(i,2) = fevd_kp(i-1,2) + (resp_b_kp(i)^2)*(sig_eb^2);
    fevd_kp(i,3) = fevd_kp(i-1,3) + (resp_g_kp(i)^2)*(sig_eg^2);
    fevd_kp(i,4) = fevd_kp(i-1,4) + (resp_qs_kp(i)^2)*(sig_eqs^2);
    fevd_kp(i,5) = fevd_kp(i-1,5) + (resp_ms_kp(i)^2)*(sig_em^2);
    fevd_kp(i,6) = fevd_kp(i-1,6) + (resp_spinf_kp(i)^2)*(sig_epinf^2);
    fevd_kp(i,7) = fevd_kp(i-1,7) + (resp_sw_kp(i)^2)*(sig_ew^2);
    fevd_kp(i,8) = fevd_kp(i,1) + fevd_kp(i,2) + fevd_kp(i,3) + fevd_kp(i,4) + fevd_kp(i,5) + fevd_kp(i,6) + fevd_kp(i,7);
    
    fevd_kp_prop(i,1) = fevd_kp(i,1)/fevd_kp(i,8);
    fevd_kp_prop(i,2) = fevd_kp(i,2)/fevd_kp(i,8);
    fevd_kp_prop(i,3) = fevd_kp(i,3)/fevd_kp(i,8);
    fevd_kp_prop(i,4) = fevd_kp(i,4)/fevd_kp(i,8);
    fevd_kp_prop(i,5) = fevd_kp(i,5)/fevd_kp(i,8);
    fevd_kp_prop(i,6) = fevd_kp(i,6)/fevd_kp(i,8);
    fevd_kp_prop(i,7) = fevd_kp(i,7)/fevd_kp(i,8);
    fevd_kp_prop(i,8) = sum(fevd_kp_prop(i,1:7));

end