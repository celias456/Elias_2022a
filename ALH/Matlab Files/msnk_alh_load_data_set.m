function [data,Sigma_hat,c,first_observation,Sigma_u_sd,theta,prior_information] = msnk_alh_load_data_set(data_set_identifier)
%Loads the data set and associated characteristics for the Medium Scale New
%Keynesian model with heterogeneous expectations

%Input
%data_set_identifier: Data set used for estimation
    %1 = pre-great moderation
    %2 = great moderation
    %3 = pre-great recession
    %4 = full data set

%Parameters
% First entry is initial value
% Second entry is prior distribution number
% Third entry is prior distribution hyperparameter 1
% Fourth entry is prior distribution hyperparameter 2

%Output:
%data: variables in the data set
%Sigma_hat: variance of the jumping distribution used in the M-H algorithm
%c: scaling parameter used in the M-H algorithm
%first_observation: first observation used in the data set
%Sigma_v_sd: standard deviation of the measurement error term in measurement equation
%theta: vector of parameters
%prior_information: prior information for parameters in "theta" vector
%persistence: persistence parameter used in the projection facility in learning

if data_set_identifier == 1 %pre-great moderation
    
    load('msnk_data_pre_great_moderation.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
    data = [dy,dc,dinve,dw,labobs,pinfobs,robs]; 
    
    %Covariance matrix used for the variance of the proposal distribution in the MH algorithm (obtained from Dynare)
    Sigma_hat = load('msnk_alh_sigma_hat_pre_great_moderation.csv'); 
    Sigma_hat = nearestSPD(Sigma_hat); %The "nearestSPD" function ensures that the matrix P will be symmetric and postive semi-definite so that it can be used in the "mvnrnd" function
    
    %Load the parameter information
    parameter_info = load('msnk_alh_parameters_pre_great_moderation.csv');
    prior_value = parameter_info(:,1);
    prior_dist = parameter_info(:,2);
    hyp_1 = parameter_info(:,3);
    hyp_2 = parameter_info(:,4);
    
    sig_ea = [prior_value(1),prior_dist(1),hyp_1(1),hyp_2(1)];        
    sig_eb = [prior_value(2),prior_dist(2),hyp_1(2),hyp_2(2)];       
    sig_eg = [prior_value(3),prior_dist(3),hyp_1(3),hyp_2(3)]; 
    sig_eqs = [prior_value(4),prior_dist(4),hyp_1(4),hyp_2(4)];         
    sig_em = [prior_value(5),prior_dist(5),hyp_1(5),hyp_2(5)];              
    sig_epinf = [prior_value(6),prior_dist(6),hyp_1(6),hyp_2(6)];       
    sig_ew = [prior_value(7),prior_dist(7),hyp_1(7),hyp_2(7)];            
    csadjcost = [prior_value(8),prior_dist(8),hyp_1(8),hyp_2(8)];       
    csigma = [prior_value(9),prior_dist(9),hyp_1(9),hyp_2(9)];       
    chabb = [prior_value(10),prior_dist(10),hyp_1(10),hyp_2(10)];
    cprobw = [prior_value(11),prior_dist(11),hyp_1(11),hyp_2(11)];
    csigl = [prior_value(12),prior_dist(12),hyp_1(12),hyp_2(12)];
    cprobp = [prior_value(13),prior_dist(13),hyp_1(13),hyp_2(13)];
    cindw = [prior_value(14),prior_dist(14),hyp_1(14),hyp_2(14)];
    cindp = [prior_value(15),prior_dist(15),hyp_1(15),hyp_2(15)];
    czcap = [prior_value(16),prior_dist(16),hyp_1(16),hyp_2(16)];
    cfc = [prior_value(17),prior_dist(17),hyp_1(17),hyp_2(17)];
    crpi = [prior_value(18),prior_dist(18),hyp_1(18),hyp_2(18)];
    crr = [prior_value(19),prior_dist(19),hyp_1(19),hyp_2(19)];
    cry = [prior_value(20),prior_dist(20),hyp_1(20),hyp_2(20)];
    crdy = [prior_value(21),prior_dist(21),hyp_1(21),hyp_2(21)];
    constepinf = [prior_value(22),prior_dist(22),hyp_1(22),hyp_2(22)];
    constebeta = [prior_value(23),prior_dist(23),hyp_1(23),hyp_2(23)];
    constelab = [prior_value(24),prior_dist(24),hyp_1(24),hyp_2(24)];
    ctrend = [prior_value(25),prior_dist(25),hyp_1(25),hyp_2(25)];
    calfa = [prior_value(26),prior_dist(26),hyp_1(26),hyp_2(26)];
    crhoa = [prior_value(27),prior_dist(27),hyp_1(27),hyp_2(27)];
    crhob = [prior_value(28),prior_dist(28),hyp_1(28),hyp_2(28)];
    crhog = [prior_value(29),prior_dist(29),hyp_1(29),hyp_2(29)];
    crhoqs = [prior_value(30),prior_dist(30),hyp_1(30),hyp_2(30)];
    crhoms = [prior_value(31),prior_dist(31),hyp_1(31),hyp_2(31)];
    crhopinf = [prior_value(32),prior_dist(32),hyp_1(32),hyp_2(32)];
    crhow = [prior_value(33),prior_dist(33),hyp_1(33),hyp_2(33)];
    gna = [prior_value(34),prior_dist(34),hyp_1(34),hyp_2(34)];
    gnb = [prior_value(35),prior_dist(35),hyp_1(35),hyp_2(35)];
    omegaa = [prior_value(36),prior_dist(36),hyp_1(36),hyp_2(36)];
    
    %Scaling parameter; increase this to have a lower acceptance rate in the MH algorithm
    c = 0.02;
    
elseif data_set_identifier == 2 %great moderation
    load('msnk_data_great_moderation.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
    data = [dy,dc,dinve,dw,labobs,pinfobs,robs]; 
    
    %Covariance matrix used for the variance of the proposal distribution in the MH algorithm (obtained from Dynare)
    Sigma_hat = load('msnk_alh_sigma_hat_great_moderation.csv'); 
    Sigma_hat = nearestSPD(Sigma_hat); %The "nearestSPD" function ensures that the matrix P will be symmetric and postive semi-definite so that it can be used in the "mvnrnd" function
    
    %Load the parameter information
    parameter_info = load('msnk_alh_parameters_great_moderation.csv');
    prior_value = parameter_info(:,1);
    prior_dist = parameter_info(:,2);
    hyp_1 = parameter_info(:,3);
    hyp_2 = parameter_info(:,4);
    
    sig_ea = [prior_value(1),prior_dist(1),hyp_1(1),hyp_2(1)];        
    sig_eb = [prior_value(2),prior_dist(2),hyp_1(2),hyp_2(2)];       
    sig_eg = [prior_value(3),prior_dist(3),hyp_1(3),hyp_2(3)]; 
    sig_eqs = [prior_value(4),prior_dist(4),hyp_1(4),hyp_2(4)];         
    sig_em = [prior_value(5),prior_dist(5),hyp_1(5),hyp_2(5)];              
    sig_epinf = [prior_value(6),prior_dist(6),hyp_1(6),hyp_2(6)];       
    sig_ew = [prior_value(7),prior_dist(7),hyp_1(7),hyp_2(7)];            
    csadjcost = [prior_value(8),prior_dist(8),hyp_1(8),hyp_2(8)];       
    csigma = [prior_value(9),prior_dist(9),hyp_1(9),hyp_2(9)];       
    chabb = [prior_value(10),prior_dist(10),hyp_1(10),hyp_2(10)];
    cprobw = [prior_value(11),prior_dist(11),hyp_1(11),hyp_2(11)];
    csigl = [prior_value(12),prior_dist(12),hyp_1(12),hyp_2(12)];
    cprobp = [prior_value(13),prior_dist(13),hyp_1(13),hyp_2(13)];
    cindw = [prior_value(14),prior_dist(14),hyp_1(14),hyp_2(14)];
    cindp = [prior_value(15),prior_dist(15),hyp_1(15),hyp_2(15)];
    czcap = [prior_value(16),prior_dist(16),hyp_1(16),hyp_2(16)];
    cfc = [prior_value(17),prior_dist(17),hyp_1(17),hyp_2(17)];
    crpi = [prior_value(18),prior_dist(18),hyp_1(18),hyp_2(18)];
    crr = [prior_value(19),prior_dist(19),hyp_1(19),hyp_2(19)];
    cry = [prior_value(20),prior_dist(20),hyp_1(20),hyp_2(20)];
    crdy = [prior_value(21),prior_dist(21),hyp_1(21),hyp_2(21)];
    constepinf = [prior_value(22),prior_dist(22),hyp_1(22),hyp_2(22)];
    constebeta = [prior_value(23),prior_dist(23),hyp_1(23),hyp_2(23)];
    constelab = [prior_value(24),prior_dist(24),hyp_1(24),hyp_2(24)];
    ctrend = [prior_value(25),prior_dist(25),hyp_1(25),hyp_2(25)];
    calfa = [prior_value(26),prior_dist(26),hyp_1(26),hyp_2(26)];
    crhoa = [prior_value(27),prior_dist(27),hyp_1(27),hyp_2(27)];
    crhob = [prior_value(28),prior_dist(28),hyp_1(28),hyp_2(28)];
    crhog = [prior_value(29),prior_dist(29),hyp_1(29),hyp_2(29)];
    crhoqs = [prior_value(30),prior_dist(30),hyp_1(30),hyp_2(30)];
    crhoms = [prior_value(31),prior_dist(31),hyp_1(31),hyp_2(31)];
    crhopinf = [prior_value(32),prior_dist(32),hyp_1(32),hyp_2(32)];
    crhow = [prior_value(33),prior_dist(33),hyp_1(33),hyp_2(33)];
    gna = [prior_value(34),prior_dist(34),hyp_1(34),hyp_2(34)];
    gnb = [prior_value(35),prior_dist(35),hyp_1(35),hyp_2(35)];
    omegaa = [prior_value(36),prior_dist(36),hyp_1(36),hyp_2(36)];
    
    %Scaling parameter; increase this to have a lower acceptance rate in the MH algorithm
    c = 0.015;
    
elseif data_set_identifier == 3 %pre-great recession
    
    load('msnk_data_pre_great_recession.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
    data = [dy,dc,dinve,dw,labobs,pinfobs,robs];  
    
    %Covariance matrix used for the variance of the proposal distribution in the MH algorithm (obtained from Dynare)
    Sigma_hat = load('msnk_alh_sigma_hat_pre_great_recession.csv'); 
    Sigma_hat = nearestSPD(Sigma_hat); %The "nearestSPD" function ensures that the matrix P will be symmetric and postive semi-definite so that it can be used in the "mvnrnd" function
    
    %Load the parameter information
    parameter_info = load('msnk_alh_parameters_pre_great_recession.csv');
    prior_value = parameter_info(:,1);
    prior_dist = parameter_info(:,2);
    hyp_1 = parameter_info(:,3);
    hyp_2 = parameter_info(:,4);
    
    sig_ea = [prior_value(1),prior_dist(1),hyp_1(1),hyp_2(1)];        
    sig_eb = [prior_value(2),prior_dist(2),hyp_1(2),hyp_2(2)];       
    sig_eg = [prior_value(3),prior_dist(3),hyp_1(3),hyp_2(3)]; 
    sig_eqs = [prior_value(4),prior_dist(4),hyp_1(4),hyp_2(4)];         
    sig_em = [prior_value(5),prior_dist(5),hyp_1(5),hyp_2(5)];              
    sig_epinf = [prior_value(6),prior_dist(6),hyp_1(6),hyp_2(6)];       
    sig_ew = [prior_value(7),prior_dist(7),hyp_1(7),hyp_2(7)];            
    csadjcost = [prior_value(8),prior_dist(8),hyp_1(8),hyp_2(8)];       
    csigma = [prior_value(9),prior_dist(9),hyp_1(9),hyp_2(9)];       
    chabb = [prior_value(10),prior_dist(10),hyp_1(10),hyp_2(10)];
    cprobw = [prior_value(11),prior_dist(11),hyp_1(11),hyp_2(11)];
    csigl = [prior_value(12),prior_dist(12),hyp_1(12),hyp_2(12)];
    cprobp = [prior_value(13),prior_dist(13),hyp_1(13),hyp_2(13)];
    cindw = [prior_value(14),prior_dist(14),hyp_1(14),hyp_2(14)];
    cindp = [prior_value(15),prior_dist(15),hyp_1(15),hyp_2(15)];
    czcap = [prior_value(16),prior_dist(16),hyp_1(16),hyp_2(16)];
    cfc = [prior_value(17),prior_dist(17),hyp_1(17),hyp_2(17)];
    crpi = [prior_value(18),prior_dist(18),hyp_1(18),hyp_2(18)];
    crr = [prior_value(19),prior_dist(19),hyp_1(19),hyp_2(19)];
    cry = [prior_value(20),prior_dist(20),hyp_1(20),hyp_2(20)];
    crdy = [prior_value(21),prior_dist(21),hyp_1(21),hyp_2(21)];
    constepinf = [prior_value(22),prior_dist(22),hyp_1(22),hyp_2(22)];
    constebeta = [prior_value(23),prior_dist(23),hyp_1(23),hyp_2(23)];
    constelab = [prior_value(24),prior_dist(24),hyp_1(24),hyp_2(24)];
    ctrend = [prior_value(25),prior_dist(25),hyp_1(25),hyp_2(25)];
    calfa = [prior_value(26),prior_dist(26),hyp_1(26),hyp_2(26)];
    crhoa = [prior_value(27),prior_dist(27),hyp_1(27),hyp_2(27)];
    crhob = [prior_value(28),prior_dist(28),hyp_1(28),hyp_2(28)];
    crhog = [prior_value(29),prior_dist(29),hyp_1(29),hyp_2(29)];
    crhoqs = [prior_value(30),prior_dist(30),hyp_1(30),hyp_2(30)];
    crhoms = [prior_value(31),prior_dist(31),hyp_1(31),hyp_2(31)];
    crhopinf = [prior_value(32),prior_dist(32),hyp_1(32),hyp_2(32)];
    crhow = [prior_value(33),prior_dist(33),hyp_1(33),hyp_2(33)];
    gna = [prior_value(34),prior_dist(34),hyp_1(34),hyp_2(34)];
    gnb = [prior_value(35),prior_dist(35),hyp_1(35),hyp_2(35)];
    omegaa = [prior_value(36),prior_dist(36),hyp_1(36),hyp_2(36)];
    
    %Scaling parameter; increase this to have a lower acceptance rate in the MH algorithm
    c = 0.002;
    
else %full data set
    
    load('msnk_data_full_data_set.mat','dy','dc','dinve','dw','labobs','pinfobs','robs');
    data = [dy,dc,dinve,dw,labobs,pinfobs,robs];  
    
    %Covariance matrix used for the variance of the proposal distribution in the MH algorithm (obtained from Dynare)
    Sigma_hat = load('msnk_alh_sigma_hat_full_data_set.csv'); 
    Sigma_hat = nearestSPD(Sigma_hat); %The "nearestSPD" function ensures that the matrix P will be symmetric and postive semi-definite so that it can be used in the "mvnrnd" function
    
    %Load the parameter information
    parameter_info = load('msnk_alh_parameters_full_data_set.csv');
    prior_value = parameter_info(:,1);
    prior_dist = parameter_info(:,2);
    hyp_1 = parameter_info(:,3);
    hyp_2 = parameter_info(:,4);
    
    sig_ea = [prior_value(1),prior_dist(1),hyp_1(1),hyp_2(1)];        
    sig_eb = [prior_value(2),prior_dist(2),hyp_1(2),hyp_2(2)];       
    sig_eg = [prior_value(3),prior_dist(3),hyp_1(3),hyp_2(3)]; 
    sig_eqs = [prior_value(4),prior_dist(4),hyp_1(4),hyp_2(4)];         
    sig_em = [prior_value(5),prior_dist(5),hyp_1(5),hyp_2(5)];              
    sig_epinf = [prior_value(6),prior_dist(6),hyp_1(6),hyp_2(6)];       
    sig_ew = [prior_value(7),prior_dist(7),hyp_1(7),hyp_2(7)];            
    csadjcost = [prior_value(8),prior_dist(8),hyp_1(8),hyp_2(8)];       
    csigma = [prior_value(9),prior_dist(9),hyp_1(9),hyp_2(9)];       
    chabb = [prior_value(10),prior_dist(10),hyp_1(10),hyp_2(10)];
    cprobw = [prior_value(11),prior_dist(11),hyp_1(11),hyp_2(11)];
    csigl = [prior_value(12),prior_dist(12),hyp_1(12),hyp_2(12)];
    cprobp = [prior_value(13),prior_dist(13),hyp_1(13),hyp_2(13)];
    cindw = [prior_value(14),prior_dist(14),hyp_1(14),hyp_2(14)];
    cindp = [prior_value(15),prior_dist(15),hyp_1(15),hyp_2(15)];
    czcap = [prior_value(16),prior_dist(16),hyp_1(16),hyp_2(16)];
    cfc = [prior_value(17),prior_dist(17),hyp_1(17),hyp_2(17)];
    crpi = [prior_value(18),prior_dist(18),hyp_1(18),hyp_2(18)];
    crr = [prior_value(19),prior_dist(19),hyp_1(19),hyp_2(19)];
    cry = [prior_value(20),prior_dist(20),hyp_1(20),hyp_2(20)];
    crdy = [prior_value(21),prior_dist(21),hyp_1(21),hyp_2(21)];
    constepinf = [prior_value(22),prior_dist(22),hyp_1(22),hyp_2(22)];
    constebeta = [prior_value(23),prior_dist(23),hyp_1(23),hyp_2(23)];
    constelab = [prior_value(24),prior_dist(24),hyp_1(24),hyp_2(24)];
    ctrend = [prior_value(25),prior_dist(25),hyp_1(25),hyp_2(25)];
    calfa = [prior_value(26),prior_dist(26),hyp_1(26),hyp_2(26)];
    crhoa = [prior_value(27),prior_dist(27),hyp_1(27),hyp_2(27)];
    crhob = [prior_value(28),prior_dist(28),hyp_1(28),hyp_2(28)];
    crhog = [prior_value(29),prior_dist(29),hyp_1(29),hyp_2(29)];
    crhoqs = [prior_value(30),prior_dist(30),hyp_1(30),hyp_2(30)];
    crhoms = [prior_value(31),prior_dist(31),hyp_1(31),hyp_2(31)];
    crhopinf = [prior_value(32),prior_dist(32),hyp_1(32),hyp_2(32)];
    crhow = [prior_value(33),prior_dist(33),hyp_1(33),hyp_2(33)];
    gna = [prior_value(34),prior_dist(34),hyp_1(34),hyp_2(34)];
    gnb = [prior_value(35),prior_dist(35),hyp_1(35),hyp_2(35)];
    omegaa = [prior_value(36),prior_dist(36),hyp_1(36),hyp_2(36)];
    
    %Scaling parameter; increase this to have a lower acceptance rate in the MH algorithm
    c = 0.001;
    
end

%First observation in the data set
first_observation = 3;

%Standard deviation of measurement error
Sigma_u_sd = 0.0;

%Stack the parameters into a column vector
theta = [sig_ea(1);sig_eb(1);sig_eg(1);sig_eqs(1);sig_em(1);sig_epinf(1);sig_ew(1);csadjcost(1);csigma(1);chabb(1);cprobw(1);csigl(1);cprobp(1);cindw(1);cindp(1);czcap(1);cfc(1);crpi(1);crr(1);cry(1);crdy(1);constepinf(1);constebeta(1);constelab(1);ctrend(1);calfa(1);crhoa(1);crhob(1);crhog(1);crhoqs(1);crhoms(1);crhopinf(1);crhow(1);gna(1);gnb(1);omegaa(1)];

%Stack the prior information in a matrix
prior_information = [sig_ea(2:4);sig_eb(2:4);sig_eg(2:4);sig_eqs(2:4);sig_em(2:4);sig_epinf(2:4);sig_ew(2:4);csadjcost(2:4);csigma(2:4);chabb(2:4);cprobw(2:4);csigl(2:4);cprobp(2:4);cindw(2:4);cindp(2:4);czcap(2:4);cfc(2:4);crpi(2:4);crr(2:4);cry(2:4);crdy(2:4);constepinf(2:4);constebeta(2:4);constelab(2:4);ctrend(2:4);calfa(2:4);crhoa(2:4);crhob(2:4);crhog(2:4);crhoqs(2:4);crhoms(2:4);crhopinf(2:4);crhow(2:4);gna(2:4);gnb(2:4);omegaa(2:4)];
end