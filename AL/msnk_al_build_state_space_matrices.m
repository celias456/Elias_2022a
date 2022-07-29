function [H_1,hee,K_1,K_2,K_3,K_4,Sigma_epsilon,t,Psi_0,Psi_1,Psi_2,S_c] = msnk_al_build_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_observed_variables,number_state_variables,theta)
%This function determines the heterogeneous expectations equilibrium and builds
%the state-space matrices that aren't a function of agent beliefs of the
%medium scale New Keynesian model with heterogeneous expectations

%Input:
%number_endogenous_variables: number of endogenous variables
%number_exogenous_variables: number of exogenous variables
%number_aux_variables: number of auxiliary variables
%number_observed_variables: number of observed variables
%number_state_variables: number of state variables
%theta: vector of parameters

%Output:
%H_1: matrix in condensed form equation (identical to F_1) 
%hee: structure with four elements:
    %hee.A_1_bar: HEE expressions for the coefficients on agent-type A's endogenous variable regressors
    %hee.A_2_bar: HEE expressions for the coefficients on agent-type A's exogenous variable regressors
    %hee.B_bar: HEE expressions for the coefficients on agent-type B's exogenous variable regressor
    %hee.solution: equals 1 if solution is unique and stable, 0 otherwise
%K_1: matrix in condensed form equation
%K_2: matrix in condensed form equation
%K_3: matrix in condensed form equation
%K_4: matrix in condensed form equation
%Sigma_epsilon: sigma_epsilon matrix
%t: vector representing trend values in the measurement equation
%Psi_0: matrix in measurement equation
%Psi_1: matrix in measurement equation
%Psi_2: matrix in measurement equation
%S_c: state equation constants

%% Initialize all values (necessary in case the HEE solution does not exist)
H_1 = 0;
hee.solution = 0;
K_1 = 0;
K_2 = 0;
K_3 = 0;
K_4 = 0;
Sigma_epsilon = 0;
t = 0;
Psi_0 = 0;
Psi_1 = 0;
Psi_2 = 0;
S_c = 0;

%% Parameters

%Fixed Parameters
ctou = 0.025;
cg = 0.18;
clandaw = 1.5;
curvp = 10;
curvw = 10;
omegaa = 1;

%Parameters
sig_ea = theta(1);
sig_eb = theta(2);
sig_eg = theta(3);
sig_eqs = theta(4);
sig_em = theta(5);
sig_epinf = theta(6);
sig_ew = theta(7);
csadjcost = theta(8);
csigma = theta(9);
chabb = theta(10);
cprobw = theta(11);
csigl = theta(12);
cprobp = theta(13);
cindw = theta(14);
cindp = theta(15);
czcap = theta(16);
cfc = theta(17);
crpi = theta(18);
crr = theta(19);
cry = theta(20);
crdy = theta(21);
constepinf = theta(22);
constebeta = theta(23);
constelab = theta(24);
ctrend = theta(25);
calfa = theta(26);
crhoa = theta(27);
crhob = theta(28);
crhog = theta(29);
crhoqs = theta(30);
crhoms = theta(31);
crhopinf = theta(32);
crhow = theta(33);

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

%% Build Matrices

%Endogenous Variable Equation Numbers
eq_x_1 = 1;     %FOC labor with MPL expressed as a function of r^k and w (flexible price economy)
eq_x_2 = 2;     %FOC capacity utilization (flexible price economy)
eq_x_3 = 3;     %Firm FOC capital (flexible price economy)
eq_x_4 = 4;     %Definitiion capital services (flexible price economy)
eq_x_5 = 5;     %Investment Euler equation (flexible price economy)
eq_x_6 = 6;     %Arbitrage equation value of capital (flexible price economy)
eq_x_7 = 7;     %Consumption Euler equation (flexible price economy)
eq_x_8 = 8;     %Aggregate resource constraint (flexible price economy)
eq_x_9 = 9;     %Aggregate production function (flexible price economy)
eq_x_10 = 10;   %Wage equation (flexible price economy)
eq_x_11 = 11;   %Law of motion for capital (flexible price economy)
eq_x_12 = 12;   %FOC labor with MPL expressed as a function r^k and w 
eq_x_13 = 13;   %FOC capacity utilization
eq_x_14 = 14;   %Firm FOC capital
eq_x_15 = 15;   %Definition capital services
eq_x_16 = 16;   %Investment Euler equation
eq_x_17 = 17;   %Arbitrage equation value of capital
eq_x_18 = 18;   %Consumption Euler equation
eq_x_19 = 19;   %Aggregate resource constraint
eq_x_20 = 20;   %Aggregate production function
eq_x_21 = 21;   %New Keynesian Phillips curve
eq_x_22 = 22;   %Wage Phillips curve
eq_x_23 = 23;   %Taylor rule
eq_x_24 = 24;   %Law of motion of capital

% Endogenous variable numbers
rrf = 1;            %Real interest rate flexible price enonomy
zcapf = 2;          %Capital utilization rate flexible price economy
rkf = 3;            %Rental rate of capital flexible price economy
kf = 4;             %Capital services flexible price economy
invef = 5;          %Investment flexible price economy
pkf = 6;            %Real value of existing capital stock flexible price economy
cf = 7;             %Consumption flexible price economy
yf = 8;             %Output flexible price economy
labf = 9;           %Hours worked flexible price economy
wf = 10;            %Real wage flexible price economy
kpf = 11;           %Capital stock flexible price economy
mc = 12;            %Gross price markup
zcap = 13;          %Capital utilization rate
rk = 14;            %Rental rate of capital
k = 15;             %Capital services
inve = 16;          %Investment
pk = 17;            %Real value of existing capital stock
c = 18;             %Consumption
y = 19;             %Output
lab = 20;           %Hours worked
pinf = 21;          %Inflation
w = 22;             %Real wage
r = 23;             %Nominal interest rate
kp = 24;            %Capital stock

% Exogenous variable numbers
a = 1;     %Productivity shock
b = 2;     %Risk premium shock
g = 3;     %Exogenous spending shock
qs = 4;    %Investment specific technology shock
ms = 5;     %Monetary policy shock
spinf = 6;  %Price markup shock
sw = 7;     %Wage markup shock

% Initialize Matrices
D_0 = zeros(number_endogenous_variables,number_endogenous_variables);
D_1 = zeros(number_endogenous_variables,number_endogenous_variables);
D_2 = zeros(number_endogenous_variables,number_endogenous_variables);
D_3 = zeros(number_endogenous_variables,number_endogenous_variables);
D_4 = zeros(number_endogenous_variables,number_exogenous_variables);
F_0 = zeros(number_exogenous_variables,number_exogenous_variables);
F_1 = zeros(number_exogenous_variables,number_exogenous_variables);
F_2 = zeros(number_exogenous_variables,number_exogenous_variables);

%Endogenous Variable Equations
%Equation 1
D_0(eq_x_1,rkf) = -calfa;
D_0(eq_x_1,wf) = -(1-calfa);

D_4(eq_x_1,a) = -1;

%Equation 2
D_0(eq_x_2,zcapf) = 1;
D_0(eq_x_2,rkf) = -(1/(czcap/(1-czcap)));

%Equation 3
D_0(eq_x_3,rkf) = 1;
D_0(eq_x_3,kf) = 1;
D_0(eq_x_3,labf) = -1;
D_0(eq_x_3,wf) = -1;

%Equation 4
D_0(eq_x_4,zcapf) = -1;
D_0(eq_x_4,kf) = 1;

D_3(eq_x_4,kpf) = 1;

%Equation 5
D_0(eq_x_5,invef) = 1;
D_0(eq_x_5,pkf) = -(1/(1+cbetabar*cgamma))*(1/(cgamma^2*csadjcost));

D_1(eq_x_5,invef) = (1/(1+cbetabar*cgamma))*cbetabar*cgamma*omegaa;

D_2(eq_x_5,invef) = (1/(1+cbetabar*cgamma))*cbetabar*cgamma*omegab;

D_3(eq_x_5,invef) = (1/(1+cbetabar*cgamma));

D_4(eq_x_5,qs) = 1;

%Equation 6
D_0(eq_x_6,rrf) = 1;
D_0(eq_x_6,pkf) = 1;

D_1(eq_x_6,rkf) = (crk/(crk+(1-ctou)))*omegaa;
D_1(eq_x_6,pkf) = ((1-ctou)/(crk+(1-ctou)))*omegaa;

D_2(eq_x_6,rkf) = (crk/(crk+(1-ctou)))*omegab;
D_2(eq_x_6,pkf) = ((1-ctou)/(crk+(1-ctou)))*omegab;

D_4(eq_x_6,b) = (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))));

%Equation 7
D_0(eq_x_7,rrf) = (1-chabb/cgamma)/(csigma*(1+chabb/cgamma));
D_0(eq_x_7,cf) = 1;
D_0(eq_x_7,labf) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)));

D_1(eq_x_7,cf) = (1/(1+chabb/cgamma))*omegaa;
D_1(eq_x_7,labf) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*omegaa;

D_2(eq_x_7,cf) = (1/(1+chabb/cgamma))*omegab;
D_2(eq_x_7,labf) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*omegab;

D_3(eq_x_7,cf) = (chabb/cgamma)/(1+chabb/cgamma);

D_4(eq_x_7,b) = 1;

%Equation 8
D_0(eq_x_8,zcapf) = -crkky;
D_0(eq_x_8,invef) = -ciy;
D_0(eq_x_8,cf) = -ccy;
D_0(eq_x_8,yf) = 1;

D_4(eq_x_8,g) = 1;

%Equation 9
D_0(eq_x_9,kf) = -cfc*calfa;
D_0(eq_x_9,yf) = 1;
D_0(eq_x_9,labf) = -cfc*(1-calfa);

D_4(eq_x_9,a) = cfc;

%Equation 10
D_0(eq_x_10,cf) = -(1/(1-chabb/cgamma));
D_0(eq_x_10,labf) = -csigl;
D_0(eq_x_10,wf) = 1;

D_3(eq_x_10,cf) = -(chabb/cgamma)/(1-chabb/cgamma);

%Equation 11
D_0(eq_x_11,invef) = -cikbar;
D_0(eq_x_11,kpf) = 1;

D_3(eq_x_11,kpf) = (1-cikbar);

D_4(eq_x_11,qs) = (cikbar)*(1+cbetabar*cgamma)*(cgamma^2*csadjcost);

%Equation 12
D_0(eq_x_12,mc) = 1;
D_0(eq_x_12,rk) = -calfa;
D_0(eq_x_12,w) = -(1-calfa);

D_4(eq_x_12,a) = -1;

%Equation 13
D_0(eq_x_13,zcap) = 1;
D_0(eq_x_13,rk) = -(1/(czcap/(1-czcap)));

%Equation 14
D_0(eq_x_14,rk) = 1;
D_0(eq_x_14,k) = 1;
D_0(eq_x_14,lab) = -1;
D_0(eq_x_14,w) = -1;

%Equation 15
D_0(eq_x_15,zcap) = -1;
D_0(eq_x_15,k) = 1;

D_3(eq_x_15,kp) = 1;

%Equation 16
D_0(eq_x_16,inve) = 1;
D_0(eq_x_16,pk) = -(1/(1+cbetabar*cgamma))*(1/(cgamma^2*csadjcost));

D_1(eq_x_16,inve) = (1/(1+cbetabar*cgamma))*(cbetabar*cgamma)*omegaa;

D_2(eq_x_16,inve) = (1/(1+cbetabar*cgamma))*(cbetabar*cgamma)*omegab;

D_3(eq_x_16,inve) = (1/(1+cbetabar*cgamma));

D_4(eq_x_16,qs) = 1;

%Equation 17
D_0(eq_x_17,pk) = 1;
D_0(eq_x_17,r) = 1;

D_1(eq_x_17,rk) = (crk/(crk+(1-ctou)))*omegaa;
D_1(eq_x_17,pk) = ((1-ctou)/(crk+(1-ctou)))*omegaa;
D_1(eq_x_17,pinf) = omegaa;

D_2(eq_x_17,rk) = (crk/(crk+(1-ctou)))*omegab;
D_2(eq_x_17,pk) = ((1-ctou)/(crk+(1-ctou)))*omegab;
D_2(eq_x_17,pinf) = omegab;

D_4(eq_x_17,b) = (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))));

%Equation 18
D_0(eq_x_18,c) = 1;
D_0(eq_x_18,lab) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)));
D_0(eq_x_18,r) = (1-chabb/cgamma)/(csigma*(1+chabb/cgamma));

D_1(eq_x_18,c) = (1/(1+chabb/cgamma))*omegaa;
D_1(eq_x_18,lab) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*omegaa;
D_1(eq_x_18,pinf) = (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*omegaa;

D_2(eq_x_18,c) = (1/(1+chabb/cgamma))*omegab;
D_2(eq_x_18,lab) = -((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*omegab;
D_2(eq_x_18,pinf) = (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*omegab;

D_3(eq_x_18,c) = (chabb/cgamma)/(1+chabb/cgamma);

D_4(eq_x_18,b) = 1;

%Equation 19
D_0(eq_x_19,zcap) = -crkky;
D_0(eq_x_19,inve) = -ciy;
D_0(eq_x_19,c) = -ccy;
D_0(eq_x_19,y) = 1;

D_4(eq_x_19,g) = 1;

%Equation 20
D_0(eq_x_20,k) = -cfc*calfa;
D_0(eq_x_20,y) = 1;
D_0(eq_x_20,lab) = -cfc*(1-calfa);

D_4(eq_x_20,a) = cfc;

%Equation 21
D_0(eq_x_21,mc) = -(1/(1+cbetabar*cgamma*cindp))*((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1);
D_0(eq_x_21,pinf) = 1;

D_1(eq_x_21,pinf) = (1/(1+cbetabar*cgamma*cindp))*(cbetabar*cgamma)*omegaa;

D_2(eq_x_21,pinf) = (1/(1+cbetabar*cgamma*cindp))*(cbetabar*cgamma)*omegab;

D_3(eq_x_21,pinf) = (1/(1+cbetabar*cgamma*cindp))*cindp;

D_4(eq_x_21,spinf) = 1;

%Equation 22
D_0(eq_x_22,c) = -(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*(1/(1-chabb/cgamma));
D_0(eq_x_22,lab) = -(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*csigl;
D_0(eq_x_22,pinf) = (1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma);
D_0(eq_x_22,w) = 1+(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1));

D_1(eq_x_22,pinf) = (cbetabar*cgamma)/(1+cbetabar*cgamma)*omegaa;
D_1(eq_x_22,w) = (cbetabar*cgamma/(1+cbetabar*cgamma))*omegaa;

D_2(eq_x_22,pinf) = (cbetabar*cgamma)/(1+cbetabar*cgamma)*omegab;
D_2(eq_x_22,w) = (cbetabar*cgamma/(1+cbetabar*cgamma))*omegab;

D_3(eq_x_22,c) = -(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*((chabb/cgamma)/(1-chabb/cgamma));
D_3(eq_x_22,pinf) = (cindw/(1+cbetabar*cgamma));
D_3(eq_x_22,w) = (1/(1+cbetabar*cgamma));

D_4(eq_x_22,sw) = 1;

%Equation 23
D_0(eq_x_23,yf) = cry*(1-crr)+crdy;
D_0(eq_x_23,y) = -(cry*(1-crr)+crdy);
D_0(eq_x_23,pinf) = -crpi*(1-crr);
D_0(eq_x_23,r) = 1;

D_3(eq_x_23,yf) = crdy;
D_3(eq_x_23,y) = -crdy;
D_3(eq_x_23,r) = crr;

D_4(eq_x_23,ms) = 1;

%Equation 24
D_0(eq_x_24,inve) = -cikbar;
D_0(eq_x_24,kp) = 1;

D_3(eq_x_24,kp) = (1-cikbar);

D_4(eq_x_24,qs) = cikbar*(1+cbetabar*cgamma)*(cgamma^2)*csadjcost;

%Exogenous Variable Equations

%Equation 1 
F_0(a,a) = 1;
F_1(a,a) = crhoa;
F_2(a,a) = sig_ea;

%Equation 2
F_0(b,b) = 1;
F_1(b,b) = crhob;
F_2(b,b) = sig_eb;

%Equation 3
F_0(g,g) = 1;
F_1(g,g) = crhog;
F_2(g,g) = sig_eg;

%Equation 4
F_0(qs,qs) = 1;
F_1(qs,qs) = crhoqs;
F_2(qs,qs) = sig_eqs;

%Equation 5
F_0(ms,ms) = 1;
F_1(ms,ms) = crhoms;
F_2(ms,ms) = sig_em;

%Equation 6
F_0(spinf,spinf) = 1;
F_1(spinf,spinf) = crhopinf;
F_2(spinf,spinf) = sig_epinf;

%Equation 7
F_0(sw,sw) = 1;
F_1(sw,sw) = crhow;
F_2(sw,sw) = sig_ew;

%Determine if the D_0 and F_0 matrices are invertible
[test_result_singular_D_0] = test_matrix_singular(D_0);
[test_result_singular_F_0] = test_matrix_singular(F_0);

if test_result_singular_D_0 == 0 && test_result_singular_F_0 == 0 %D_0 and F_0 are invertible; Now calculate the G and H matrices and solve for the HEE

    G_1 = D_0\D_1;
    G_2 = D_0\D_2;
    G_3 = D_0\D_3;
    G_4 = D_0\D_4;
    H_1 = F_0\F_1;

    %% Solve for HE equilibrium
    %Identity Matrices
    I_n = eye(number_endogenous_variables);
    I_m = eye(number_exogenous_variables);
    I_nm = eye(number_endogenous_variables*number_exogenous_variables);

    % Solve for A_bar using the matrix quadratic solution function
    [A_1_bar,solution] = matrix_quadratic_solve(G_1,I_n,-G_3);

    if solution == 1 %A solution to the matrix quadratic equation was found
        %Solve for A_2 and B
        %Calculate the Sigma_v matrix
        Sigma_v = zeros(number_exogenous_variables,number_exogenous_variables);
        Sigma_v(1,1) = (sig_ea^2)/(1-crhoa^2);
        Sigma_v(2,2) = (sig_eb^2)/(1-crhob^2);
        Sigma_v(3,3) = (sig_eg^2)/(1-crhog^2);
        Sigma_v(4,4) = (sig_eqs^2)/(1-crhoqs^2);
        Sigma_v(5,5) = (sig_em^2)/(1-crhoms^2);
        Sigma_v(6,6) = (sig_epinf^2)/(1-crhopinf^2);
        Sigma_v(7,7) = (sig_ew^2)/(1-crhow^2);

        [Sigma_v_inv,~] = matrix_inverse(Sigma_v);

        %Calculate the P matrix
        P_1_hee = kron(G_3+G_1*A_1_bar*A_1_bar,Sigma_v_inv);
        P_2_hee = (I_nm - kron(G_3+G_1*A_1_bar*A_1_bar,H_1));
        [P_2_inv_hee,~] = matrix_inverse(P_2_hee);
        P_3_hee = kron(I_n,H_1*Sigma_v);
        P = P_1_hee*P_2_inv_hee*P_3_hee;

        %Calculate the Gamma matrix
        Gamma_1_hee = kron(G_1*A_1_bar,I_m) + kron(G_1,H_1') - I_nm;
        Gamma_2_hee = kron(G_2,H_1');
        Gamma_3_hee = (I_nm + P)*(kron(G_1*A_1_bar,I_m) + kron(G_1,H_1')); 
        Gamma_4_hee = (I_nm + P)*(kron(G_2,H_1')) - I_nm;
        Gamma_hee = [Gamma_1_hee,Gamma_2_hee;Gamma_3_hee,Gamma_4_hee];

        %Test for uniqueness of solution (i.e., invertibility of the Gamma matrix)
        test_unique = test_matrix_singular(Gamma_hee);
        if test_unique == 0 %Solution is unique; now get A_2_bar and B_bar
            Psi_1_hee = reshape(G_4',number_endogenous_variables*number_exogenous_variables,1);
            Psi_2_hee = (I_nm + P)*reshape(G_4',number_endogenous_variables*number_exogenous_variables,1);
            Psi_hee = [Psi_1_hee;Psi_2_hee];
            solution_1 = -Gamma_hee\Psi_hee;
            solution_2 = solution_1(1:number_endogenous_variables*number_exogenous_variables);
            A_2_bar = reshape(solution_2,number_exogenous_variables,number_endogenous_variables)';
            solution_3 = solution_1((number_endogenous_variables*number_exogenous_variables)+1:end);
            B_bar = reshape(solution_3,number_exogenous_variables,number_endogenous_variables)';

            %Now test for E-stability
            jac = kron(A_1_bar',G_1) + kron(I_n,G_1*A_1_bar-I_n);
            jac_estab = test_estability(jac);
            Gamma_hee_estab = test_estability(Gamma_hee);
            if jac_estab == 1 && Gamma_hee_estab == 1 %Solution is e-stable; now record the HEE solution
                %Get HEE solution
                hee.A_1_bar = A_1_bar;
                hee.A_2_bar = A_2_bar;
                hee.B_bar = B_bar;
                hee.solution = 1;

                % Now build the adaptive learning state-space matrices

                %Zero matrices
                Zero_m = zeros(number_exogenous_variables);
                Zero_k = zeros(number_aux_variables);
                Zero_nm = zeros(number_endogenous_variables,number_exogenous_variables);
                Zero_mn = zeros(number_exogenous_variables,number_endogenous_variables);
                Zero_nk = zeros(number_endogenous_variables,number_aux_variables);
                Zero_kn = zeros(number_aux_variables,number_endogenous_variables);
                Zero_mk = zeros(number_exogenous_variables,number_aux_variables);
                Zero_km = zeros(number_aux_variables,number_exogenous_variables);

                %Aux_0 and Aux_3 matrices
                Aux_0 = zeros(number_aux_variables,number_aux_variables);
                Aux_3 = zeros(number_aux_variables,number_endogenous_variables);

                % Auxiliary variable numbers
                ylaga = 1;
                claga = 2;
                invelaga = 3;
                wlaga = 4;

                %Auxiliary Variable Equations
                %Equation 1
                Aux_0(ylaga,ylaga) = 1;
                Aux_3(ylaga,y) = 1;

                %Equation 2
                Aux_0(claga,claga) = 1;
                Aux_3(claga,c) = 1;

                %Equation 3
                Aux_0(invelaga,invelaga) = 1;
                Aux_3(invelaga,inve) = 1;

                %Equation 4
                Aux_0(wlaga,wlaga) = 1;
                Aux_3(wlaga,w) = 1;

                %Calculate the adaptive learning state-space matrices
                J_0 = [D_0,-D_4,Zero_nk;Zero_mn,F_0,Zero_mk;Zero_kn,Zero_km,Aux_0];
                J_1 = [D_1,Zero_nm,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
                J_2 = [D_2,Zero_nm,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
                J_3 = [D_3,Zero_nm,Zero_nk;Zero_mn,F_1,Zero_mk;Aux_3,Zero_km,Zero_k];
                J_4 = [Zero_nm;F_2;Zero_km];

                K_1 = J_0\J_1;
                K_2 = J_0\J_2;
                K_3 = J_0\J_3;
                K_4 = J_0\J_4;

                %Now build the Sigma_epsilon matrix (i.e., the covariance matrix for the
                %stochastic shocks)
                Sigma_epsilon = eye(number_exogenous_variables);

                %Now build the measurement equation matrices

                %Vector representing trend values in the measurement equation
                t = ones(number_observed_variables,1);

                %Measurement equation numbers
                dy = 1;  %Output growth
                dc = 2;  %Consumption growth
                dinve = 3;  %Investment growth
                dw = 4;  %Real wage growth
                labobs = 5;  %Hours worked
                pinfobs = 6;  %Inflation
                robs = 7;  %Federal funds rate

                %Measurement equation variable numbers
                ylaga_s = 32;
                claga_s = 33;
                invelaga_s = 34;
                wlaga_s = 35;

                %Psi matrices
                Psi_0 = zeros(number_observed_variables,1);
                Psi_1 = zeros(number_observed_variables,number_observed_variables);
                Psi_2 = zeros(number_observed_variables,number_endogenous_variables+number_exogenous_variables+number_aux_variables);

                %Equation 1
                Psi_0(dy,1) = ctrend;

                Psi_2(dy,y) = 1;
                Psi_2(dy,ylaga_s) = -1;

                %Equation 2
                Psi_0(dc,1) = ctrend;

                Psi_2(dc,c) = 1;
                Psi_2(dc,claga_s) = -1;

                %Equation 3
                Psi_0(dinve,1) = ctrend;

                Psi_2(dinve,inve) = 1;
                Psi_2(dinve,invelaga_s) = -1;

                %Equation 4
                Psi_0(dw,1) = ctrend;

                Psi_2(dw,w) = 1;
                Psi_2(dw,wlaga_s) = -1;

                %Equation 5
                Psi_0(labobs,1) = constelab;

                Psi_2(labobs,lab) = 1;

                %Equation 6
                Psi_0(pinfobs,1) = constepinf;

                Psi_2(pinfobs,pinf) = 1;

                %Equation 7
                Psi_0(robs,1) = conster;

                Psi_2(robs,r) = 1;

                %State equation constants
                S_c = zeros(number_state_variables,1);
            end
        end
    end
end

end

