function [S_1_t,S_2_t,test,A_new,B_new,R_A_new,R_B_new] = msnk_alh_update_state_space_matrices(number_endogenous_variables,number_exogenous_variables,number_aux_variables,number_state_variables,K_1,K_2,K_3,K_4,H_1,A_t,B_t,R_A_t,R_B_t,hee,A_tm1,B_tm1,R_A_tm1,R_B_tm1,S_1_tm1,S_2_tm1)
%This function updates the state-space matrices with new values of agent beliefs for the medium-scale New Keynesian model with heterogeneous expectations 
%Input:
%number_endogenous_variables: number of endogenous variables
%number_exogenous_variables: number of exogenous variables
%number_aux_variables: number of auxiliary variables
%number_state_variables: number of state variables
%theta: vector of parameters
%K_1: matrix in condensed form equation
%K_2: matrix in condensed form equation
%K_3: matrix in condensed form equation
%K_4: matrix in condensed form equation
%H_1: matrix in condensed form equation
%A_t: agent-type A beliefs in current period
%B_t: agent-type B beliefs in current period
%R_A_t: agent-type A moment matrix in current period
%R_B_t: agent-type B moment matrix in current period
%hee: heterogeneous expectations equilibrium solution

%Output:
%S_1: matrix in state equation
%S_2: matrix in state equation
%test: flags used for NaN/inf,non-uniquensess, and non-stationarity
%A_new: agent-type A updated beliefs
%B_new: agent-type B updated beliefs
%R_A_new: agent-type A updated moment matrix
%R_B_new: agent-type B updated moment matrix

%% Build the belief matrices with the current period beliefs
%Zero matrices
Zero_n = zeros(number_endogenous_variables);
Zero_m = zeros(number_exogenous_variables);
Zero_k = zeros(number_aux_variables);
Zero_nm = zeros(number_endogenous_variables,number_exogenous_variables);
Zero_mn = zeros(number_exogenous_variables,number_endogenous_variables);
Zero_nk = zeros(number_endogenous_variables,number_aux_variables);
Zero_kn = zeros(number_aux_variables,number_endogenous_variables);
Zero_mk = zeros(number_exogenous_variables,number_aux_variables);
Zero_km = zeros(number_aux_variables,number_exogenous_variables);

%Build the A_1, A_2, and B matrices with the current period beliefs
A_1 = A_t(:,1:number_endogenous_variables);
A_2 = A_t(:,25:end);
B = B_t;

%Build the A_1_s, A_2_s, and B_s matrices
A_1_s = [A_1*A_1,Zero_nm,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
A_2_s = [Zero_n,A_1*A_2+A_2*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
B_s = [Zero_n,B*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];

%% State Equation Matrices

test_naninf = 0; %testing flag for NaN or infinity elements
test_noexist = 0; %testing flag for nonexistence of the solution
test_nonstationary = 0; %testing flag for nonstationarity of the solution

%Identity matrix
I_nmk = eye(number_state_variables);

%Calculate the S matrix with the current period beliefs
S = I_nmk - K_1*A_2_s - K_2*B_s;

%Test the S matrix to see if it has any NaN or infinity elements
[test_result_nan] = test_matrix_nan(S);
[test_result_inf] = test_matrix_inf(S);


if test_result_nan == 1 || test_result_inf == 1 %S has NaN and/or infinity elements; set the agent beliefs and moment matrices to the initial values and recalculate the S_1 and S_2 matrices
    
    %Set the "naninf" flag
    test_naninf = 1;

    %Set the agent beliefs and moment matrices to the initial values
    A_1 = hee.A_1_bar;
    A_2 = hee.A_2_bar;
    A_new = [A_1,A_2];
    B_new = hee.B_bar; 
    R_A_new = 0.1*eye(number_endogenous_variables+number_exogenous_variables);
    R_B_new = 0.1*eye(number_exogenous_variables);

    %Build the A_1_s, A_2_s, and B_s matrices
    A_1_s = [A_1*A_1,Zero_nm,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
    A_2_s = [Zero_n,A_1*A_2+A_2*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
    B_s = [Zero_n,B*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];

    %Re-calculate the S matrix
    S = I_nmk - K_1*A_2_s - K_2*B_s;

    %Calculate the S_1 and S_2 matrices
    S_1_t = S\(K_1*A_1_s + K_3);
    S_2_t = S\K_4;
  
else %S does not have any NaN or infinity elements
    
    %Determine if the solution exists (i.e., if S is invertible)
    [test_result_singular] = test_matrix_singular(S); 
    
    if test_result_singular == 1 %The solution does not exist (i.e., S is not invertible); set the agent beliefs and moment matrices to the initial values and recalculate the S_1 and S_2 matrices
        
        %Set the "no exist" flag
        test_noexist = 1;

        %Set the agent beliefs and moment matrices to the initial values
        A_1 = hee.A_1_bar;
        A_2 = hee.A_2_bar;
        A_new = [A_1,A_2];
        B_new = hee.B_bar; 
        R_A_new = 0.1*eye(number_endogenous_variables+number_exogenous_variables);
        R_B_new = 0.1*eye(number_exogenous_variables);

        %Build the A_1_s, A_2_s, and B_s matrices
        A_1_s = [A_1*A_1,Zero_nm,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
        A_2_s = [Zero_n,A_1*A_2+A_2*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];
        B_s = [Zero_n,B*H_1,Zero_nk;Zero_mn,Zero_m,Zero_mk;Zero_kn,Zero_km,Zero_k];

        %Re-calculate the S matrix
        S = I_nmk - K_1*A_2_s - K_2*B_s;

        %Calculate the S_1 and S_2 matrices
        S_1_t = S\(K_1*A_1_s + K_3);
        S_2_t = S\K_4;
        
    else %The solution exists (i.e., S is invertible); keep the S matrix and now calculate the S_1 matrix with current period beliefs
        
        S_1_t = S\(K_1*A_1_s + K_3);
        
        %Now test the S_1 matrix for stationarity       
        [test_result_nonstationarity] = test_matrix_nonstationarity(S_1_t);
        
        if test_result_nonstationarity == 1 %S_1 matrix is not stationary; set the agent beliefs, moment matrices, and S_1 and S_2 matrices to the previous period values
            
            %Set the "nonstationary" flag
            test_nonstationary = 1;

            %Set the agent beliefs, moment matrices, and S_1 and S_2 matrices to the previous period values
            A_new = A_tm1;
            B_new = B_tm1;
            R_A_new = R_A_tm1;
            R_B_new = R_B_tm1;
            S_1_t = S_1_tm1;
            S_2_t = S_2_tm1;
            
        else %The S_1 matrix is stationary; now calculate the S_2 matrix with the current period beliefs and keep the current period beliefs
            
            S_2_t = S\K_4;
            A_new = A_t;
            B_new = B_t;
            R_A_new = R_A_t;
            R_B_new = R_B_t;
            
        end
        
    end
    
end

test = [test_naninf,test_noexist,test_nonstationary];

end

