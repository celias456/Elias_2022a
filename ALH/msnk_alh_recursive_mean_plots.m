function [] = msnk_alh_recursive_mean_plots(theta)
%Plots the recursive means of each parameter for the medium scale new
%keynesian model with adaptive learning

rmean = recursive_mean(theta);

n = size(theta,1);
x = (1:n);

figure(7)

subplot(4,4,1), plot(x,rmean(:,1));%Plot of recursive mean of sigma_a
title('Recursive mean of $\sigma_a$','interpreter','latex')
subplot(4,4,2), plot(x,rmean(:,2));%Plot of recursive mean of sigma_b
title('Recursive mean of $\sigma_b$','interpreter','latex')
subplot(4,4,3), plot(x,rmean(:,3));%Plot of recursive mean of sigma_g
title('Recursive mean of $\sigma_g$','interpreter','latex')
subplot(4,4,4), plot(x,rmean(:,4));%Plot of recursive mean of sigma_i
title('Recursive mean of $\sigma_i$','interpreter','latex')
subplot(4,4,5), plot(x,rmean(:,5));%Plot of recursive mean of sigma_r
title('Recursive mean of $\sigma_r$','interpreter','latex')
subplot(4,4,6), plot(x,rmean(:,6));%Plot of recursive mean of sigma_p
title('Recursive mean of $\sigma_p$','interpreter','latex')
subplot(4,4,7), plot(x,rmean(:,7));%Plot of recursive mean of sigma_w
title('Recursive mean of $\sigma_w$','interpreter','latex')
subplot(4,4,8), plot(x,rmean(:,8));%Plot of recursive mean of varphi
title('Recursive mean of $\varphi$','interpreter','latex')
subplot(4,4,9), plot(x,rmean(:,9));%Plot of recursive mean of sigma_c
title('Recursive mean of $\sigma_c$','interpreter','latex')
subplot(4,4,10), plot(x,rmean(:,10));%Plot of recursive mean of lambda
title('Recursive mean of $\lambda$','interpreter','latex')
subplot(4,4,11), plot(x,rmean(:,11));%Plot of recursive mean of xi_w
title('Recursive mean of $\xi_w$','interpreter','latex')
subplot(4,4,12), plot(x,rmean(:,12));%Plot of recursive mean of sigma_l
title('Recursive mean of $\sigma_l$','interpreter','latex')
subplot(4,4,13), plot(x,rmean(:,13));%Plot of recursive mean of xi_p
title('Recursive mean of $\xi_p$','interpreter','latex')
subplot(4,4,14), plot(x,rmean(:,14));%Plot of recursive mean of iota_w
title('Recursive mean of $\iota_w$','interpreter','latex')
subplot(4,4,15), plot(x,rmean(:,15));%Plot of recursive mean of iota_p
title('Recursive mean of $\iota_p$','interpreter','latex')
subplot(4,4,16), plot(x,rmean(:,16));%Plot of recursive mean of psi
title('Recursive mean of $\psi$','interpreter','latex')



figure(8)

subplot(4,4,1), plot(x,rmean(:,17));%Plot of recursive mean of phi_p
title('Recursive mean of $\phi_p$','interpreter','latex')
subplot(4,4,2), plot(x,rmean(:,18));%Plot of recursive mean of r_pi
title('Recursive mean of $r_{\pi}$','interpreter','latex')
subplot(4,4,3), plot(x,rmean(:,19));%Plot of recursive mean of rho
title('Recursive mean of $\rho$','interpreter','latex')
subplot(4,4,4), plot(x,rmean(:,20));%Plot of recursive mean of r_y
title('Recursive mean of $r_y$','interpreter','latex')
subplot(4,4,5), plot(x,rmean(:,21));%Plot of recursive mean of r_delta_y
title('Recursive mean of $r_{\Delta_y}$','interpreter','latex')
subplot(4,4,6), plot(x,rmean(:,22));%Plot of recursive mean of pi_bar
title('Recursive mean of $\bar{\pi}$','interpreter','latex')
subplot(4,4,7), plot(x,rmean(:,23));%Plot of recursive mean of beta_bar
title('Recursive mean of $\bar{\beta}$','interpreter','latex')
subplot(4,4,8), plot(x,rmean(:,24));%Plot of recursive mean of l_bar
title('Recursive mean of $\bar{l}$','interpreter','latex')
subplot(4,4,9), plot(x,rmean(:,25));%Plot of recursive mean of gamma_bar
title('Recursive mean of $\bar{\gamma}$','interpreter','latex')
subplot(4,4,10), plot(x,rmean(:,26));%Plot of recursive mean of alpha
title('Recursive mean of $\alpha$','interpreter','latex')
subplot(4,4,11), plot(x,rmean(:,27));%Plot of recursive mean of rho_a
title('Recursive mean of $\rho_a$','interpreter','latex')
subplot(4,4,12), plot(x,rmean(:,28));%Plot of recursive mean of rho_b
title('Recursive mean of $\rho_b$','interpreter','latex')
subplot(4,4,13), plot(x,rmean(:,29));%Plot of recursive mean of rho_g
title('Recursive mean of $\rho_g$','interpreter','latex')
subplot(4,4,14), plot(x,rmean(:,30));%Plot of recursive mean of rho_i
title('Recursive mean of $\rho_i$','interpreter','latex')
subplot(4,4,15), plot(x,rmean(:,31));%Plot of recursive mean of rho_r
title('Recursive mean of $\rho_r$','interpreter','latex')
subplot(4,4,16), plot(x,rmean(:,32));%Plot of recursive mean of rho_p
title('Recursive mean of $\rho_p$','interpreter','latex')




figure(9)

subplot(4,4,1), plot(x,rmean(:,33));%Plot of recursive mean of rho_w
title('Recursive mean of $\rho_w$','interpreter','latex')
subplot(4,4,2), plot(x,rmean(:,34));%Plot of recursive mean of gn_A
title('Recursive mean of $gna$','interpreter','latex')
subplot(4,4,3), plot(x,rmean(:,35));%Plot of recursive mean of gn_B
title('Recursive mean of $gnb$','interpreter','latex')
subplot(4,4,4), plot(x,rmean(:,36));%Plot of recursive mean of omega_A
title('Recursive mean of $\omega_A$','interpreter','latex')


end

