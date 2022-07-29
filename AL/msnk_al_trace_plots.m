function [] = msnk_al_trace_plots(theta,log_posterior)
%Plots the post burn-in draws of each parameter for the medium scale new
%keynesian model with adaptive learning

n = size(theta,1);
x = (1:n);

figure(4)

subplot(4,4,1), plot(x,theta(:,1));%Plot of draws of sigma_a
title('Draws of $\sigma_a$','interpreter','latex')
subplot(4,4,2), plot(x,theta(:,2));%Plot of draws of sigma_b
title('Draws of $\sigma_b$','interpreter','latex')
subplot(4,4,3), plot(x,theta(:,3));%Plot of draws of sigma_g
title('Draws of $\sigma_g$','interpreter','latex')
subplot(4,4,4), plot(x,theta(:,4));%Plot of draws of sigma_i
title('Draws of $\sigma_i$','interpreter','latex')
subplot(4,4,5), plot(x,theta(:,5));%Plot of draws of sigma_r
title('Draws of $\sigma_r$','interpreter','latex')
subplot(4,4,6), plot(x,theta(:,6));%Plot of draws of sigma_p
title('Draws of $\sigma_p$','interpreter','latex')
subplot(4,4,7), plot(x,theta(:,7));%Plot of draws of sigma_w
title('Draws of $\sigma_w$','interpreter','latex')
subplot(4,4,8), plot(x,theta(:,8));%Plot of draws of varphi
title('Draws of $\varphi$','interpreter','latex')
subplot(4,4,9), plot(x,theta(:,9));%Plot of draws of sigma_c
title('Draws of $\sigma_c$','interpreter','latex')
subplot(4,4,10), plot(x,theta(:,10));%Plot of draws of lambda
title('Draws of $\lambda$','interpreter','latex')
subplot(4,4,11), plot(x,theta(:,11));%Plot of draws of xi_w
title('Draws of $\xi_w$','interpreter','latex')
subplot(4,4,12), plot(x,theta(:,12));%Plot of draws of sigma_l
title('Draws of $\sigma_l$','interpreter','latex')
subplot(4,4,13), plot(x,theta(:,13));%Plot of draws of xi_p
title('Draws of $\xi_p$','interpreter','latex')
subplot(4,4,14), plot(x,theta(:,14));%Plot of draws of iota_w
title('Draws of $\iota_w$','interpreter','latex')
subplot(4,4,15), plot(x,theta(:,15));%Plot of draws of iota_p
title('Draws of $\iota_p$','interpreter','latex')
subplot(4,4,16), plot(x,theta(:,16));%Plot of draws of psi
title('Draws of $\psi$','interpreter','latex')



figure(5)

subplot(4,4,1), plot(x,theta(:,17));%Plot of draws of phi_p
title('Draws of $\phi_p$','interpreter','latex')
subplot(4,4,2), plot(x,theta(:,18));%Plot of draws of r_pi
title('Draws of $r_{\pi}$','interpreter','latex')
subplot(4,4,3), plot(x,theta(:,19));%Plot of draws of rho
title('Draws of $\rho$','interpreter','latex')
subplot(4,4,4), plot(x,theta(:,20));%Plot of draws of r_y
title('Draws of $r_y$','interpreter','latex')
subplot(4,4,5), plot(x,theta(:,21));%Plot of draws of r_delta_y
title('Draws of $r_{\Delta_y}$','interpreter','latex')
subplot(4,4,6), plot(x,theta(:,22));%Plot of draws of pi_bar
title('Draws of $\bar{\pi}$','interpreter','latex')
subplot(4,4,7), plot(x,theta(:,23));%Plot of draws of beta_bar
title('Draws of $\bar{\beta}$','interpreter','latex')
subplot(4,4,8), plot(x,theta(:,24));%Plot of draws of l_bar
title('Draws of $\bar{l}$','interpreter','latex')
subplot(4,4,9), plot(x,theta(:,25));%Plot of draws of gamma_bar
title('Draws of $\bar{\gamma}$','interpreter','latex')
subplot(4,4,10), plot(x,theta(:,26));%Plot of draws of alpha
title('Draws of $\alpha$','interpreter','latex')
subplot(4,4,11), plot(x,theta(:,27));%Plot of draws of rho_a
title('Draws of $\rho_a$','interpreter','latex')
subplot(4,4,12), plot(x,theta(:,28));%Plot of draws of rho_b
title('Draws of $\rho_b$','interpreter','latex')
subplot(4,4,13), plot(x,theta(:,29));%Plot of draws of rho_g
title('Draws of $\rho_g$','interpreter','latex')
subplot(4,4,14), plot(x,theta(:,30));%Plot of draws of rho_i
title('Draws of $\rho_i$','interpreter','latex')
subplot(4,4,15), plot(x,theta(:,31));%Plot of draws of rho_r
title('Draws of $\rho_r$','interpreter','latex')
subplot(4,4,16), plot(x,theta(:,32));%Plot of draws of rho_p
title('Draws of $\rho_p$','interpreter','latex')




figure(6)

subplot(4,4,1), plot(x,theta(:,33));%Plot of draws of rho_w
title('Draws of $\rho_w$','interpreter','latex')
subplot(4,4,2), plot(x,theta(:,34));%Plot of draws of gn_A
title('Draws of $gna$','interpreter','latex')
subplot(4,4,3), plot(x,log_posterior);%Plot of draws of log posterior
title('Draws of log posterior','interpreter','latex')



end

