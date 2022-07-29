function [] = msnk_al_posterior_distributions(theta)
%Plots the posterior distributions for the medium scale new keynesian model
%with adaptive learning

figure(1)

subplot(4,4,1),hist(theta(:,1));%Histogram of the posterior dist. of sigma_a
title('Posterior of $\sigma_a$','interpreter','latex')
subplot(4,4,2),hist(theta(:,2));%Histogram of the posterior dist. of sigma_b
title('Posterior of $\sigma_b$','interpreter','latex')
subplot(4,4,3),hist(theta(:,3));%Histogram of the posterior dist. of sigma_g
title('Posterior of $\sigma_g$','interpreter','latex')
subplot(4,4,4),hist(theta(:,4));%Histogram of the posterior dist. of sigma_i
title('Posterior of $\sigma_i$','interpreter','latex')
subplot(4,4,5),hist(theta(:,5));%Histogram of the posterior dist. of sigma_r
title('Posterior of $\sigma_r$','interpreter','latex')
subplot(4,4,6),hist(theta(:,6));%Histogram of the posterior dist. of sigma_p
title('Posterior of $\sigma_p$','interpreter','latex')
subplot(4,4,7),hist(theta(:,7));%Histogram of the posterior dist. of sigma_w
title('Posterior of $\sigma_w$','interpreter','latex')
subplot(4,4,8),hist(theta(:,8));%Histogram of the posterior dist. of varphi
title('Posterior of $\varphi$','interpreter','latex')
subplot(4,4,9),hist(theta(:,9));%Histogram of the posterior dist. of sigma_c
title('Posterior of $\sigma_c$','interpreter','latex')
subplot(4,4,10),hist(theta(:,10));%Histogram of the posterior dist. of lambda
title('Posterior of $\lambda$','interpreter','latex')
subplot(4,4,11),hist(theta(:,11));%Histogram of the posterior dist. of xi_w
title('Posterior of $\xi_w$','interpreter','latex')
subplot(4,4,12),hist(theta(:,12));%Histogram of the posterior dist. of sigma_l
title('Posterior of $\sigma_l$','interpreter','latex')
subplot(4,4,13),hist(theta(:,13));%Histogram of the posterior dist. of xi_p
title('Posterior of $\xi_p$','interpreter','latex')
subplot(4,4,14),hist(theta(:,14));%Histogram of the posterior dist. of iota_w
title('Posterior of $\iota_w$','interpreter','latex')
subplot(4,4,15),hist(theta(:,15));%Histogram of the posterior dist. of iota_p
title('Posterior of $\iota_p$','interpreter','latex')
subplot(4,4,16),hist(theta(:,16));%Histogram of the posterior dist. of psi
title('Posterior of $\psi$','interpreter','latex')

figure(2)

subplot(4,4,1),hist(theta(:,17));%Histogram of the posterior dist. of phi_p
title('Posterior of $\phi_p$','interpreter','latex')
subplot(4,4,2),hist(theta(:,18));%Histogram of the posterior dist. of r_pi
title('Posterior of $r_{\pi}$','interpreter','latex')
subplot(4,4,3),hist(theta(:,19));%Histogram of the posterior dist. of rho
title('Posterior of $\rho$','interpreter','latex')
subplot(4,4,4),hist(theta(:,20));%Histogram of the posterior dist. of r_y
title('Posterior of $r_y$','interpreter','latex')
subplot(4,4,5),hist(theta(:,21));%Histogram of the posterior dist. of r_delta_y
title('Posterior of $r_{\Delta_y}$','interpreter','latex')
subplot(4,4,6),hist(theta(:,22));%Histogram of the posterior dist. of pi_bar
title('Posterior of $\bar{\pi}$','interpreter','latex')
subplot(4,4,7),hist(theta(:,23));%Histogram of the posterior dist. of beta_bar
title('Posterior of $\bar{\beta}$','interpreter','latex')
subplot(4,4,8),hist(theta(:,24));%Histogram of the posterior dist. of l_bar
title('Posterior of $\bar{l}$','interpreter','latex')
subplot(4,4,9),hist(theta(:,25));%Histogram of the posterior dist. of gamma_bar
title('Posterior of $\bar{\gamma}$','interpreter','latex')
subplot(4,4,10),hist(theta(:,26));%Histogram of the posterior dist. of alpha
title('Posterior of $\alpha$','interpreter','latex')
subplot(4,4,11),hist(theta(:,27));%Histogram of the posterior dist. of rho_a
title('Posterior of $\rho_a$','interpreter','latex')
subplot(4,4,12),hist(theta(:,28));%Histogram of the posterior dist. of rho_b 
title('Posterior of $\rho_b$','interpreter','latex')
subplot(4,4,13),hist(theta(:,29));%Histogram of the posterior dist. of rho_g 
title('Posterior of $\rho_g$','interpreter','latex')
subplot(4,4,14),hist(theta(:,30));%Histogram of the posterior dist. of rho_i
title('Posterior of $\rho_i$','interpreter','latex')
subplot(4,4,15),hist(theta(:,31));%Histogram of the posterior dist. of rho_r
title('Posterior of $\rho_r$','interpreter','latex')
subplot(4,4,16),hist(theta(:,32));%Histogram of the posterior dist. of rho_p
title('Posterior of $\rho_p$','interpreter','latex')


figure(3)
subplot(4,4,1),hist(theta(:,33));%Histogram of the posterior dist. of rho_w
title('Posterior of $\rho_w$','interpreter','latex')
subplot(4,4,2),hist(theta(:,34));%Histogram of the posterior dist. of gn_A
title('Posterior of $gna$','interpreter','latex')


end

