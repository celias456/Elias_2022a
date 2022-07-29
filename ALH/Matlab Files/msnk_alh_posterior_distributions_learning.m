function [] = msnk_alh_posterior_distributions_learning(theta)
%Plots the posterior distributions for the medium scale new keynesian model
%with adaptive learning

figure(1)
subplot(3,1,1),histogram(theta(:,34));%Histogram of the posterior dist. of gn_A
title('$gn_A$','interpreter','latex')
subplot(3,1,2),histogram(theta(:,35));%Histogram of the posterior dist. of gn_B
title('$gn_B$','interpreter','latex')
subplot(3,1,3),histogram(theta(:,36));%Histogram of the posterior dist. of omega_A
title('$\omega_A$','interpreter','latex')


end

