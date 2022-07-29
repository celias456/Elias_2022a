function [mu_t] = recursive_mean(data)
%This function computes the recursive mean for a matrix of data.

%Input:
    %data:  each column represents a variable and each row represents an
    %observation of the corresponding variable

%Output:
    %mu_t: each column represents a variable, each row is the recursive
    %average for the observation corresponding with that row

%References: 
%Metropolis-Hastings and DSGE models slides Slides_DSGE_MH.pdf); slide numbers 58 through 62
%https://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Average.htm

number_observations = size(data,1);
number_variables = size(data,2);

mu_t = zeros(number_observations,number_variables);
mu_t(1,:) = data(1,:);

for t = 2:number_observations
    mu_t(t,:) = ((t-1)/t)*mu_t(t-1,:) + (1/t)*data(t,:);
end

end

