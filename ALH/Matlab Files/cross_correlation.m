function [ccor] = cross_correlation(series1,series2,lag_length)
%This function computes the cross-correlation between two series

%Input:
    %series1: a column vector of data
    %series2: a column vector of data
    %lag_length: the lag length
    
%Output:
    %ccor: the cross-correlation coefficient

T = size(series1,1);

lag = abs(lag_length);

if lag_length <= 0
    series1_mod = series1(1+lag:end);
    series2_mod = series2(1:T-lag);
else
    series1_mod = series1(1:T-lag);
    series2_mod = series2(1+lag:end);
end

ccor_matrix = corrcoef(series1_mod,series2_mod);

ccor = ccor_matrix(2,1);

end

