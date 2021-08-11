function CI95 = ConfInt95(data, dim)
% calculates the 95% confidence interval of input data along dimension dim
% So, if observations are rows and variables are columns, then call
% ConfInf95(data, 1)
stdev = std(data, 0, dim);
N = size(data, dim);
CI95 = stdev/sqrt(N)*tinv(0.975,N-1);

end