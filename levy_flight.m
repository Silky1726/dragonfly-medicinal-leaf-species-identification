function levy_step = levy_flight(Lambda)
    sigma = (gamma(1 + Lambda) * sin(pi * Lambda / 2) / ...
            (gamma((1 + Lambda) / 2) * Lambda * 2^((Lambda - 1) / 2)))^(1 / Lambda);
    u = randn * sigma;
    v = randn;
    levy_step = u / abs(v)^(1 / Lambda);
end
