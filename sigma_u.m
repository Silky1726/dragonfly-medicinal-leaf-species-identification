function sigma=sigma_u(Lambda)
    sigma = (gamma(1 + Lambda) * sin(pi * Lambda / 2) / (gamma((1 + Lambda) / 2) * Lambda * 2^((Lambda - 1) / 2)))^(1 / Lambda);
end