using Optim

# Minimising a multivariate function
f(u) = u[1]^2 + u[2]^2 # Let f = f(u(x,t)) be the function we want to minimise
u_0 = [0.0, 0.1] # and u_0 be the initial guess

result = optimize(f, u_0) # Defaults to Nelder-Mead method (direct method but can converge to non-stationary points)
optimize(f, u_0, LBFGS()) # Can also specify other methods, such as LBFGS (quasi-Newton method, approximates inverse Hessian and stores it implicitly)

# Access results
method = Optim.summary(result)
minimum = Optim.minimum(result)
minimizer = Optim.minimizer(result)

# Minimising a univariate function
f(x) = x^2 + x + 2
start_interval, end_interval = -5, 5
result = optimize(f, start_interval, end_interval) # Defaults to Brent's method, can specify GoldenSection()
Optim.lower_bound(result)
Optim.upper_bound(result)
