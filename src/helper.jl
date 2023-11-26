using Plots

# Define the compute_g! subroutine
function compute_g(ncorrel, x, y, N, Lx, Ly, gcum, nbin, dr)
    function separation(a, b, L)
        return (a - b + L / 2) % L - L / 2
    end

    for i in 1:N - 1
        for j in i + 1:N
            dx = separation(x[i], x[j], Lx)
            dy = separation(y[i], y[j], Ly)
            r_squared = dx^2 + dy^2
            r = sqrt(r_squared)

            ibin = trunc(Int, r / dr)

            if ibin <= nbin
                gcum[ibin] += 1
            end
        end
    end
    return gcum
    ncorrel += 1
end

function normalize_g(ncorrel, N, Lx, Ly, gcum, dr, nbin)
    density = N / (Lx * Ly)
    rmax = min(Lx / 2, Ly / 2)
    normalization = density * ncorrel
    bin = 1
    r = 0.0

    result = []
    while r <= rmax
        area_shell = π * ((r + dr)^2 - r^2)
        gcum_bin = normalization * area_shell
        println("$r $(gcum[bin]) $(gcum_bin)")
        push!(result, (r, gcum[bin], gcum_bin))

        bin += 1
        r += dr
    end

    return result
end


function compute_equilibrium_prob(N, vx, vy, nbin, dvel, num_steps, m, kB, T)
    v_squared = zeros(Float64, N)
    prob = zeros(Float64, nbin)

    for step in 1:num_steps
        for i in 1:N
            v_squared[i] = vx[i]^2 + vy[i]^2
            v = sqrt(v_squared[i])

            ibin = min(nbin, Int(floor(v / dvel)) + 1)
            prob[ibin] += 1
        end
    end

    # Normalize by dividing by the number of particles and MD time steps
    prob /= (N * num_steps * dvel)
    v_values_random = collect(0:dvel:(nbin - 1) * dvel)
    # Compute the theoretical Maxwell-Boltzmann distribution
    A = (m / (2π * kB * T))^(3/2)
    theoretical_prob = A * exp.(-m * v_values_random.^2 / (2 * kB * T))
    return v_squared, prob, theoretical_prob,v_values_random
end

# Example 
N = 64 
vx = randn(N)  
vy = randn(N) 
nbin = 32
dvel = 0.16
num_steps = 1000
m = 1.0  # Particle mass
kB = 1.0  # Boltzmann constant
T = 300.0  # Temperature

# Call the function to get the data
v_squared, prob, theoretical_prob ,v_values= compute_equilibrium_prob(N, vx, vy, nbin, dvel, num_steps, m, kB, T)

# Plot the results
v_values = collect(0:dvel:(nbin - 1) * dvel)
plot(v_values, prob, xlabel="v", ylabel="P(v)", label="Simulation")
plot(v_values, theoretical_prob, label="Maxwell-Boltzmann Distribution")

Lx = Ly = 10.0  
nbin = 1000 
gcum = zeros(Float64, nbin)  
ncorrel = 1000 
dr = 0.025 


normalize_g(ncorrel, N, Lx, Ly, gcum, dr,nbin)
N=10
x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]  
y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
compute_g(ncorrel, x, y, N, Lx, Ly, gcum, nbin, dr)
