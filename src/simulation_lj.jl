using Random
using Statistics
using Plots
using LinearAlgebra
using Optim
include("square_lattice.jl")
function lj_force(r::Float64)
    epsilon = 1.0
    sigma = 1.0
    r_inv = 1.0 / r
    r6 = r_inv^6
    r12 = r6^2
    force = 24 * epsilon * (2 * r12 - r6) * r_inv
    return force
end


function radial_distribution_function(positions, Lx, Ly, nbins=50)
    dr = min(Lx, Ly) / (2 * nbins)
    bins = collect(0:dr:min(Lx, Ly)/2)
    rdf = zeros(nbins)

    for i in 1:N
        for j in i+1:N
            rij = positions[j, :] - positions[i, :]
            rij .= rij .- Lx * round.(Int, rij ./ Lx)
            r = norm(rij)

            if r < min(Lx, Ly) / 2
                bin_index = floor(Int, r / dr) + 1
                rdf[bin_index] += 2  # Count each particle pair only once
            end
        end
    end

    rdf /= (N * (N - 1) / 2) * π * (bins[2] - bins[1])^2

    return bins, rdf
end

function speed_probability_density(velocities, nbin=32, dvel=0.16)
    speeds = sqrt.(sum(velocities .^ 2, dims=2))
    max_speed = maximum(speeds)
    bins = collect(0:dvel:max_speed)

    # Initialize histogram array
    hist = zeros(nbin)

    for v in speeds
        ibin = floor(Int, v[1] / dvel) + 1
        if ibin > nbin
            ibin = nbin
        end
        hist[ibin] += 1  # Adjust index to start from 1
    end

    # Normalize by the number of particles and the number of MD time steps
    hist /= (N * total_steps * dvel)

    return bins, hist
end

function maxwell_boltzmann_distribution(v, A, m, k_B, T)
    return A * exp(-m * v^2 / (2 * k_B * T)) * v
end

function apply_periodic_boundary_conditions!(positions, L)
    for i in 1:N
        for j in 1:2
            positions[i, j] %= L
        end
    end
end
function verlet_integration!(positions, velocities, forces, L, tau)
    positions .+= tau * velocities + 0.5 * tau^2 * forces
    apply_periodic_boundary_conditions!(positions, L)
    new_forces = calculate_forces(positions, L)
    velocities .+= 0.5 * tau * (forces + new_forces)
    return positions, velocities, new_forces
end

function calculate_forces(positions, L)
    forces = zeros(N, 2)
    for i in 1:N
        for j in i+1:N
            rij = [positions[j][1] - positions[i][1], positions[j][2] - positions[i][2]]
            rij .-= L .* round.(Int, rij ./ L) # Apply periodic boundary conditions
            r = norm(rij)
            force_magnitude = lj_force(r)
            force = force_magnitude * rij / r
            forces[i] = force[1]
            forces[j] -= force[2]
        end
    end
    display(forces)
    println(typeof(forces))
    return forces
end


function calculate_energy(positions, velocities, L)
    kinetic_energy = 0.5 * sum(sum(velocities .^ 2, dims=2))
    potential_energy = 0.0
    for i in 1:N
        for j in i+1:N
            rij = [positions[j][1] - positions[i][1], positions[j][2] - positions[i][2]]
            # Apply periodic boundary conditions
            rij .-= L .* round.(Int, rij ./ L)
            
            r = norm(rij)
            potential_energy += 4 * ((1 / r)^12 - (1 / r)^6)
        end
    end
    return kinetic_energy, potential_energy
end



function md_simulation(positions, velocities, L, tau, steps,  temperature_interval, equilibrium_steps, nbin=50, dr=0.025, gcum=zeros(nbin))
    temperature_values = []
    pressure_values = []
    mean_temperature_values = []
    mean_pressure_values = []

    for step in 1:equilibrium_steps
        forces = calculate_forces(positions, L)
        positions, velocities, forces = verlet_integration!(positions, velocities, forces, L, tau)
    

        if step % steps == 0
            total_momentum = sum(velocities, dims=1)
            velocities .-= total_momentum / N
        end

        kinetic_energy, potential_energy = calculate_energy(hcat(positions...), velocities, Lx)

        temperature = (2 / (3 * N)) * kinetic_energy
        pressure = (N / (Lx * Ly)) * kinetic_energy + (1 / (Lx * Ly)) * sum(forces .* velocities)

        push!(temperature_values, temperature)
        push!(pressure_values, pressure)
        # Call compute_g! function
        compute_g!(step, positions, N, Lx, Ly, gcum, nbin, dr)
        if step % temperature_interval == 0
            mean_temperature = mean(temperature_values[end-temperature_interval+1:end])
            mean_pressure = mean(pressure_values[end-temperature_interval+1:end])

            push!(mean_temperature_values, mean_temperature)
            push!(mean_pressure_values, mean_pressure)

            println("Step $step: Mean Temperature = $mean_temperature, Mean Pressure = $mean_pressure")
        end
    end
    return positions, velocities, forces
end
# Example usage
# Example usage
N = 64
Lx = 9
Ly = 9
L=Lx
tau = 0.01
steps = 200
total_steps = 200

# Generate initial positions and velocities
lattice = SquareLattice(9, 9, 64,1.0)

positions, velocities = lattice.positions, lattice.velocities
positions, velocities, forces = md_simulation(positions, velocities, L, tau, steps, 10, 1000)
kinetic_energy, potential_energy = calculate_energy(positions, velocities, L)
total_energy = kinetic_energy + potential_energy
temperature = (2 / (3 * N)) * kinetic_energy

println("Equilibrium reached in steps.")
println("Kinetic Energy: $kinetic_energy")
println("Potential Energy: $potential_energy")
println("Total Energy: $total_energy")
println("Temperature: $temperature")

# Visualization (plotting particle positions)
using Plots
scatter([pos[1] for pos in positions], [pos[2] for pos in positions], xlabel="X", ylabel="Y", title="Particle Positions", legend=false)
