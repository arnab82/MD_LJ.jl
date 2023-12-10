using Random
using Statistics
using Plots
using LinearAlgebra
using Optim
function lj_force(r::Float64)
    epsilon = 1.0
    sigma = 1.0
    r_inv = 1.0 / r
    r6 = r_inv^6
    r12 = r6^2
    force = 24 * epsilon * (2 * r12 - r6) * r_inv
    return force
end

function maxwell_boltzmann_distribution(v, A, m, k_B, T)
    return A * exp(-m * v^2 / (2 * k_B * T)) * v
end

function apply_periodic_boundary_conditions!(positions, L)
    for i in 1:N
        if positions[i][1] > L
            positions[i][1] -= L
        elseif positions[i][1] < 0
            positions[i][1] += L
        end
        if positions[i][2] > L
            positions[i][2] -= L
        elseif positions[i][2] < 0
            positions[i][2] += L
        end
    end
end
function verlet_integration!(positions, velocities, forces, L, tau,fc)
    N = length(positions)

    new_positions = deepcopy(positions)
    # println(typeof(new_positions))
    for i in 1:N
        for dim in 1:2
            # Update positions
            new_positions[i][dim] = new_positions[i][dim] + tau * velocities[i, dim] + 0.5 * tau^2 * forces[i][dim]

            # Apply periodic boundary conditions
            if new_positions[i][dim] > L / 2
                new_positions[i][dim] = (L - abs(new_positions[i][dim])) * (-new_positions[i][dim]) / abs(new_positions[i][dim])
            end
        end
    end

    apply_periodic_boundary_conditions!(new_positions, L)

    # Update forces and velocities
    new_forces = calculate_forces(new_positions, L, fc)
    new_velocities = deepcopy(velocities)  # Make a copy to avoid modifying the original velocities array
    for i in 1:N
        for dim in 1:2
            # println(i)
            # println(dim)
            # println((forces[i][dim] .+ new_forces[i][dim]))
            # println(new_velocities[i, dim])
            # Update velocities
            new_velocities[i, dim] += 0.5 * tau * (forces[i][dim] .+ new_forces[i][dim])
        end
    end

    return new_positions, new_velocities, new_forces
end



function calculate_forces(positions, L, fc, rcut=2.5)
    forces = []
    N = length(positions)
    for i in 1:N
        fx = 0.0
        fy = 0.0
        for j in i+1:N
            rij = [positions[j][1] - positions[i][1], positions[j][2] - positions[i][2]]
            if rij[1] > L / 2
                rij[1] = (L - abs(rij[1])) * (-rij[1]) / abs(rij[1])
            end
            if rij[2] > L / 2
                rij[2] = (L - abs(rij[2])) * (-rij[2]) / abs(rij[2])
            end
            # rij .-= L .* round.(Int, rij ./ L) # Apply periodic boundary conditions
            r = norm(rij)
            if r <= rcut
                force_magnitude = lj_force(r)
                force = force_magnitude * rij / r
                force[1] -= fc
                force[2] -= fc
                # println("debug")
                # println(force)
                fx += force[1]
                fy += force[2]

            end
        end
        push!(forces, (fx, fy))
    end

    # display(forces)
    # println(typeof(forces))
    return forces
end


function calculate_energy(positions, velocities, L, Vc, fc, rc)
    kinetic_energy = 0.5 * sum(sum(velocities .^ 2, dims=2))
    potential_energy = 0.0
    for i in 1:N
        for j in i+1:N
            rij = [positions[j][1] - positions[i][1], positions[j][2] - positions[i][2]]
            # Apply periodic boundary conditions
            if rij[1] > L / 2
                rij[1] = (L - abs(rij[1])) * (-rij[1]) / abs(rij[1])
            end
            if rij[2] > L / 2
                rij[2] = (L - abs(rij[2])) * (-rij[2]) / abs(rij[2])
            end
            # rij .-= L .* round.(Int, rij ./ L)

            r = norm(rij)
            V = 4 * ((1 / r)^12 - (1 / r)^6)
            Vm = V - Vc + r * fc - rc * fc
            potential_energy += Vm
        end
    end
    return kinetic_energy, potential_energy
end



function md_simulation(positions, velocities, L, tau, steps, temperature_interval, equilibrium_steps, fc, Vc, rc,tolerance=1e-3, stability_threshold=5)
    temperature_values = []
    pressure_values = []
    mean_temperature_values = []
    mean_pressure_values = []
    N = length(positions)
    V= L^2
    dV = 1.0
    kB = 1.0
    prev_mean_temperature = Inf
    prev_mean_pressure = Inf
    stability_counter = 0
    
    for step in 1:equilibrium_steps
        println("Step $step")
        forces = calculate_forces(positions, L, fc)
        positions, velocities, forces = verlet_integration!(positions, velocities, forces, L, tau,fc)

        if step % steps == 0
            total_momentum = sum(velocities, dims=1)
            velocities .-= total_momentum / N
        end

        kinetic_energy, potential_energy = calculate_energy(positions, velocities, Lx, Vc, fc, rc)

        temperature = (2 / (3 * N)) * kinetic_energy
        pairwise_distances = [norm(positions[i] - positions[j]) for i in 1:N for j in 1:N if i != j]
        pairwise_forces = Float64[]

        for i in 1:N
            for j in 1:N
                if i != j
                    force = 0.0
                    for k in 1:2
                        force += forces[i][k] * (positions[i][k] - positions[j][k])
                    end
                    push!(pairwise_forces, force / pairwise_distances[(i-1)*(N-1)+(j-1)])
                end
            end
        end
        pressure_term = 0.5 * sum(pairwise_distances .* pairwise_forces)

        pressure = (N / V) * kB * temperature + (1 / dV) * pressure_term
        push!(temperature_values, temperature)
        push!(pressure_values, pressure)
        if step % temperature_interval == 0
            mean_temperature = mean(temperature_values[end-temperature_interval+1:end])
            mean_pressure = mean(pressure_values[end-temperature_interval+1:end])

            push!(mean_temperature_values, mean_temperature)
            push!(mean_pressure_values, mean_pressure)

            println("Step $step: Mean Temperature = $mean_temperature, Mean Pressure = $mean_pressure")
            
            # Check for stability
            if abs(mean_temperature - prev_mean_temperature) < tolerance && abs(mean_pressure - prev_mean_pressure) < tolerance
                stability_counter += 1
                println("System is stable. Stability counter = $stability_counter")
            else
                stability_counter = 0
            end
            
            if stability_counter >= stability_threshold
                println("System has reached equilibrium. Stopping simulation.")
                break
            end
            
            prev_mean_temperature = mean_temperature
            prev_mean_pressure = mean_pressure
        end
    end
    return positions, velocities
end
