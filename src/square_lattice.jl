mutable struct SquareLattice
    Lx::Int
    Ly::Int
    N::Int
    velocities::Array{Float64, 2}
    positions::Array{Array{Float64, 1}, 1}
end

# Default constructor
SquareLattice(Lx::Int, Ly::Int, N::Int) = begin
    positions = [[i, j] for i in 1:Lx-1 for j in 1:Ly-1]
    velocities = zeros(Float64, N, 2)
    return SquareLattice(Lx, Ly, N, velocities, positions)
end

# Constructor with initial velocities based on Maxwell-Boltzmann distribution
SquareLattice(Lx::Int, Ly::Int, N::Int, temperature::Float64) = begin
    positions = [[i, j] for i in 1:Lx-1 for j in 1:Ly-1]
    
    # Initialize velocities based on Maxwell-Boltzmann distribution
    mass = 1.0  # Assuming unit mass for simplicity
    velocities = sqrt(temperature) * randn(N, 2)
    
    # Ensure conservation of total linear momentum
    total_momentum = sum(velocities, dims=1)
    velocities .-= total_momentum / N
    
    # Assign opposite velocity to the last particle to make total linear momentum zero
    velocities[end, 1] = -sum(velocities[1:end-1, 1])
    velocities[end, 2] = -sum(velocities[1:end-1, 2])
    println("Total Linear Momentum: ", sum(velocities, dims=1))
    
    return SquareLattice(Lx, Ly, N, velocities, positions)
end

lattice = SquareLattice(9,9, 64,1.0)

initial_positions, initial_velocities = lattice.positions, lattice.velocities
println("Initial Positions:")

display(initial_positions)
println("Initial Velocities:")
display(initial_velocities)

println(typeof(initial_positions))
println(typeof(initial_velocities))