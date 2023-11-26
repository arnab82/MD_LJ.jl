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


# Example:ß
N = 10 
Lx = Ly = 10.0  
nbin = 1000 
gcum = zeros(Float64, nbin)  
ncorrel = 1000 
dr = 0.025 


normalize_g(ncorrel, N, Lx, Ly, gcum, dr,nbin)
x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]  
y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
compute_g(ncorrel, x, y, N, Lx, Ly, gcum, nbin, dr)