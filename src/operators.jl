@inline function kernel1d(k, d, h)
    return -h*1im*exp(1im*k*d)/(2*k)
end

@inline function kernel2d(k, d, h)
    if d < h/2
        return -(1im*h/k)*hankelh1(1, k*h/2) + 4/(pi*k^2)
    end
    return -h^2*(.25im)*hankelh1(0, k*d) # TODO: Maybe change to a better approximation ? 
end

# Green's function
function generate_op(grid_out, grid_in; k, h, kernel::Function)
    G = zeros(ComplexF64, length(grid_out), length(grid_in))

    for (j, y_loc) in enumerate(grid_in)
        for (i, x_loc) in enumerate(grid_out)
            d = norm(x_loc .- y_loc) # ||x_i - x_j||
            G[i, j] = kernel(k, d, h)
        end
    end

    return G
end
 
function source_field(grid_out, grid_in, src; k, h, kernel::Function)
    b = zeros(ComplexF64, length(grid_out))

    for (i, x_loc) in enumerate(grid_out)
        for (j, y_loc) in enumerate(grid_in)
            d = norm(x_loc .- y_loc) # ||x_i - x_j||
            b[i] += kernel(k, d, h)*src[j]
        end
    end

    return b
end