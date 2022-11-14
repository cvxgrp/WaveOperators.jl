import Base.CartesianIndices
export Grid, pos, solve, add_current!, add_mode!, generate_G, set_contrast!, compute_modes, getindices

struct Grid
    h::Float64
    k::Float64
    dims::Tuple{Int64, Int64}

    contrast::Matrix{Float64}
    input::Matrix{ComplexF64}
end

"""
    Grid(height::Real, width::Real, h::Real, k::Real)

Creates a `height` by `width` grid with discretization `h` and wavenumber `k`.
This grid will define our problem domain.
"""
function Grid(height::Real, width::Real, h::Real, k::Real)
    dims = (floor(Int64, height/h)+1, floor(Int64, width/h)+1)
    contrast = zeros(dims)
    input = zeros(ComplexF64, dims)

    return Grid(float(h), float(k), dims, contrast, input)
end

# Make Grid behave like a scalar for `.` broacasting purposes
Base.broadcastable(g::Grid) = Ref(g)

# Some grid methods
CartesianIndices(g::Grid) = CartesianIndices(g.contrast)
generate_op(g::Grid, grid_out, grid_in) = generate_op(grid_out, grid_in; g.k, g.h, kernel=kernel2d)
source_field(g::Grid, grid_out, grid_in, src) = source_field(grid_out, grid_in, src; g.k, g.h, kernel=kernel2d)

"""
    pos(g::Grid, x)

Returns the position of `x` in grid 
"""
function pos end

pos(g::Grid, x::N) where N <: Number = g.h*(x-1)
pos(g::Grid, x::CartesianIndex) = pos.(g, x.I)


"""
    getindices(g::Grid, s::Shape) where S <: Shape

Gets indices (wrt the grid `g`) of the shape `s`
"""
getindices(g::Grid, s::S) where S <: Shape = filter(p -> is_in(s, pos(g, p)), CartesianIndices(g))
pos(g::Grid, s::S) where S <: Shape = pos.(g, getindices(g, s))

function getindices(g::Grid, l::VerticalLine) 
    slice_idx = floor(Int64, l.x/g.h) + 1
    return CartesianIndices(g)[:, slice_idx]
end

function getindices(g::Grid, l::HorizontalLine) 
    slice_idx = floor(Int64, l.y/g.h) + 1
    return CartesianIndices(g)[slice_idx, :]
end

function getindices(g::Grid, slice::ShapeSlice{S}) where {S <: Shape}
    sliceidxs = getindices(g, slice.line)
    return filter(p -> is_in(slice.shape, pos(g, p)), sliceidxs)
end

"""
    set_contrast!(g::Grid, s::S, val) where S <: Shape

Sets the contrast of the shape `s` in grid `g` to be `val`.
"""
function set_contrast!(g::Grid, s::S, val) where S <: Shape
    for I in getindices(g, s)
        g.contrast[I] = val
    end
    
    return nothing
end

"""
    solve(g::Grid)

Computes the field on the grid `g` (an excitation should be added first).
"""
function solve(g::Grid)
    idx_contrast = findall(!iszero, g.contrast)
    idx_exterior = findall(iszero, g.contrast)

    pos_contrast = pos.(g, idx_contrast)
    pos_exterior = pos.(g, idx_exterior)

    output = zero(g.input)

    G_dd = generate_op(g, pos_contrast, pos_contrast)
    output[idx_contrast] .= (I + g.k^2*G_dd*Diagonal(g.contrast[idx_contrast])) \ g.input[idx_contrast]

    aux_field = -g.k^2*source_field(g, pos_exterior, pos_contrast, Diagonal(g.contrast[idx_contrast])*output[idx_contrast])
    output[idx_exterior] = aux_field + g.input[idx_exterior]

    return reshape(output, size(g.input))
end

"""
    compute_modes(g::Grid, x)

Returns the modes of the excitation on the grid `g`.

"""
function compute_modes(g::Grid, idxs::V) where V <: AbstractVector{<: CartesianIndex}
    slice_contrast = g.contrast[idxs]
    slice_pos = pos.(g, idxs)

    G = generate_op(slice_pos, slice_pos; g.k, g.h, kernel=kernel1d)

    return eigen!(-G\(I + g.k^2*G*Diagonal(slice_contrast)))
end

function compute_modes(g::Grid, s::S) where S <: Shape
    I = getindices(g, s)

    return compute_modes(g, I)
end

"""
    add_current!(g::Grid, s::S, input::N) where {S <: Shape, N <: Number}

Adds excitation `input` to the grid `g` at position `s`.
"""
function add_current!(g::Grid, s::S, input::N) where {S <: Shape, N <: Number}
    p = pos.(g, CartesianIndices(g))
    dims = size(g.contrast)

    s = pos.(g, getindices(g, s))

    g.input .+= reshape(source_field(g, p, s, input*ones(length(s))), dims...)

    return nothing 
end

"""
    add_current!(g::Grid, positions, input)

Adds excitation `input` to the grid `g` at positions positions.
"""
function add_current!(g::Grid, positions, input)
    idx = CartesianIndices(g)
    p = pos.(g, idx)

    dims = size(g.contrast)

    g.input .+= reshape(source_field(g, p, positions, input), dims...)

    return nothing
end

"""
    add_mode!(g::Grid, x, n_mode::Int)

Adds mode number `n_mode` to the grid `g` at position `x`, which may be an x
axis position or a [`VerticalLine`](@ref).
"""
function add_mode!(g::Grid, line::Line, n_mode::Int) 
    indices = getindices(g, line)

    modes = compute_modes(g, indices)

    slice_pos = pos.(g, indices)
    add_current!(g, slice_pos, modes.vectors[:, n_mode])

    return modes.vectors[:, n_mode]
end

add_mode!(g::Grid, x::Real, n_mode::Int=1) = add_mode!(g, VerticalLine(x), n_mode)

@doc raw"""
    function generate_G(g::Grid, design, target)

Generates the [`IntegralEquation`](@ref) 
(Green's function `G` and excitations `b`) for the equation
```math
z + G\mathbf{diag}(\theta)z = Gb,
```
The inputs `design` and `target` may be indices or shapes.
"""
function generate_G(g::Grid, design_idx::U, target_idx::V) where {U <: AbstractArray{<: CartesianIndex}, V <: AbstractArray{<: CartesianIndex}}
    # Function assumes that constrast in design_idx is set to the upper limits of contrast
    d_idx = Set(design_idx)
    t_idx = Set(target_idx)
    length(intersect(d_idx, t_idx)) != 0 && error("Overlapping design and target regions not yet supported")

    contrast_or_target(I) = !iszero(g.contrast[I]) || I in t_idx

    nonzero_contrast = filter(contrast_or_target, CartesianIndices(g))
    in_pos = pos.(g, nonzero_contrast)

    G_d = generate_op(g, in_pos, in_pos)
    G = g.k^2*G_d*Diagonal(g.contrast[nonzero_contrast])
    for (i, I) in enumerate(nonzero_contrast)
        if I in d_idx
            continue
        end
        G[i, i] += 1.0
    end

    out_idx = union(d_idx, t_idx)
    contrast_bool = [I in out_idx for I in nonzero_contrast]
    G_new, b_new, _ = schur_complement(G, g.input[nonzero_contrast], contrast_bool)

    new_indices = nonzero_contrast[contrast_bool]

    design_bool = [I in d_idx for I in new_indices]
    G_design, b_design, G_target, b_target, _ = schur_complement_full(G_new, b_new, design_bool)

    return IntegralEquation(
        G_design,
        G_target,
        b_design,
        b_target
    )
end

function generate_G(g::Grid, design_shape::S, target_shape::T) where {S <: Shape, T <: Shape}
    design_idx = getindices(g, design_shape)
    target_idx = getindices(g, target_shape)

    return generate_G(g, design_idx, target_idx)
end