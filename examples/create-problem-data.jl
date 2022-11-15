#=
# Mode converter design setup
We show how to setup the problem data to design a mode converter, as we do in
[Bounds on Efficiency Metrics in Photonics](https://arxiv.org/abs/2204.05243).
=#

## Import the packages
using WaveOperators
using Plots


#=
## Constructing the design region
First we set basic parameters of our design region and the
simulation grid.
=#
n = 50                                  # Number of gridpoints per unit length
width, height = 2, 1                    # Width and height of the domain
k = 2π                                  # Wavenumber of the domain

## Construct grid
g = Grid(height, width, 1/n, k);


#=
Next, we add the waveguide slab, centered, of length `height/2` and width 
`height/3`. We set the material contrast to be 5.
=#
contrast = 5
waveguide = Slab(height/2, height/3)
set_contrast!(g, waveguide, contrast)


#=
Finally we define the design region itself, centered on the domain, of size `2/3*height`
by `height`.
=#
d_height = 2*height/3
d_width = width/2
x_pos = height/2 - d_height/2
y_pos = width/2 - d_width/2

design_region = Rectangle((x_pos, y_pos), (d_height, d_width))

## Set the maximum allowable contrast for the design region
set_contrast!(g, design_region, contrast)
heatmap(g.contrast, title="Design Region")


#=
## Adding an input field
With our design region defined, we now add the input mode to the design region.
=#
mode_input_position = 0
input_mode = 1
input_line = VerticalLine(mode_input_position)
mode = add_mode!(g, input_line, input_mode);

#=
We can visualize the amplitude of this mode.
=#
modes = compute_modes(g, input_line)
mode_in = modes.vectors[:, input_mode]
plot(abs.(mode_in), grid=false, legend=false,lw=3, color=:black)
#=
And the resulting field. (Note, we haven't designed anything yet!)
=#
solution = solve(g)
heatmap(abs.(solution))


#=
For a nicer looking plot (like the ones we used in 
[the paper](https://arxiv.org/abs/2204.05243)), we can interpolate:
=#
using Interpolations: LinearInterpolation
function interpolate(img; factor=8, ylims=(0,2), xlims=(0,2))
    xx = range(xlims..., size(img, 1))
    yy = range(ylims..., size(img, 2))
    itp = LinearInterpolation((xx,yy), img)
    x2 = range(xlims..., size(img, 1)*factor)
    y2 = range(ylims..., size(img, 2)*factor)
    return [itp(x, y) for x in x2, y in y2]
end
heatmap(interpolate(abs.(solution)))


#=
## Defining the objective
For a mode converter, we want to measure correlation with some predefined 
output mode at a target (in this case, the RHS of the domain.) See section 1.1 
of [our paper](https://arxiv.org/abs/2204.05243). 
=#

## Make a target line at the right hand side of the domain
target_line = VerticalLine(width)

## Generate the Green's function for the problem
g_functions = generate_G(g, design_region, target_line)

## Compute the desired output mode vector
modes = compute_modes(g, target_line)
output_mode = 2
c = modes.vectors[:, output_mode];

#=
The desired output mode is plotted below:
=#
out_fig = plot(abs.(c), grid=false, legend=false, lw=3, color=:black)


#=
## Making a design
We can use the tools in `PhysicalBounds.jl` to create a design. Below, we show
how to use the optimal value of the optimization problem to set the design from
our variable vector $\theta$.
=#
# Usually θ would be replaced by a real design. Here we use a 'random' design.
# See PhysicalBounds.jl for a complete example.
designidx = getindices(g, design_region)
θ = rand(length(designidx))

# Set the contrast to be the design's
g.contrast[designidx] .*= θ
heatmap(g.contrast, title="Design Region")

# Solve for the new design's fields
solution_new = solve(g)
heatmap(abs.(solution_new))
