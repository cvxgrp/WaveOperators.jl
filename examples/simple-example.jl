#=
# Simple Example
This is a very simple example showing how to create a design region
and add an input field.
=#

## Import the packages
using WaveOperators
using Plots

n = 50                                  # Number of gridpoints per unit length
width, height = 2, 1                    # Width and height of the domain
k = 2π                                  # Wavenumber of the domain

## Construct grid
g = Grid(height, width, 1/n, k);

#=
First, we add the waveguide and visualize the design region.
=#
waveguide = Slab(height/2, height/3)
set_contrast!(g, waveguide, 5)
heatmap(g.contrast, title="Design Region")


#=
Next, we add an excitation and compute the resulting field.
=#
## Add excitation
mode_input_position = 0
mode_number = 2
m = add_mode!(g, mode_input_position, mode_number)

## Compute field
solution = solve(g)
heatmap(abs.(solution))