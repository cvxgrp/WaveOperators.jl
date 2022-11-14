@testset "Schur complement" begin
    n = 50
    h, w = 1, 2
    g = Grid(h, w, 1/n, 2π)

    waveguide = Slab(h/2, h/2)
    set_contrast!(g, waveguide, 5)

    design_region = Rectangle((h/4, w/4), (h/2, w/2))

    add_mode!(g, 0, 2)

    solution = solve(g)

    design_idx = findall(p -> is_in(design_region, p), pos.(g, CartesianIndices(g.contrast)))
    target_idx = [CartesianIndex(1, 1), CartesianIndex(2, 2)]
    int_eq = generate_G(g, design_idx, target_idx)
    G, b = int_eq.G_design, int_eq.b_design

    new_sol = (I + G)\b
    new_target = int_eq.G_target*new_sol + int_eq.b_target

    @test maximum(abs.(new_sol - solution[design_idx])) < 1e-10
    @test maximum(abs.(new_target - solution[target_idx])) < 1e-10
end