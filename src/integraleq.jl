@doc raw"""
    IntegralEquation

This object represents the physics equation
```math
z + G\mathbf{diag}(\theta)z = Gb,
```
in compact form. Specifically, `G_design` maps from the design region to the
design region, and `G_target` maps from the design region to the target. By
only dealing with the reduced form of the system, we avoid unnecessary
computation. The reduced form of the physics is
```math
\begin{aligned}
z^\mathrm{design} + G^\mathrm{design}z^\mathrm{design} &= b^\mathrm{design}\\
z^\mathrm{target} + G^\mathrm{target}z^\mathrm{design} &= b^\mathrm{target}
\end{aligned}
```
"""
struct IntegralEquation
    G_design::Matrix{ComplexF64}
    G_target::Matrix{ComplexF64}
    b_design::Vector{ComplexF64}
    b_target::Vector{ComplexF64}
end