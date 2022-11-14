export Rectangle, Circle, Slab, VerticalLine, HorizontalLine, ShapeSlice, is_in 

abstract type Shape end

Base.broadcastable(x::S) where S <: Shape = Ref(x)

abstract type Item end

Position = Tuple{Float64, Float64}
Dimensions = Tuple{Float64, Float64}

"""
    is_in(s::Shape, p::Position)

Checks if `p` (a tuple `(x::Float64, y::Float64)`) is in the shape `s`.
"""
function is_in end

struct Rectangle <: Shape
    position::Position
    shape::Dimensions
end
"""
    Rectangle(w::N, h::M) where {N <: Real, M <: Real}

Creates a `Rectangle` at (0,0) with dimensions `w` by `h`.
"""
Rectangle(w::N, h::M) where {N <: Real, M <: Real} = Rectangle((0.0, 0.0), (float(w), float(h)))

"""
    Rectangle(position::T, shape::U) where {T <: Tuple, U <: Tuple}

Creates a `Rectangle` at `position`` with dimensions `shape`. Both inputs
are tuples of floats (`Tuple{Float64, Float64}`).
"""
Rectangle(position::T, shape::U) where {T <: Tuple, U <: Tuple} = Rectangle(float.(position), float.(shape))

is_in(r::Rectangle, p::Position) = all(@. 0 <= (p - r.position) <= r.shape)

struct Circle <: Shape
    position::Position
    radius::Float64
end
"""
    Circle(position::T, radius) where T <: Tuple

Creates a `Circle` at `position` with radius `radius`.
"""
Circle(position::T, radius) where T <: Tuple = Circle(float.(position), float(radius))

is_in(c::Circle, p::Position) = norm(c.position .- p) <= r.radius

# Maybe we should add possible orientations?
"""
    Slab(position::Float64, width::Float64)

Creates a `Slab` (horizontal rectangle that strectes across the entire domain)
at y position `position` with (vertical) width `width`.
"""
struct Slab <: Shape
    position::Float64
    width::Float64
end

is_in(s::Slab, p::Position) = -s.width/2 <= p[1] - s.position <= s.width/2

abstract type Line <: Shape end

"""
    VerticalLine(x::Float64)

Creates a vertical line (across the entire domain) at x position `x`.
"""
struct VerticalLine <: Line
    x::Float64
end

"""
    HoriztonalLine(y2::Float64)

Creates a horiztonal line (across the entire domain) at y position `y`.
"""
struct HorizontalLine <: Line
    y::Float64
end


"""
    ShapeSlice(shape::Shape, line::VerticalLine)

Creates a vertical slice (defined by `line`) of `shape`.
"""
struct ShapeSlice{S} <: Shape where {S <: Shape}
    shape::S
    line::VerticalLine
end