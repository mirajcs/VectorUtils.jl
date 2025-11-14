# Getting Started

This guide will help you get started with VectorUtils.jl for analyzing parametric curves.

## Installation

Install VectorUtils.jl using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/mirajcs/VectorUtils")
```

Then load the package:

```julia
using VectorUtils
```

## Basic Concepts

### Parametric Curves

A parametric curve in 3D space is defined as a vector-valued function:

**r**(t) = [x(t), y(t), z(t)]

where t is the parameter (often representing time or arc length).

### The Frenet-Serret Frame

For any point on a curve, the Frenet-Serret frame consists of three orthonormal vectors:

- **T** (Tangent): Points in the direction of motion
- **N** (Normal): Points toward the center of curvature
- **B** (Binormal): Perpendicular to both T and N

These vectors form a right-handed orthonormal basis at each point.

## Your First Computation

Let's compute properties of a simple circle:

```julia
using VectorUtils

# Define a circle in the xy-plane with radius 1
circle(t) = [cos(t), sin(t), 0]

# Evaluate at t = π/4
t₀ = π/4

# Compute tangent vector
T = tangent(circle, t₀)
println("Tangent: ", T)

# Compute curvature
κ = curvature(circle, t₀)
println("Curvature: κ = ", κ)  # Should be 1.0

# Compute full Frenet frame
frame = frenet_frame(circle, t₀)
println("Normal: ", frame.normal)
println("Binormal: ", frame.binormal)
```

## Working with Different Curves

### 2D Curves

For planar curves, set z = 0:

```julia
# Parabola
parabola(t) = [t, t^2, 0]

κ = curvature(parabola, 1.0)
```

### 3D Space Curves

The most common example is a helix:

```julia
# Circular helix with radius a and pitch 2πb
a, b = 2.0, 0.5
helix(t) = [a*cos(t), a*sin(t), b*t]

# At one complete turn
t₁ = 2π
frame = frenet_frame(helix, t₁)
κ = curvature(helix, t₁)
τ = torsion(helix, t₁)

println("Curvature: ", κ)
println("Torsion: ", τ)
```

## Numeric vs Symbolic

### Numeric Computation

Use numeric functions for specific parameter values:

```julia
# Define curve
r(t) = [cos(t), sin(t), t]

# Evaluate at specific point
T = tangent(r, π/2)
κ = curvature(r, π/2)
```

### Symbolic Computation

Use symbolic functions for analytical expressions:

```julia
using Symbolics

# Define symbolic parameter
@variables t

# Define curve symbolically
r_sym = [cos(t), sin(t), t]

# Get symbolic expressions
T_sym = tangent_symbolic(r_sym, t)
κ_sym = curvature_symbolic(r_sym, t)

# Simplify
κ_simplified = simplify(κ_sym)
```

## Plotting Your Results

While VectorUtils focuses on computation, you can easily visualize results:

```julia
using Plots

# Generate curve points
ts = range(0, 2π, length=100)
curve(t) = [cos(t), sin(t), t]
points = [curve(t) for t in ts]

xs = [p[1] for p in points]
ys = [p[2] for p in points]
zs = [p[3] for p in points]

# Plot the curve
plot3d(xs, ys, zs, label="Helix", linewidth=2)

# Add Frenet frame at a point
t₀ = π
frame = frenet_frame(curve, t₀)
p₀ = curve(t₀)

# Plot frame vectors
quiver!([p₀[1]], [p₀[2]], [p₀[3]], 
        quiver=([frame.tangent[1]], [frame.tangent[2]], [frame.tangent[3]]),
        color=:red, label="T")
quiver!([p₀[1]], [p₀[2]], [p₀[3]], 
        quiver=([frame.normal[1]], [frame.normal[2]], [frame.normal[3]]),
        color=:green, label="N")
quiver!([p₀[1]], [p₀[2]], [p₀[3]], 
        quiver=([frame.binormal[1]], [frame.binormal[2]], [frame.binormal[3]]),
        color=:blue, label="B")
```

## Common Workflows

### Analyzing a New Curve

1. Define your parametric curve as a function
2. Choose parameter values of interest
3. Compute desired properties (tangent, curvature, torsion)
4. Verify results (check orthonormality, signs, etc.)

### Finding Special Points

```julia
# Find where curvature is maximum
using Optim

curve(t) = [t, t^2, t^3]
κ_func(t) = curvature(curve, t)

# Maximize curvature in interval [0, 2]
result = optimize(t -> -κ_func(t), 0.0, 2.0)
t_max = result.minimizer
κ_max = κ_func(t_max)

println("Maximum curvature κ = $κ_max at t = $t_max")
```

## Next Steps

- Explore the [API Reference](api/core.md) for all available functions
- Check out detailed [Examples](examples/basic.md)
- Read the [Mathematical Background](theory.md) for theory
- Try the example curves in the next sections!