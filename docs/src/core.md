# Core Functions API

This page documents the main functions provided by VectorUtils.jl.

## Tangent Vectors

```@docs
tangent
tangent_symbolic
unit_tangent
```

### Description

The tangent vector **T**(t) points in the direction of increasing parameter along the curve.

**Formula**: **T**(t) = **r**'(t) / ||**r**'(t)||

### Example

```julia
# Circle tangent
circle(t) = [cos(t), sin(t), 0]
T = tangent(circle, π/4)
# T ≈ [-0.707, 0.707, 0]
```

---

## Frenet Frame

```@docs
frenet_frame
frenet_frame_symbolic
FrenetFrame
```

### Description

The Frenet-Serret frame provides a moving coordinate system along the curve:

- **T**: Unit tangent vector
- **N**: Principal normal vector (points toward center of curvature)
- **B**: Binormal vector (perpendicular to osculating plane)

These vectors satisfy:
- **T** · **N** = 0
- **T** · **B** = 0
- **N** · **B** = 0
- **B** = **T** × **N**

### Example

```julia
helix(t) = [cos(t), sin(t), t]
frame = frenet_frame(helix, π/2)

println("T: ", frame.tangent)
println("N: ", frame.normal)
println("B: ", frame.binormal)

# Verify orthonormality
@assert abs(dot(frame.tangent, frame.normal)) < 1e-10
@assert abs(norm(frame.tangent) - 1.0) < 1e-10
```

---

## Curvature

```@docs
curvature
curvature_symbolic
radius_of_curvature
```

### Description

Curvature κ(t) measures how quickly the curve changes direction. It is the magnitude of the rate of change of the unit tangent vector with respect to arc length.

**Formula**: κ(t) = ||**r**'(t) × **r**''(t)|| / ||**r**'(t)||³

For a circle of radius R, κ = 1/R (constant).

### Properties

- κ ≥ 0 always (non-negative)
- κ = 0 for straight lines
- Larger κ means tighter curve
- Units: 1/length

### Example

```julia
# Circle of radius 2
circle(t) = [2*cos(t), 2*sin(t), 0]
κ = curvature(circle, 0.0)
# κ = 0.5 = 1/R

# Radius of curvature
R = radius_of_curvature(circle, 0.0)
# R = 2.0
```

---

## Torsion

```@docs
torsion
torsion_symbolic
```

### Description

Torsion τ(t) measures how much the curve twists out of its osculating plane. It is the rate at which the binormal vector rotates.

**Formula**: τ(t) = (**r**'(t) × **r**''(t)) · **r**'''(t) / ||**r**'(t) × **r**''(t)||²

### Properties

- τ can be positive, negative, or zero
- τ = 0 for planar curves
- τ > 0: right-handed twist
- τ < 0: left-handed twist
- Units: 1/length

### Example

```julia
# Planar curve (circle)
circle(t) = [cos(t), sin(t), 0]
τ = torsion(circle, 0.0)
# τ = 0 (no twisting)

# Helix
helix(t) = [cos(t), sin(t), t]
τ = torsion(helix, 0.0)
# τ ≈ 0.5 (constant twist)
```

---

## Derivative Functions

```@docs
derivative
derivative_numeric
second_derivative
third_derivative
```

### Description

These functions compute numerical derivatives of vector-valued functions using finite differences.

### Example

```julia
r(t) = [t^2, t^3, t^4]

r_prime = derivative(r, 1.0)
r_double_prime = second_derivative(r, 1.0)
r_triple_prime = third_derivative(r, 1.0)
```

---

## Utility Functions

```@docs
normalize_vector
cross_product
arc_length
arc_length_numeric
```

### normalize_vector

Returns a unit vector in the same direction.

```julia
v = [3.0, 4.0, 0.0]
v_unit = normalize_vector(v)
# v_unit = [0.6, 0.8, 0.0]
# norm(v_unit) = 1.0
```

### cross_product

Computes the cross product of two 3D vectors.

```julia
a = [1.0, 0.0, 0.0]
b = [0.0, 1.0, 0.0]
c = cross_product(a, b)
# c = [0.0, 0.0, 1.0]
```

### arc_length

Computes the arc length of a curve between two parameter values.

```julia
circle(t) = [cos(t), sin(t), 0]

# Arc length from 0 to π
s = arc_length(circle, 0.0, π)
# s = π (half circumference of unit circle)
```

---

## Type Definitions

```@docs
FrenetFrame
CurveProperties
```

### FrenetFrame

A structure containing the Frenet-Serret frame vectors.

**Fields**:
- `tangent::Vector`: Unit tangent vector T
- `normal::Vector`: Principal normal vector N
- `binormal::Vector`: Binormal vector B

### CurveProperties

A structure containing all geometric properties at a point.

**Fields**:
- `position::Vector`: Position vector r(t)
- `tangent::Vector`: Unit tangent T
- `normal::Vector`: Principal normal N
- `binormal::Vector`: Binormal B
- `curvature::Real`: Curvature κ
- `torsion::Real`: Torsion τ

### Example

```julia
curve(t) = [cos(t), sin(t), t]
props = curve_properties(curve, π/4)

println("Position: ", props.position)
println("Curvature: ", props.curvature)
println("Torsion: ", props.torsion)
```