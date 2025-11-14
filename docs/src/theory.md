# Mathematical Background

This section provides the mathematical theory behind VectorUtils.jl.

## Parametric Curves

A parametric curve in 3D space is a vector-valued function:

**r**(t) = [x(t), y(t), z(t)]

where t ∈ ℝ is the parameter. The curve traces out a path in space as t varies.

### Examples

**Circle**: **r**(t) = [R cos(t), R sin(t), 0]

**Helix**: **r**(t) = [a cos(t), a sin(t), bt]

**Parabola**: **r**(t) = [t, t², 0]

## The Frenet-Serret Frame

The Frenet-Serret frame is a moving orthonormal coordinate system that follows the curve.

### Unit Tangent Vector

The unit tangent vector points in the direction of motion:

**T**(t) = **r**'(t) / ||**r**'(t)||

This is the unit vector in the direction of the velocity vector.

### Curvature

Curvature measures how rapidly the curve changes direction:

κ(t) = ||**T**'(t)|| / ||**r**'(t)|| = ||**r**'(t) × **r**''(t)|| / ||**r**'(t)||³

Alternative forms:
- κ = ||d**T**/ds|| where s is arc length
- For 2D: κ = |x'y'' - x''y'| / (x'² + y'²)^(3/2)

**Properties**:
- κ ≥ 0 (always non-negative)
- κ = 0 if and only if the curve is a straight line
- κ = 1/R for a circle of radius R

### Principal Normal Vector

The principal normal vector points toward the center of curvature:

**N**(t) = **T**'(t) / ||**T**'(t)||

When κ ≠ 0, this can also be written as:

**N**(t) = (**r**''(t) - (**r**''(t) · **T**(t))**T**(t)) / ||**r**''(t) - (**r**''(t) · **T**(t))**T**(t)||

### Binormal Vector

The binormal vector is perpendicular to both **T** and **N**:

**B**(t) = **T**(t) × **N**(t)

It is normal to the osculating plane (the plane containing **T** and **N**).

### Torsion

Torsion measures how much the curve twists out of the osculating plane:

τ(t) = -(**N**'(t) · **B**(t)) = (**r**'(t) × **r**''(t)) · **r**'''(t) / ||**r**'(t) × **r**''(t)||²

**Properties**:
- τ can be positive, negative, or zero
- τ = 0 if and only if the curve is planar
- τ > 0: right-handed twist
- τ < 0: left-handed twist

## The Frenet-Serret Formulas

The frame vectors satisfy the following differential equations:

d**T**/ds = κ**N**

d**N**/ds = -κ**T** + τ**B**

d**B**/ds = -τ**N**

where s is arc length. These are the Frenet-Serret formulas.

## Geometric Interpretations

### Curvature

- **Physical meaning**: Rate of change of direction per unit arc length
- **Geometric meaning**: Reciprocal of the radius of the osculating circle
- **Osculating circle**: The circle that best approximates the curve at a point
  - Center: **r**(t) + (1/κ)**N**(t)
  - Radius: R = 1/κ

### Torsion

- **Physical meaning**: Rate at which the osculating plane rotates
- **Geometric meaning**: Measure of how non-planar the curve is
- **Sign convention**: 
  - Positive torsion: right-hand screw motion
  - Negative torsion: left-hand screw motion

## Special Cases

### Straight Lines

For a straight line **r**(t) = **a** + t**b**:
- κ = 0
- τ = undefined (or 0)
- **T** = constant
- **N**, **B** undefined

### Circles

For a circle of radius R:
- κ = 1/R (constant)
- τ = 0 (planar)
- **T** rotates in the plane
- **N** points toward center
- **B** is constant (perpendicular to plane)

### Helices

For a circular helix **r**(t) = [a cos(t), a sin(t), bt]:
- κ = a/(a² + b²) (constant)
- τ = b/(a² + b²) (constant)
- All three vectors rotate as you move along the curve

## Arc Length

The arc length from t = a to t = b is:

s(b) - s(a) = ∫ₐᵇ ||**r**'(t)|| dt

The arc length parameterization satisfies:

ds/dt = ||**r**'(t)||

When a curve is parameterized by arc length, ||**r**'(s)|| = 1.

## Osculating, Normal, and Rectifying Planes

At each point on the curve, three important planes exist:

**Osculating plane**: Contains **T** and **N**, perpendicular to **B**
- The plane in which the curve is "most nearly planar"
- Normal vector: **B**

**Normal plane**: Contains **N** and **B**, perpendicular to **T**
- Perpendicular to the direction of motion
- Normal vector: **T**

**Rectifying plane**: Contains **T** and **B**, perpendicular to **N**
- Normal vector: **N**

## Computation Methods

### Numeric Differentiation

VectorUtils uses finite differences for numeric derivatives:

**r**'(t) ≈ (**r**(t + h) - **r**(t - h)) / (2h)

**r**''(t) ≈ (**r**(t + h) - 2**r**(t) + **r**(t - h)) / h²

where h is a small step size (typically 10⁻⁵).

### Symbolic Differentiation

For symbolic computations, VectorUtils uses Symbolics.jl to compute exact derivatives:

d/dt [x(t), y(t), z(t)] = [dx/dt, dy/dt, dz/dt]

This provides exact analytical expressions.

## Applications

### Computer Graphics

- Camera paths for smooth animations
- Curve interpolation
- Surface generation by sweeping frames

### Robotics

- Path planning for smooth trajectories
- End-effector orientation along paths
- Twist and bending analysis

### Differential Geometry

- Study of space curves
- Surface theory (curves on surfaces)
- Minimal surfaces

### Physics

- Particle trajectories
- Electromagnetic field lines
- String theory

## References

1. Do Carmo, M. P. (1976). *Differential Geometry of Curves and Surfaces*
2. Pressley, A. (2010). *Elementary Differential Geometry*
3. Struik, D. J. (1961). *Lectures on Classical Differential Geometry*
4. O'Neill, B. (2006). *Elementary Differential Geometry*

## Further Reading

For deeper understanding:
- Differential geometry textbooks
- Vector calculus resources
- Computer graphics literature on curve design
- Robotics texts on trajectory planning