# Contributing to VectorUtils.jl

Thank you for your interest in contributing to VectorUtils.jl! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Julia version and VectorUtils.jl version
- Minimal code example demonstrating the issue

### Suggesting Enhancements

Enhancement suggestions are welcome! Please open an issue with:
- A clear description of the enhancement
- Use cases and motivation
- Possible implementation approach (optional)

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature-name`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass (`julia --project=. test/runtests.jl`)
6. Update documentation as needed
7. Commit your changes (`git commit -am 'Add new feature'`)
8. Push to your branch (`git push origin feature/your-feature-name`)
9. Open a Pull Request

## Development Setup

### Clone the Repository

```bash
git clone https://github.com/mirajcs/VectorUtils.git
cd VectorUtils
```

### Install Dependencies

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Run Tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

### Build Documentation Locally

```julia
using Pkg
Pkg.activate("docs")
Pkg.instantiate()

include("docs/make.jl")
```

Then open `docs/build/index.html` in your browser.

## Coding Standards

### Style Guidelines

- Follow [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 92 characters
- Use descriptive variable names
- Add docstrings for all public functions

### Docstring Format

Use the following format for docstrings:

```julia
"""
    function_name(arg1, arg2; kwarg1=default)

Brief description of what the function does.

# Arguments
- `arg1::Type`: Description of arg1
- `arg2::Type`: Description of arg2
- `kwarg1::Type`: Description of kwarg1 (default: $default)

# Returns
- `ReturnType`: Description of return value

# Examples
```jldoctest
julia> function_name(1, 2)
3
```

# See also
- [`related_function`](@ref)
"""
```

### Testing

- Add tests for all new functionality
- Aim for high test coverage
- Use `@testset` to organize related tests
- Include edge cases and error conditions

Example test structure:

```julia
@testset "Tangent Vectors" begin
    @testset "Numeric" begin
        circle(t) = [cos(t), sin(t), 0]
        T = tangent(circle, 0.0)
        @test T ≈ [0, 1, 0] atol=1e-6
    end
    
    @testset "Symbolic" begin
        @variables t
        r = [cos(t), sin(t)]
        T = tangent_symbolic(r, t)
        @test !isnothing(T)
    end
end
```

## Documentation

### Adding Documentation

- Update relevant `.md` files in `docs/src/`
- Add docstrings to new functions
- Include examples in docstrings
- Update the API reference if adding new public functions

### Documentation Structure

```
docs/
├── make.jl                    # Documentation build script
└── src/
    ├── index.md              # Main landing page
    ├── getting_started.md    # Tutorial for new users
    ├── theory.md             # Mathematical background
    ├── api/
    │   ├── core.md          # Core API reference
    │   ├── symbolic.md      # Symbolic computation API
    │   └── numeric.md       # Numeric computation API
    └── examples/
        ├── basic.md         # Basic examples
        ├── curves.md        # Example curves
        └── frenet.md        # Frenet frame examples
```

## Types of Contributions

We welcome various types of contributions:

### Code Contributions
- New features (e.g., new curve types, algorithms)
- Performance improvements
- Bug fixes
- Code refactoring

### Documentation Contributions
- Fixing typos or unclear explanations
- Adding examples
- Improving API documentation
- Writing tutorials

### Testing Contributions
- Adding test cases
- Improving test coverage
- Performance benchmarks

### Other Contributions
- Reporting bugs
- Suggesting features
- Helping other users in issues
- Reviewing pull requests

## Code Review Process

1. All submissions require review before merging
2. Maintainers will review your PR and may request changes
3. Address review comments by pushing new commits
4. Once approved, a maintainer will merge your PR

## Questions?

If you have questions about contributing:
- Open an issue with the "question" label
- Check existing issues and discussions
- Reach out to the maintainers

## License

By contributing to VectorUtils.jl, you agree that your contributions will be licensed under the MIT License.

## Recognition

Contributors will be acknowledged in:
- The README.md file
- Release notes
- The documentation

Thank you for contributing to VectorUtils.jl! 🎉