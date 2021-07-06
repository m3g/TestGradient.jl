# TestGradient.jl

A simple function to compare the analytical implementation of a gradient with the finite-difference computation given the implemetation of the function. To be used for testing purposes, as it is very slow, but it might help finding errors in the implementation of manual gradients and sometimes on the function evaluations.

### Installation

```julia
julia> ] add https://github.com/m3g/TestGradient.jl
```

### Examples

```julia
julia> using TestGradient

julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> g!(g,x) = 2*x
g! (generic function with 1 method)

julia> x = rand(20);

julia> test_gradient(f,g!,x)

Comparing analytical and finite-difference gradients... may take a while.

Function value on test point: 8.671516638852811

Component   Analytical     Discrete        Error    Best_Step
        1  1.38356e+00  1.38356e+00  8.66637e-15  1.00000e-02
        2  1.97347e+00  1.97347e+00  4.20806e-14  1.00000e-02
        3  1.53930e-01  1.53930e-01  3.17351e-14  1.00000e-02
                            ⋮ 
       18  1.89896e-01  1.89896e-01  1.12252e-13  1.00000e-02
       19  1.37622e+00  1.37622e+00  5.84064e-14  1.00000e-02
       20  1.62172e+00  1.62172e+00  2.19071e-15  1.00000e-02

Maximum error = 1.1225244035991744e-13 at component: 18
Component   Analytical     Discrete        Error    Best_Step
       18  1.89896e-01  1.89896e-01  1.12252e-13  1.00000e-02

```

If the variables are vectors of something different than numbers, a method tries to reshape the array to to a vector of numbers. This is particularly useful for vectors of static arrays:


```julia
julia> using TestGradient 

julia> x = rand(SVector{2,Float64},10);

julia> f(x) = sum(v[1]^2 + v[2]^2 for v in x)
f (generic function with 1 method)

julia> g!(g,x) = @. g = (SVector{2,Float64}(2*v[1],2*v[2]) for v in x)
g! (generic function with 1 method)

julia> test_gradient(f,g!,x)

Comparing analytical and finite-difference gradients... may take a while.

Function value on test point: 6.5919481328739655

Computing the 1th of 20 components. Worst error: 0.0
Component   Analytical     Discrete        Error    Best_Step
        1  1.57869e+00  1.57869e+00  4.78215e-15  1.00000e-02
        2  1.04203e+00  1.04203e+00  1.91779e-14  1.00000e-02
        3  1.72248e+00  1.72248e+00  1.03128e-15  1.00000e-02
                            ⋮ 
       18  3.32460e-01  3.32460e-01  2.67153e-15  1.00000e-02
       19  1.03341e+00  1.03341e+00  8.16492e-15  1.00000e-02
       20  9.99455e-01  9.99455e-01  8.88663e-16  1.00000e-02

Maximum error = 2.2306785702197458e-12 at component: 7
Component   Analytical     Discrete        Error    Best_Step
        7  1.21440e-02  1.21440e-02  2.23068e-12  1.00000e-02

```



