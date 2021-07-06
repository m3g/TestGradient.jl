module TestGradient

using Printf

export test_gradient

struct GradientComponent
  i::Int 
  g::Float64 
  gbest::Float64 
  error::Float64 
  stepbest::Float64
end
Base.show(io::IO, c::GradientComponent) =
  @printf("%9i %12.5e %12.5e %12.5e %12.5e", c.i, c.g, c.gbest, c.error, c.stepbest)

title() = println("Component   Analytical     Discrete        Error    Best_Step")
function Base.show(io::IO, ::MIME"text/plain", c::AbstractVector{GradientComponent} )
  title()
  for i in 1:min(length(c),3)
    println(c[i])
  end
  if length(c) > 7
    @printf("%30s\n","⋮ ")
  end
  for i in max(4,length(c)-2):length(c)
    println(c[i])
  end
  max_error = findmax(map(x->x.error,c))
  print("""

  Maximum error = $(max_error[1]) at component: $(max_error[2])
  """)
  title()
  println(c[max_error[2]])
end

function discret(icomp,x,step,f)
  save = x[icomp] 
  x[icomp] = save + step 
  fplus = f(x) 
  x[icomp] = save - step 
  fminus = f(x) 
  x[icomp] = save 
  return (fplus - fminus) / (2 * step) 
end

"""

```
test_gradient(f::Function,g!::Function,x::AbstractVector{T}) where T<:Real
```

Function that performs finite difference and analytical gradient
comparision. Used only for test purpouses. The output of `test_gradient` 
is a vector of type `GradientComponent`, containing the index of the
component, the analytical and discrete gradientes, and the error associated
to the best step found for that component. 

### Example

```julia-repl
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

"""
function test_gradient(f::Function,g!::Function,x::AbstractVector{T}) where T<:Real
  n = length(x)

  println("""

  Comparing analytical and finite-difference gradients... may take a while.
  """)

  fx = f(x) 
  println("Function value on test point: ", fx)
  println()

  g = similar(x)
  g = g!(g,x)
  output = Vector{GradientComponent}(undef,0)
  t0 = time()
  eworst = 0.
  iworst = 1
  for i in 1:n
    step = 1.e-2
    stepbest = step
    error = +Inf
    gbest = 0.
    if mod(time() - t0,10) == 0
      println("Computing the $(i)th of $n components. Worst error: $eworst")
    end
    while error > 1.e-6 && step >= 1.e-20
      gcomp = discret(i,x,max(abs(g[i]*step),1e-20),f)
      if g[i] ≈ 0.  
        steperror = abs(gcomp - g[i])
      else
        steperror = abs( ( gcomp - g[i] ) / g[i] )
      end
      if steperror < error
        error = steperror
        gbest = gcomp
        stepbest = step
      end
      step = step / 10
    end
    push!(output,GradientComponent(i, g[i], gbest, error, stepbest))
    if error > eworst
      iworst = i
      eworst = error
    end
  end
  return output
end

"""

```
test_gradient(f::Function,g!::Function,x::AbstractVector{T}) where T
```

Try to reinterpret the type of input into an array of `Float64`, and then
call `test_gradient(f,g,x::AbstractVector{<:Real})`

### Example

```julia-repl
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

"""
function test_gradient(f::Function, g!::Function, x::AbstractVector{T}) where T
  test_gradient(
    x -> f(reinterpret(T,x)),
    (g,x) -> reinterpret(
      Float64,
      g!(
        reinterpret(T,g),
        reinterpret(T,x)
      ),
    ),
    reinterpret(Float64,x)
  )
end

end # module
