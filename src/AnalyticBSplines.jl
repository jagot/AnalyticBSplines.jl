module AnalyticBSplines

using Polynomials
using BSplines
import BSplines: AbstractKnotSet

include("intervals.jl")
include("ppoly.jl")

function gen_basis(t::AbstractKnotSet{T}) where T
    d(a,b) = b == 0 ? 0 : a/b
    B = [[PPoly([((t[j],t[j+1]) => [one(T)])])
          for j = 1:length(t)-1]]
    for k = 2:order(t)
        Bₖ = [d(Poly([-t[i],1]),(t[i+k-1]-t[i]))*B[end][i] +
              d(Poly([t[i+k],-1]),(t[i+k]-t[i+1]))*B[end][i+1]
              for i = 1:length(t)-k]
        push!(B,Bₖ)
    end
    B
end

function scalar_op(Bₖ, t::AbstractKnotSet{T}, f::Poly=Poly([one(T)])) where T
    B = Bₖ[end]
    n = order(t)+1
    O = Matrix{T}(undef, n,n)
    l(v::Rational) = denominator(v) == 1 ? numerator(v) : 1.0*v
    l(v) = v
    for i = 1:n
        for j = 1:n
            O[i,j] = polyint(B[i]*f*B[j], l(first(t)),l(last(t)))
        end
    end
    O
end

function derop(Bₖ, t, o)
    n = order(t)+1
    O = Matrix{Any}(undef, n,n)
    l(v::Rational) = denominator(v) == 1 ? numerator(v) : 1.0*v
    l(v) = v
    for i = 1:n
        for j = 1:n
            O[i,j] = polyint(Bₖ[end][i]*polyder(Bₖ[end][j],o), l(first(t)),l(last(t)))
        end
    end
    O
end

end # module
