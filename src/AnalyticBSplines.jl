module AnalyticBSplines

using Polynomials
using BSplinesQuasi
import BSplinesQuasi: AbstractKnotSet

include("intervals.jl")
include("ppoly.jl")

function gen_basis(t::AbstractKnotSet{k,ml,mr,T}, bl=one(T), br=one(T)) where {k,ml,mr,T}
    d(a,b) = b == 0 ? 0 : a/b
    B = [[PPoly([((t[j],t[j+1]) => [one(T)])])
          for j = 1:length(t)-1]]
    for k = 2:order(t)
        Bₖ = [d(Poly([-t[i],1]),(t[i+k-1]-t[i]))*B[end][i] +
              d(Poly([t[i+k],-1]),(t[i+k]-t[i+1]))*B[end][i+1]
              for i = 1:length(t)-k]
        push!(B,Bₖ)
    end
    B[end][1] *= bl
    B[end][end] *= br
    B
end

function scalar_op(Bₖ, t::AbstractKnotSet{k,ml,mr,T}, f::Poly=Poly([one(T)])) where {k,ml,mr,T}
    B = Bₖ[end]
    n = length(B)
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

function derop(L, R, t::AbstractKnotSet{k,ml,mr,T}, o::Integer) where {k,ml,mr,T}
    m = length(L[end])
    n = length(R[end])
    O = Matrix{T}(undef, m,n)
    l(v::Rational) = denominator(v) == 1 ? numerator(v) : 1.0*v
    l(v) = v
    d,r = divrem(o,2)
    a,b = d,d+r
    for i = 1:m
        for j = 1:n
            O[i,j] = (iseven(a) ? 1 : -1)*polyint(polyder(L[end][i],a)*polyder(R[end][j],b), l(first(t)),l(last(t)))
        end
    end
    O
end

derop(B, t, o) = derop(B, B, t, o)

export gen_basis, scalar_op, derop

end # module
