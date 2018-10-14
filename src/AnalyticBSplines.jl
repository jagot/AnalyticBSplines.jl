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

end # module
