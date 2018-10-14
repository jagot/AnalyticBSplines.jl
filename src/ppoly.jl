import Base: +, -, *, /, ==
import Polynomials: polyder, polyint

struct PPoly{T}
    pis::Vector{Pair{Interval{T},Poly{T}}}
    function PPoly(pis::Vector{Pair{Interval{T},Poly{T}}}) where T
        for (j,i) in enumerate(first.(pis))
            for i′ in first.(pis[j+1:end])
                g = i ∩ i′
                !isempty(g) && error("Overlapping intervals: $(i) ∩ $(i′) = $(g) ≠ ∅")
            end
        end
        sort!(pis, by = ip -> ip[1].a)
        new{T}(pis)
    end
end

PPoly() = PPoly{Int64}(Vector{Pair{Interval{Int64},Poly{Int64}}}())

# Convenience constructor
PPoly(pis::Vector{Pair{Tuple{T,T},Vector{T}}}) where T =
    PPoly([Interval(i...)=>Poly(p)
           for (i,p) in pis])

# Constructor promoting all arguments to a common type
function PPoly(pis::Vector)
    isempty(pis) && PPoly()
    T = typeof(pis[1][1][1])
    for (i,p) in pis
        T = promote_type(T, typeof(i[1]), typeof(i[2]), eltype(p))
    end
    PPoly([Interval(T.(i)...)=>Poly(T.(p))
           for (i,p) in pis])
end

function show(io::IO, pp::PPoly{T}) where T
    write(io, "Piecewise polynomial{$(T)}")
    for (i,p) in pp.pis
        write(io, "\n $(p) on $(i)")
    end
end

# Evaluation
function (pp::PPoly)(x::Number)
    S = 0
    for (i,p) ∈ pp.pis
        if x ∈ i
            S += polyval(p, x)
        end
    end
    S
end

# Arithmetic
function +(p₁::PPoly{T}, p₂::PPoly{T}) where T
    pis = Vector{Pair{Interval{T},Poly{T}}}()
    intersections = Matrix{Bool}(undef, length(p₁.pis), length(p₂.pis))
    isub = Vector{Vector{Interval{T}}}()
    i′sub = Vector{Vector{Interval{T}}}()
    for (k,(i,p)) ∈ enumerate(p₁.pis)
        push!(isub, [i])
        for (l,(i′,p′)) ∈ enumerate(p₂.pis)
            ii′ = i ∩ i′
            intersections[k,l] = !isempty(ii′)
        end
    end
    for (i′,p′) ∈ p₂.pis
        push!(i′sub, [i′])
    end

    irows = reduce(|, intersections, dims=2)
    icols = reduce(|, intersections, dims=1)
    for (k,(i,p)) ∈ zip(eachindex(irows), p₁.pis)
        !irows[k] && continue
        for (l,(i′,p′)) ∈ zip(eachindex(icols), p₂.pis)
            !icols[l] && continue
            ii′ = i ∩ i′
            push!(pis, ii′ => p+p′)
            isub[k] = collapse!(Vector{Interval{T}}(vcat([(is \ ii′)[1] for is ∈ isub[k]]...)))
            i′sub[l] = collapse!(Vector{Interval{T}}(vcat([(is \ ii′)[1] for is ∈ i′sub[l]]...)))
        end
    end
    for (is,p) ∈ zip(isub,last.(p₁.pis))
        for i ∈ is
            push!(pis, i => p)
        end
    end
    for (i′s,p′) ∈ zip(i′sub,last.(p₂.pis))
        for i′ ∈ i′s
            push!(pis, i′ => p′)
        end
    end
    filter!(ip -> !isempty(ip[1]) && ip[2] != zero(T), pis)
    PPoly(pis)
end

-(p::PPoly) = PPoly([i=>-p for (i,p) in p.pis])
-(p₁::PPoly{T}, p₂::PPoly{T}) where T = p₁ + (-p₂)

==(p₁::PPoly{T}, p₂::PPoly{T}) where T = isempty((p₁-p₂).pis)

*(pp::PPoly, n::Number) = PPoly([i=>p*n for (i,p) in pp.pis])
*(n::Number, pp::PPoly) = pp*n
/(pp::PPoly, n::Number) = PPoly([i=>p/n for (i,p) in pp.pis])

*(pp::PPoly, p′::Poly) = PPoly([i=>p*p′ for (i,p) in pp.pis])
*(p′::Poly,pp::PPoly) = pp*p′

function *(p₁::PPoly{T}, p₂::PPoly{T}) where T
    pis = Vector{Pair{Interval{T},Poly{T}}}()
    intersections = Matrix{Bool}(undef, length(p₁.pis), length(p₂.pis))
    isub = Vector{Vector{Interval{T}}}()
    i′sub = Vector{Vector{Interval{T}}}()
    for (k,(i,p)) ∈ enumerate(p₁.pis)
        push!(isub, [i])
        for (l,(i′,p′)) ∈ enumerate(p₂.pis)
            ii′ = i ∩ i′
            intersections[k,l] = !isempty(ii′)
        end
    end
    for (i′,p′) ∈ p₂.pis
        push!(i′sub, [i′])
    end

    irows = reduce(|, intersections, dims=2)
    icols = reduce(|, intersections, dims=1)
    for (k,(i,p)) ∈ zip(eachindex(irows), p₁.pis)
        !irows[k] && continue
        for (l,(i′,p′)) ∈ zip(eachindex(icols), p₂.pis)
            !icols[l] && continue
            ii′ = i ∩ i′
            push!(pis, ii′ => p*p′)
        end
    end
    filter!(ip -> !isempty(ip[1]) && ip[2] != zero(T), pis)
    PPoly(pis)
end

# /(pp::PPoly, p′::Poly) = PPoly([i=>p/p′ for (i,p) in pp.pis])

function polyder(pp::PPoly, k=1)
    pis = [i => polyder(p,k)
           for (i,p) in pp.pis]
    filter!(ip -> ip[2]!=0, pis)
    PPoly(pis)
end

polyint(pp::PPoly) = PPoly([i=>polyint(p) for (i,p) in pp.pis])

function polyint(pp::PPoly{T}, a::Number, b::Number) where T
    ab = Interval(a, b)
    S = 0
    for (i,p) ∈ pp.pis
        if !(isempty(i ∩ ab))
            S += polyint(p, max(i.a,a), min(i.b,b))
        end
    end
    S
end

export PPoly, polyder, polyint
