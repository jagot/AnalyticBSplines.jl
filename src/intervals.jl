struct Interval{T}
    a::T
    b::T
end
Interval(a::T,b::U) where {T,U} = Interval(promote(a,b)...)

import Base: ∪, ∩, \, ⊆, ∈, isempty, show, isless

isless(i₁::Interval,i₂::Interval) = i₁.a < i₂.a
adjacent(i₁::Interval,i₂::Interval) = i₁.b == i₂.a
isempty(i::Interval) = i.a >= i.b

⊂(i₁::Interval,i₂::Interval) = i₂.a < i₁.a && i₁.b < i₂.b
⊆(i₁::Interval,i₂::Interval) = i₂.a <= i₁.a && i₁.b <= i₂.b

# This imposes right-continuity, as de Boor uses.
∈(e::Number, i::Interval) = i.a < e && e ≤ i.b

function ∩(i₁::Interval,i₂::Interval)
    a = max(i₁.a,i₂.a)
    b = min(i₁.b,i₂.b)
    Interval(a,b)
end

function ∪(i₁::Interval,i₂::Interval)
    a = min(i₁.a,i₂.a)
    b = max(i₁.b,i₂.b)
    if !isempty(i₁ ∩ i₂) || adjacent(i₁,i₂) || adjacent(i₂,i₁)
        Interval(a,b)
    else
        Interval(1,0)
    end
end

function ∪(i₁::Vector{Interval{T}}, i₂::Interval{T}) where T
    sort!(vcat([i ∪ i₂ for i ∈ i₁]...))
    # Should merge possibly adjacent intervals
end

function \(i₁::Interval,i₂::Interval)
    ii = i₁ ∩ i₂
    ia = Interval(i₁.a,i₂.a)
    ib = Interval(i₂.b,i₁.b)
    if isempty(ii)
        [i₁],[i₂]
    elseif i₂ ⊆ i₁
        if isempty(ia) && isempty(ib)
            Interval[],Interval[]
        elseif isempty(ia)
            [ib],Interval[]
        elseif isempty(ib)
            [ia],Interval[]
        else
            [ia,ib],Interval[]
        end
    elseif i₁ ⊂ i₂
        Interval[],[Interval(i₂.a,i₁.a),Interval(i₁.b,i₂.b)]
    else # i₁ and i₂ partially overlap
        if i₁.a <= i₂.a
            isempty(ia) ? Interval[] : [ia],[Interval(i₁.b,i₂.b)]
        else
            isempty(ib) ? Interval[] : [ib],[Interval(i₂.a,i₁.a)]
        end
    end
end

function \(i₁::Vector{Interval}, i₂::Interval)
    is = Interval[]
    for i ∈ i₁
        i′ = i \ i₂
        !isempty(i′) && append!(is, i′)
    end
    is
end

function collapse!(is::Vector{Interval{T}}) where T
    isempty(is) && return is
    sort!(is)
    nis = Interval{T}[is[1]]
    i = 1
    for j = 2:length(is)
        ij = nis[i] ∪ is[j]
        if !isempty(ij)
            nis[i] = ij
        else
            push!(nis, is[j])
            i += 1
        end
    end
    is[1:length(nis)] = nis
    resize!(is, length(nis))
    is
end

function show(io::IO, i::Interval)
    if isempty(i)
        write(io, "Empty Interval")
    else
        write(io, "Interval[$(i.a),$(i.b)]")
    end
end

export Interval, ⊂
