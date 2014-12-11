type Pedigree{T<:Integer}
    sire::Vector{T}
    dam::Vector{T}
    lap::Vector{T}
end

function Pedigree{T<:Integer}(sire::Vector{T},dam::Vector{T})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    for i in 1:n
        zero(T) ≤ sire[i] ≤ n && zero(T) ≤ dam[i] ≤ n ||
        error("row $i: sire and dam must be in 0:$n")
        sire[i] == i && error("Malformed pedigree, sire[$i] == $i")
        dam[i] == i && error("Malformed pedigree, dam[$i] == $i")
        sire[i] < i && dam[i] < i || error("order failed at $i, with sire and dam, $(sire[i]), $(dam[i])")
    end
    Pedigree(sire,dam,lap(sire,dam))
end

# evaluate the longest ancestral path for each animal in an ordered pedigree
function lap{T<:Integer}(sire::Vector{T},dam::Vector{T}) 
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    LAP = zeros(T,n)
    mone = convert(T,-1)
    @inbounds for i in 1:n
        LAP[i] = one(T)+max(sire[i]>0 ? LAP[sire[i]] : mone, dam[i]>0 ? LAP[dam[i]] : mone)
    end
    LAP
end

function orderped{T<:Integer}(sire::Vector{T},dam::Vector{T})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    for i in 1:n
        0 ≤ sire[i] ≤ n && 0 ≤ dam[i] ≤ n || error("row $i: sire and dam must be in 0:$n")
    end
    ord = sizehint(T[],n)               # empty array with space reserved
    ## animals with both sire and dam unknown are put at the end of ord before reversing
    unknownrents = IntSet((1:n)[sire .== 0 & dam .== 0]) 
    pop = setdiff(IntSet(1:n),unknownrents)
    while length(pop) > 0
        inds = collect(pop)
        parents = union(IntSet(sire[inds]),IntSet(dam[inds]))
        append!(ord,reverse!(collect(setdiff(pop,parents))))
        intersect!(pop,parents)
    end
    reverse!(append!(ord,reverse!(collect(unknownrents))))
    length(ord) == n || error("Logic error in orderped, length(ord) == $(length(ord)) != $n")

    ip = invperm(ord)                   # will check that ord is a permutation
    ss = permute!([sire[i] == zero(T) ? zero(T) : ip[sire[i]] for i in 1:n],ord)
    dd = permute!([dam[i] == zero(T) ? zero(T) : ip[dam[i]] for i in 1:n],ord)
    ord,ss,dd
end
