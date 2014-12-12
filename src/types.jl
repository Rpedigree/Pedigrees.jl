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

## Create an ordering by longest ancestral path
function laporder{T<:Integer}(sire::Vector{T},dam::Vector{T})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    for i in 1:n
        0 ≤ sire[i] ≤ n && 0 ≤ dam[i] ≤ n || error("row $i: sire and dam must be in 0:$n")
    end

    anc = IntSet(find(sire+dam .== 0))    # animals without ancestors in the pedigree
    ord = sizehint(collect(anc),n)        # anc first in ord, also reserve space
    pop = setdiff!(IntSet(1:n),anc)       # animals who have not yet been sorted
    lappt = sizehint([0,length(anc)],20)  # pointer to start of each lap level
    push!(anc,0)                          # add zero to the set of ancestors already done
    nextgen = sizehint(IntSet(),n)
    while length(pop) > 0
        empty!(nextgen)
        for i in pop
            sire[i] ∈ anc && dam[i] ∈ anc && push!(nextgen,i)
        end
        append!(ord,collect(nextgen))
        push!(lappt,lappt[end]+length(nextgen))
        union!(anc,nextgen)
        setdiff!(pop,nextgen)
    end
    length(ord) == n || error("Logic error in laporder: ord is not length n = $n")
    invp = invperm(ord)                 # checks that ord is indeed a permutation
    ss = Array(T,(n,))
    dd = Array(T,(n,))
    for i in 1:n
        j = ord[i]
        ss[i] = sire[j] > 0 ? invp[sire[j]] : zero(T)
        dd[i] = dam[j] > 0 ? invp[dam[j]] : zero(T)
    end
    ord,lappt,ss,dd
end
