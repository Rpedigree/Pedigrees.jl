type Pedigree{T<:Integer}
    sire::Vector{T}
    dam::Vector{T}
    perm::Vector{T}
    lappt::Vector{T} 
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
    Pedigree(sire,dam,T[],T[])
end

## Create an ordering by longest ancestral path
function laporder{T<:Integer}(sire::Vector{T},dam::Vector{T})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    anc = sizehint!(DataStructures.IntSet(),n)  # current set of ancestors
    for i in 1:n
        0 ≤ sire[i] ≤ n && 0 ≤ dam[i] ≤ n || error("row $i: sire and dam must be in 0:$n")
        sire[i] == 0 && dam[i] == 0 && push!(anc,i)  # first generation
    end

    ord = sizehint!(collect(anc),n)        # anc first in ord, also reserve space
    pop = setdiff!(DataStructures.IntSet(1:n), collect(anc))       # animals who have not yet been sorted # Collect
    lappt = sizehint!([0,length(anc)],20)  # pointer to start of each lap level
    push!(anc,0)                          # add zero to the set of ancestors already done
    nextgen = sizehint!(DataStructures.IntSet(),n)
    while length(pop) > 0
        empty!(nextgen)
        for i in pop
            sire[i] ∈ anc && dam[i] ∈ anc && push!(nextgen,i)
        end
        length(nextgen) > 0 || error("algorithm failure, empty nextgen")
        append!(ord,collect(nextgen))
        push!(lappt,lappt[end]+length(nextgen))
        union!(anc,collect(nextgen)) # Collect
        setdiff!(pop,collect(nextgen)) # Collect
    end
    length(ord) == n || error("Logic error in laporder: ord is not length n = $n")
    invp = invperm(ord)                 # checks that ord is indeed a permutation
    #ss = Array(T,(n,))
    #dd = Array(T,(n,))
    ss = zeros(T, n)
    dd = zeros(T, n)
    for i in 1:n
        j = ord[i]
        ss[i] = sire[j] > 0 ? invp[sire[j]] : zero(T)
        dd[i] = dam[j] > 0 ? invp[dam[j]] : zero(T)
    end
    ord,lappt,ss,dd
end
