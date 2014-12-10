module pedigree
if isfile(joinpath(Pkg.dir("pedigree"),"deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("pedigree not properly installed. Please run Pkg.build(\"pedigree\")")
end

    export Pedigree, inbreeding
           
    type Pedigree
        sire::Vector{Int}
        dam::Vector{Int}
        lap::Vector{Int} # longest ancestral path, lap[1]=-1 for unknown parent
        function Pedigree(sire::Vector{Int},dam::Vector{Int})
            (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
            ordered = true
            for i in 1:n
                0 <= sire[i] <= n && 0 <= dam[i] <= n || error("sire and dam must be in 0:n")
                sire[i] == i && error("Malformed pedigree, sire[$i] == $i")
                dam[i] == i && error("Malformed pedigree, dam[$i] == $i")
                sire[i] < i && dam[i] < i || (ordered = false)
            end
        if !ordered
            reorder_pedigree(sire,dam)
        end
        new(sire,dam,lap(sire,dam))
        end
    end

    function lap(sire::Vector{Int},dam::Vector{Int}) # longest ancestral path
        (n = length(sire)) == length(dam) || throw(DimensionMismatch())
        LAP = zeros(Int,n+1)
        LAP[1] = -1
        for i in 1:n
            LAP[i+1] = max(LAP[sire[i]+1],LAP[dam[i]+1]) + 1
        end
        LAP
    end

    function inbreeding(p::Pedigree)
        sire = p.sire
        dam = p.dam
        lap = p.lap
        n = length(sire)
        ## Arrays of length n+1 have a placeholder for missing parents in the first element
        F = zeros(n+1)          # inbreeding coefficients
        ccall((:inbreeding,libpedigree),Void,(Ptr{Int},Ptr{Int},Ptr{Int},Int,Ptr{Cdouble}),
              sire,dam,lap,n,F)
        return F

        F[1] = -1.              # initialize F for unknown parents
        L = zeros(n+1)
        B = zeros(n)            # diagonal factor in A
        Anc = zeros(Int,n+1)    # Ancestor
        maxlap = -1
        for i in 1:n            # Evaluate the longest ancestral path
            S = sire[i]
            D = dam[i]
            (LAP[i+1] = max(LAP[S+1],LAP[D+1])+1) > maxlap && (maxlap = LAP[i+1])
        end
        SI = zeros(Int,maxlap)              # start indices
        MI = zeros(Int,maxlap)              # minor indices
        for i in 1:n
            S = sire[i]
        D = dam[i]
        B[i] = 0.5 - 0.25*(F[S+1]+F[D+1])
        for j in 1:LAP[i+1]             # adjust start and minor
            SI[j] += 1
            MI[j] += 1
        end
        if S == 0 && D == 0             # both parents unknown
            F[i+1] = L[i+1] = 0.            # non-inbred
            continue
        end
        if i > 1 && S == sire[i-1] && D == dam[i-1] # full sib with last animal
            F[i+1] = F[i]
            L[i+1] = L[i]
            continue
        end
        F[i+1] = -1.
        L[i+1] = -1.
        t = LAP[i+1]     # longest ancestral path for animal i
        Anc[MI[t]] = i   # initialize Ancestors with this animal
        MI[t] += 1       # increment to compensate for decrement below
        while t ≥ 0
            @show t, MI[t]
            j = Anc[MI[t] -= 1]         # next ancestor
            @show i,j,t
            S = sire[j]                 # parents of the ancestor
            @show S,LAP[S+1]
            if S > 0
                if (L[S] != 0)
                    Anc[MI[LAP[S+1]]] = S # add sire to Ancestor array
                    MI[LAP[S+1]] += 1     # increment minor index for the group
                end
                L[S] += 0.5 * L[j]
            end
            D = dam[j]
            @show D,LAP[D+1]
            if D > 0                    # same for dam
                if (L[D] != 0)
                    Anc[MI[LAP[D+1]]] = D # add sire to Anc
                    MI[LAP[D+1]] += 1
                end
                L[D] += 0.5 * L[j]
            end
            @show B
            @show L
            F[i+1] += abs2(L[j+1]) * B[j]
            @show F
            L[j] = 0.                   # clear L[j] for evaluation of next animal
            if MI[t] == SI[t]           # move to the next LAP group
                t -= 1
            end
       end
    end
    F .+ 1.
end

function incr(LAP::Vector{Int},sire::Vector{Int},dam::Vector{Int},k::Int)
    sl = (S = sire[k]) == 0 ? -1 : (LAP[S] >= 0 ? LAP[S] : incr(LAP,sire,dam,S))
    dl = (D = dam[k]) == 0 ? -1 : (LAP[D] >= 0 ? LAP[D] : incr(LAP,sire,dam,D))
    @show k, S, sl, D, dl
    LAP[k] = max(sl,dl) + 1
end

function editped(sire::Vector{Int},dam::Vector{Int})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    LAP = fill(-1,n)                    # longest ancestor path
    for i in 1:n
        S = sire[i]
        D = dam[i]
        0 <= S <= n && 0 <= D <= n || error("all sire and dam values must be in [0,n]")
        if S == 0 && D == 0
            LAP[i] = 0
        end
    end
    for i in 1:n
        @show i,LAP[i]
        LAP[i] == -1 && incr(LAP,sire,dam,i)
        @show i,LAP[i]
    end
    LAP
end

function orderped(sire::Vector{Int},dam::Vector{Int})
    (n = length(sire)) == length(dam) || throw(DimensionMismatch(""))
    for i in 1:n
        0 ≤ sire[i] ≤ n && 0 ≤ dam[i] ≤ n || error("sire and dam must be in 0:n")
    end
    ord = Int[]
    pop = IntSet(1:n)
    ss = delete!(IntSet(sire),0)
    dd = delete!(IntSet(dam),0)
    parents = union(ss,dd)
    while length(parents) > 0
        append!(ord,collect(setdiff(pop,parents)))
        pp = collect(parents)
        ss = delete!(IntSet(sire[pp]),0)
        dd = delete!(IntSet(dam[pp]),0)
        pop = parents
        parents = union(ss,dd)
    end
    append!(ord,collect(pop))
    reverse(ord)
end

end #module
