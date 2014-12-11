# evaluate the inbreeding coefficients in a pedigree
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
    
    F[1] = -1.                      # unknown parent
    L = zeros(n+1)
    B = zeros(n+1)                  # diagonal factor in A
    Anc = zeros(Int,n+1)            # Ancestor
    maxlap = maximum(lap)
    SI = zeros(Int,maxlap+1)        # start indices
    MI = zeros(Int,maxlap+1)        # minor indices
    for i in 1:n
        S = sire[i]
        D = dam[i]
        B[i] = 0.5 - 0.25*(F[S+1]+F[D+1])
        for j in 1:lap[i+1]         # adjust start and minor
            SI[j+1] += 1
            MI[j+1] += 1
        end
        @show SI
        @show MI
        if S == 0 && D == 0         # both parents unknown
            F[i+1] = L[i+1] = 0.    # non-inbred
            continue
        end
        if i > 1 && S == sire[i-1] && D == dam[i-1]
            F[i+1] = F[i]           # full sib with last animal
            L[i+1] = L[i]
            continue
        end
        F[i+1] = -1.
        L[i+1] = 1.
        t = lap[i+1]       # longest ancestral path length
        @show i,S,D,t
        @show F
        @show L
        @show SI
        @show MI
        Anc[MI[t+1]+1] = i   # initialize Ancestors with this animal
        MI[t+1] += 1 # increment to compensate for decrement below
        while t ≥ 0
            @show t, MI[t+1]
            MI[t+1] -= 1
            j = Anc[MI[t+1]+1]   # next ancestor
            @show i,j,t
            S = sire[j]             # parents of the ancestor
            @show S,lap[S+1]
            if S > 0
                if (L[S+1] != 0)
                    @show MI[lap[S+1]+1]
                    Anc[MI[lap[S+1]+1]+1] = S # add sire to Ancestor array
                    MI[lap[S+1]+1] += 1     # increment minor index for the group
                end
                L[S+1] += 0.5 * L[j+1]
            end
            D = dam[j]
            @show D,lap[D+1]
            if D > 0                # same for dam
                if (L[D+1] != 0)
                    Anc[MI[lap[D+1]+1]+1] = D
                    MI[lap[D+1]+1] += 1
                end
                L[D+1] += 0.5 * L[j+1]
            end
            @show B
            @show L
            @show Anc
            F[i+1] += abs2(L[j+1]) * B[j]
            @show F
            L[j+1] = 0. # clear L[j] for evaluation of next animal
            if MI[t+1] == SI[t+1]           # move to the next LAP group
                t -= 1
            end
        end
        end
    F .+ 1.
end
