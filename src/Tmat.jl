## return the right Cholesky factor from the LDLt decomposition of A⁻¹
function Tinvt{T}(p::Pedigree{T})
    n = length(p.sire)
    sp = convert(Vector{T},find(p.sire .!= 0))
    dp = convert(Vector{T},find(p.dam .!= 0))
    Triangular(sparse([p.sire[sp],p.dam[dp]],[sp,dp],-0.5,n,n),:U,true)
end

## return the right Cholesky factor, Lt, from the LLt factorization of A
function Ltrans{T}(p::Pedigree{T})
    n = length(p.sire)
    sp = convert(Vector{Int32},find(p.sire .!= 0))
    dp = convert(Vector{Int32},find(p.dam .!= 0))
    tri = sparse(convert(Vector{Int32},[p.sire[sp],p.dam[dp]]),[sp,dp],1.0,n,n)
    Tv = eltype(tri)
    cpi = tri.colptr
    rvi = tri.rowval
    nvi = tri.nzval
    cpl = sizehint(Int32[1],n+1)             # column pointers for L'
    rvl = sizehint(Int32[], 8n)              # row values for L'
    nvl = sizehint(Tv[], 8n)               # nonzero values in L'
    anc = sizehint(IntSet([]),n)           # set of ancestors 
    for j in 1:n
        empty!(anc)
        l2parents = zero(Tv)   # accumulator for sum of squares of L cols of parents
        for k in cpi[j]:(cpi[j+1]-1)    # each parent in the pedigree
            for i in cpl[rvi[k]]:(cpl[rvi[k]+1]-1)
                push!(anc,rvl[i])       # add the parents ancestors to the set
                l2parents += abs2(nvl[i])
            end
        end
        push!(anc,j)                    # the animal counts as one of its ancestors
        push!(cpl,cpl[end] + length(anc))
        rv = collect(anc)
        append!(rvl,rv)
        lrv = length(rv)
        cc = zeros(eltype(tri),(lrv,))
        cc[lrv] = sqrt(one(Tv) - convert(Tv,0.25)*l2parents)
        for k in cpi[j]:(cpi[j+1]-1)    # each parent in the pedigree
            for i in cpl[rvi[k]]:(cpl[rvi[k]+1]-1)  # add 0.5*each parent's column
                cc[searchsortedfirst(rv,rvl[i])] += convert(Tv,0.5)*nvl[i]
            end
        end
        append!(nvl,cc)
    end
    SparseMatrixCSC(n,n,cpl,rvl,nvl)
end

function inbreeding(p::Pedigree)
    Lt = Ltrans(p)
    inb = zeros(Lt.n)
    nvl = Lt.nzval
    cpt = Lt.colptr
    for i in 1:Lt.n
        for k in cpt[i]:(cpt[i+1]-1)
            inb[i] += abs2(nvl[k])
        end
        inb[i] -= one(eltype(inb))
    end
    inb,Lt
end
