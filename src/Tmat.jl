## return the right Cholesky factor from the LDLt decomposition of A⁻¹
function Tinvt{T}(p::Pedigree{T})
    n = length(p.sire)
    sp = convert(Vector{T},find(p.sire .!= 0))
    dp = convert(Vector{T},find(p.dam .!= 0))
    Triangular(sparse([p.sire[sp],p.dam[dp]],[sp,dp],-0.5,n,n),:U,true)
end

## return the right unit Cholesky factor, Lt, from the LDLt decomposition of A
function Tmat{T}(p::Pedigree{T})
    tri = Tinvt(p).data                 # strict upper triangle of transpose of T⁻¹
    n = tri.n
    Tv = eltype(tri)
    cpi = tri.colptr
    rvi = tri.rowval
    nvi = tri.nzval
    cpl = sizehint(Int[1],n+1)             # column pointers for L'
    rvl = sizehint(Int[], 8n)              # row values for L'
    nvl = sizehint(eltype(tri)[], 8n)      # nonzero values in L'
    anc = sizehint(IntSet([]),n)           # set of ancestors 
    for j in 1:n
        empty!(anc)
        for k in cpi[j]:(cpi[j+1]-1)    # each parent in the pedigree
            kp = rvi[k]                 # parent number
            union!(anc,rvl[cpl[kp]:(cpl[kp+1]-1)])
        end
        push!(anc,j)                    # add
        push!(cpl,cpl[end] + length(anc))
        rv = collect(anc)
        append!(rvl,rv)
        lrv = length(rv)
        cc = zeros(eltype(tri),(lrv,))
        cc[lrv] = one(eltype(tri))
        for k in cpi[j]:(cpi[j+1]-1)    # each parent in the pedigree
            kp = rvi[k]                 # parent number
            for i in cpl[kp]:(cpl[kp+1]-1)  # add 0.5*each parent's column
                cc[searchsortedfirst(rv,rvl[i])] += 0.5*nvl[i]
            end
        end
        append!(nvl,cc)
    end
    SparseMatrixCSC(n,n,cpl,rvl,nvl)
end
