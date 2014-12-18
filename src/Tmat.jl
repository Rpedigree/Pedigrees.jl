## return the right Cholesky factor from the LDLt decomposition of A⁻¹
function Tinvt{T}(p::Pedigree{T})
    n = length(p.sire)
    cpi = Array(Int32,n + 1)            # column pointers
    cpi[1] = one(Int32)
    rvi = sizehint(Int32[],n + n>>1)    # row values of T⁻¹
    for i in 1:n
        r1 = min(p.sire[i],p.dam[i])
        r2 = max(p.sire[i],p.dam[i])
        if r2 == 0
            cpi[i+1] = cpi[i]
        elseif r1 == 0
            cpi[i+1] = cpi[i] + one(Int32)
            push!(rvi,r2)
        else
            cpi[i+1] = cpi[i] + convert(Int32,2)
            push!(rvi,r1)
            push!(rvi,r2)
        end
    end
    Triangular(SparseMatrixCSC(n,n,cpi,rvi,fill(-0.5,(length(rvi,)))),:L,true)
end

## return the right Cholesky factor, Lt, from the LLt factorization of A
function Ltrans{T}(p::Pedigree{T})
    n = length(p.sire)
    tri = Tinvt(p).data
    cpi = tri.colptr
    rvi = tri.rowval
    cpl = sizehint(Int32[1],n+1)             # column pointers for L'
    rvl = sizehint(Int32[], 8n)              # row values for L'
    nvl = sizehint(Float64[], 8n)            # nonzero values in L'
    anc = sizehint(IntSet([]),n)             # set of ancestors
    for j in 1:n
        empty!(anc)
        l2parents = 0.0   # accumulator for sum of squares of L cols of parents
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
        cc[lrv] = sqrt(1.0 - 0.25*l2parents)
        for k in cpi[j]:(cpi[j+1]-1)    # each parent in the pedigree
            for i in cpl[rvi[k]]:(cpl[rvi[k]+1]-1)  # add 0.5*each parent's column
                cc[searchsortedfirst(rv,rvl[i])] += 0.5*nvl[i]
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
