function Tinv{T}(p::Pedigree{T})
    n = length(p.sire)
    sp = convert(Vector{T},find(p.sire .!= 0))
    dp = convert(Vector{T},find(p.dam .!= 0))
    Triangular(sparse([sp,dp],[p.sire[sp],p.dam[dp]],-0.5,n,n),:L,true)
end
