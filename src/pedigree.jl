module pedigree

## if isfile(joinpath(Pkg.dir("pedigree"),"deps","deps.jl"))
##     include("../deps/deps.jl")
## else
##     error("pedigree not properly installed. Please run Pkg.build(\"pedigree\")")
## end

    export Pedigree, laporder, Tinvt, Ltrans

    include("types.jl")
    include("inbreeding.jl")
    include("Tmat.jl")

end #module
