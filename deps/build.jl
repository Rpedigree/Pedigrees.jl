using BinDeps

@BinDeps.setup

libpedigree = library_dependency("libpedigree")

depsdir = BinDeps.depsdir(libpedigree)
prefix = joinpath(depsdir,"usr")
srcdir = joinpath(depsdir,"src")
libdir = joinpath(prefix,"lib")
CreateDirectory(libdir)

provides(SimpleBuild,
         (@build_steps begin
             ChangeDirectory(srcdir)
             (@build_steps begin
                 `make`
                 `make install`
             end)
         end),[libpedigree])

@BinDeps.install Dict(:libpedigree => :libpedigree)
