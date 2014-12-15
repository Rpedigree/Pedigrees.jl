using pedigree, HDF5
using Base.Test

fid = h5open(Pkg.dir("pedigree","data","pedCows.h5"))
const sire = readmmap(fid["sire"]);
const dam = readmmap(fid["dam"]);
close(fid)

ord,lappt,ss,dd = laporder(sire,dam);
p = Pedigree(ss,dd);
tt = Tinvt(p);
@test all(tt.data.nzval .== -0.5)
@test sort(unique(diff(tt.colptr))) = [0,1,2]

Tt = Tmat(p)

