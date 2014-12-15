using pedigree, HDF5
using Base.Test

## simple example from Mrode(2004)
p = Pedigree([0,0,1,1,4,5],[0,0,2,0,3,2])
Ltrans(p)

fid = h5open(Pkg.dir("pedigree","data","pedCows.h5"))
const sire = readmmap(fid["sire"]);
const dam = readmmap(fid["dam"]);
close(fid)

ord,lappt,ss,dd = laporder(sire,dam);
p = Pedigree(ss,dd);
tt = Tinvt(p);
@test all(tt.data.nzval .== -0.5)
@test sort(unique(diff(tt.data.colptr))) == [0,1,2]

Lt = Ltrans(p);

