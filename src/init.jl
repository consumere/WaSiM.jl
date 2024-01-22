#@profile
cd("/mnt/c/Users/chs72fw/.julia/dev/WaSiM/")
using Pkg
Pkg.activate("WaSiM.jl")
Pkg.status()
Pkg.add(PackageSpec(path = "."))
