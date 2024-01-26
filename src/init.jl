cd("/mnt/c/Users/chs72fw/.julia/dev/WaSiM/")
using Pkg
Pkg.activate(".")
Pkg.status()
Pkg.test()

force_recompile(package_name::String) = Base.compilecache(Base.identify_package(package_name))
force_recompile("NCDatasets")
Pkg.test()

# ##das müsste dann in dem Modul WaSiM stehen
# include("win/smallfuncs.jl")
# @cmk
#import Pkg; Pkg.add("CairoMakie")

import WaSiM
import WaSiM.du as d
d()
WaSiM.ls()
WaSiM.tree()


cd(src)
import WaSiM.rglob as g
#g("func-w")
g("raster")
g("py")
g("rca")
g("plotly")

# copy dev files to WaSiM.jl
src = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
dst = "/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src"
files = ["./cairomakie.jl","./func-win.jl","./rasterfuncs.jl",
    "./pyjl.jl","./rcall.jl","./RCall_gof.jl","./wajs.jl"]
#"./win/smallfuncs.jl"
function copy_files_with_warnings(src::AbstractString, dst::AbstractString, files::Vector{<:AbstractString})
    for file in files
        #dst_file = joinpath(dst, basename(file))  # Remove leading "./" from file path
        
        src_file = joinpath(src,replace(file,"./"=>""))  
        dst_file = joinpath(dst, replace(file,"./"=>""))  # Remove leading "./" from file path
        
        if isfile(dst_file)
            printstyled("Warning: File $dst_file already exists and will be overwritten.\n",color=:red)
        end
        cp(src_file, dst_file; force=true)
    end
end
copy_files_with_warnings(src, dst, files)


#julia --threads auto -q --startup-file=no --project="." 
using Pkg; Pkg.add("JET")
__precompile__(false)
report_package("WaSiM")

cd("/mnt/c/Users/chs72fw/.julia/dev/WaSiM/")
using Pkg
Pkg.activate(".")
Pkg.status()
__precompile__(false)

#julia --threads auto -q --startup-file=no --project="."
using JET
report_package("WaSiM")

import WaSiM as wa
wa.lat()