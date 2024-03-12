# 
function gz(prefix::AbstractString)
#     rootdir="."
#     results = []
#     if (any(x->isdir(x),readdir()))
#         for (looproot, dirs, filenames) in walkdir(rootdir)
#             for filename in filenames
#                 #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                 if (occursin(Regex(prefix,"i"),filename))
#                     push!(results, joinpath(looproot, filename)) 
#                 end
#             end
#         end
#     else
#         printstyled("no dirs in $rootdir !\n",color=:light_red)
#         for filename in (filter(x->isfile(x),
#             readdir(;join=false)))
#             #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#             if (occursin(Regex(prefix,"i"),filename))
#                 push!(results, filename) 
#             end
#         end
#     end
#     return results
# end

# gz("s")
# cd("b1")
# gz("wind")

# ##in NEW session:
# using PackageCompiler
# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "win_ts_plots.so")
# # joinpath(pt,nm)
# create_sysimage([:StatsPlots,
#     :Plots,:PlotThemes,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))
#julia --startup-file=no --sysimage "C:/Users/chs72fw/.julia/sysimages/win_ts_plots.so" C:\Users\Public\Documents\Python_Scripts\julia\waba-rev.jl
#C:/Users/Public/Documents/Python_Scripts


# using PackageCompiler
# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "win_makie.so")
# create_sysimage([
#     :CairoMakie,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))


# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "windf.so")
# create_sysimage([
#     :Grep,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))

#######################
# sqrt(a) === a^0.5
# log(a^0.5) === 0.5*log(a)
# log(a^2) === 2*log(a)
# m,n = 20,15
# log(m * n) === log(m) + log(n)
# log(m ÷ n) === log(m) - log(n) #f
# log(m / n) === log(m) - log(n)  #but rounded, its true
# round(log(m / n),digits=3) === round(log(m) - log(n),digits=3)  #true

# cd(raw"J:\jras")
# using PackageCompiler
# create_sysimage(;sysimage_path="all.so")
