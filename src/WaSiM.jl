
#julia --startup-file=no -q --color=yes --project="/mnt/c/Users/chs72fw/.julia/dev/WaSiM"

#cd "/mnt/c/Users/chs72fw/.julia/dev/WaSiM"
#julia --startup-file=no -q --color=yes --project=.

#-e 'using Pkg; Pkg.instantiate(); Pkg.API.precompile()'
#include("src/jl")
# cd(raw"C:\Users\chs72fw\.julia\dev\WaSiM")
# pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM";cd(pt)

# using DataFrames, CSV, Statistics, Dates, StatsPlots, Distributions,DataFramesMeta
# using DelimitedFiles, Grep, Printf, PrettyTables, Rasters
# import NCDatasets
# import ArchGDAL
# import GeoInterface
# import GeoDataFrames
# import Shapefile
# import InteractiveUtils
# using Plots.PlotMeasures
# using SHA
# using PyCall
# import Conda

#import Main #for toMain

# function toMain()
#     fnames = names(Main.WaSiM, all=true)
#     for submodule in fnames
#         @eval import Main.WaSiM.$submodule
#     end
# end

#import StatsPlots.@df

# module kat
#     function xz()
#         println("hi from kat")
#     end
# end
#using StatsPlots

   
    # function toexp()
    # fnames = names(WaSiM, all=true)
	# for submodule in fnames
	    # export WaSiM.$submodule
	# end
    # end

#__precompile__(false)

module WaSiM
    #using Reexport
    #@reexport 
    using DataFrames, StatsPlots, Dates
    #import DataFrames
    #@reexport using DataFramesMeta, CSV, Statistics, Dates, StatsPlots, Distributions    
    using DataFramesMeta, CSV, Statistics, Distributions    
    import Conda
    # import ArchGDAL
    # import GeoDataFrames
    # import GeoInterface
    # import InteractiveUtils
    # import NCDatasets
    # import Shapefile
    # import StatsPlots:@df
    import DelimitedFiles:readdlm
    using Grep
    using KernelDensity
    using Plots.PlotMeasures
    using PrettyTables
    using Printf
    using Rasters
    using SHA #for python deps & condasize
    
    ##import rasterstuff
    include("rasterfuncs.jl")
    src_path = "./src"
    include("smallfuncs.jl")
    include("timeseries.jl")
    #@reexport using smfc #no Pkg!   

    function toMain()
    fnames = names(WaSiM, all=true)
        for submodule in fnames
            @eval import WaSiM.$submodule
        end
    end
    
    export toMain

end #endof module

println("used Threads: ", Threads.nthreads())
