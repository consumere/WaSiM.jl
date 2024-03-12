# 
function toMain()
        fnames = names(Main.cmk, all=true)
        for submodule in fnames
            @eval import Main.cmk.$submodule
        end
    end

    # using DataFrames, CSV, Statistics, Dates, Distributions

    """
    find LOG. R-SQUARE > .4 recursivley
    """
    