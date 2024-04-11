# 
function waba()
    #wpth="C:/Users/Public/Documents/Python_Scripts/julia/water-balance.jl"
    wpth=src_path*"/water-balance.jl"
    include(wpth)
    yd=waread("waba-input.wa")|>yrsum
    @info "waba done !"
end

"""
filters internal WaSiM stats of routed discharge files
works recursively
"""
