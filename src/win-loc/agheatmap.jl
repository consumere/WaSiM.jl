if length(ARGS) == 0
	println("need args! <ncfileregex>...")
    exit()
end

import ArchGDAL
using Plots: heatmap,annotate!,display,text,savefig

function nconly(x1::AbstractString)
    v::Vector{String} = readdir();
    v = v[broadcast(x->endswith(x,"nc"),v)];
    z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    return(z)
end

function agheat(s::AbstractString;step=30,lyr=1,
    msk=0.001,umask=10^6,roundto=2)
    """
    AG only.
    import ArchGDAL as AG errors...
    plots heatmap from raster with ArchGDAL
    """
    if !isfile(s)
        error("file not found!")
    end
    @info("processing $s ...")
    r = ArchGDAL.readraster(s)
    dx = ArchGDAL.getband(r,lyr)
    if umask==10^6
        umask = maximum(dx)
    end
    println("MIN:",minimum(r),"\nMAX:",maximum(r))
    println("masking value range to $msk and $umask ...")
    bitmat = (dx .> msk) .& (dx .<= umask)
    # Filter the dx matrix using bitmat
    dx_filtered = dx .* bitmat
    dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
    dx_output .= NaN
    dx_output[bitmat] .= dx_filtered[bitmat]
    dx = dx_output
    if !endswith(s,".nc")
        @warn("file seems not to be a NetCDF!")
        dx = permutedims(dx, (2, 1))
        dx = reverse(dx, dims=1)
    else
        @info("processing heatmap ...")
        dx = reverse(dx, dims=2)
        dx = reverse(dx, dims=1)
    end   
    heatmap_plot = heatmap(dx, c=:lightrainbow,
                    title=basename(s), 
                    xaxis="", yaxis=""; dpi=300, size=(800, 400));
    #step = Int(ceil(size(dx)[2]/4))
    #step = Int(ceil(prod(size(dx))/400))
    step = Int(ceil(sum(size(dx))/4^2))
    for i in 1:step:size(dx, 1)
        for j in 1:step:size(dx, 2)
            value = round(dx[i, j]; digits=roundto)
            color = isnan(value) ? :white : :black
            annotate!(j, i, text(string(value), 7, color, :center, 
                halign=:center, rotation=-35.0))
        end
    end
    # Save the plot
    outfile = replace(s,".nc"=>"-heatmap.png")
    savefig(heatmap_plot, outfile)
    @info( "$outfile saved!")
end


fn = ARGS[1];
x=nconly(fn)|>last
#agheat(x;step=100)
#s=raw"D:\Wasim\regio\out\c1\sb1_rcm_1500.mit.nc"
agheat(x)