# 
function mkheat(x::Union{String,Regex};msk=true,layer=0)
    #     if isa(x,Regex)
    #         #fn = filter(z->endswith(z,"nc"),glob(x))[1]
    #         #ne = join([x,"nc","\$"],"+.")
    #         fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]        
    #     else
    #         fn = x
    #     end
    #     ds = xr.open_dataset(fn)
    #     k = ds.keys()|>collect|>first
    #     A = ds[k].isel(t=layer).squeeze()[:values]
    #     A = reverse(A, dims=2) #very important!
    #     if msk
    #         A = map(x -> x <= 0 ? NaN : x, A)
    #     end
    #     fig = Figure(
    #         #size=(800, 600), 
    #         fontsize=22);
    #     ti = replace(fn,".nc"=>"","_"=>" ")
    #     axs = Axis(fig[1,1],title=ti, aspect=1, 
    #         xlabel="x", 
    #         ylabel="y")
    #     #replace NaN values with missing
    #     B = replace(A, NaN => missing)
    #     # Find rows and columns where all values are missing
    #     xm = .!all(ismissing, B, dims=2)
    #     ym = .!all(ismissing, B, dims=1)
    #     # Remove rows and columns where all values are missing
    #     A = B[xm[:], ym[:]]

    #     p1 = heatmap!(A,colormap=(:turbo, 0.9))
    #         contour!(axs, A; color=:black) #, levels=lscale
    #         Colorbar(fig[1, 2], p1, 
    #             width=20, ticksize=20, 
    #             tickalign=1)
    #         xmax,ymax = size(A)
    #         limits!(axs, 1, xmax, 1, ymax)
    #     # hideydecorations!(axs, grid=true, ticks=false)
    #     # hidexdecorations!(axs, grid=true, ticks=false)
    #     return fig
    # end


    """
    mkrheat(x::Union{String,Regex};msk=true,umask=10^6,mskval=0.001,layer=1) \n
    uses Rasters to read netcdf and makie to plot.
    turn=true flips the raster
    mkrheat(x::Union{String,Regex,Raster};msk=true,missval=0.0,umask=10^6,mskval=0.001,layer=1,crs=EPSG(25832),mappedcrs=EPSG(25832),turn=true,kw...)
    """
    