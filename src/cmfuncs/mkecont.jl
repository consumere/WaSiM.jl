# 
function mkecont(r::Raster;wasim=true,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            replace!(z, missval=>missing)
        else
            z = r.data
            reverse!(z,dims=1)
            reverse!(z,dims=2)
            replace!(z, missval=>missing)
        end
        ex = extrema(z|>skipmissing)
        lscale = ex[1]:10:ex[2] # Adjusted levels
        fig = Figure(
            #size=(800, 600), 
            fontsize=22);
        axs = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")
        p1 = heatmap!(axs, z, colormap=(:plasma, 0.87))
        contour!(axs, z; color=:black) #, levels=lscale
        Colorbar(fig[1, 2], p1, width=20, ticksize=20, tickalign=1)
        xmax = findmax(size(z))[1] #gets the index of the max value
        limits!(axs, 1, xmax, 1, xmax)
        # hideydecorations!(axs, grid=true, ticks=false)
        # hidexdecorations!(axs, grid=true, ticks=false)
        return fig
    end


    """
    contourplots from rasters 3 plots
    r::Raster;wasim=true,missval=0
    """
    