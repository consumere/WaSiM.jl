# 
function plot_raster_and_points(raster, points)    
        fig = Figure(fontsize=22); #size=(1200, 400),
        #axs = Axis(fig, xlabel="x", ylabel="y"
        fig = Figure();
        axs = Axis(fig, xlabel="", ylabel="",aspect=1)
        Mke.heatmap!(axs,r)
        for pt in ptc_tuples
            Mke.scatter!(axs,pt, markersize = 5, color = :red)
        end
        # Set explicit limits if needed
        #GeoInterface.bbox(ptc_tuples)
        # xmax,ymax = findmax(ptc_tuples)[1]
        # xmin,ymin = findmin(ptc_tuples)[1]
        # limits!(axs, xmin, xmax, ymin, ymax)
        #limits!(axs, 1, size(r)[2], 1,size(r)[1])
        hideydecorations!(axs, grid=false, ticks=true)
        return fig
    end

    """
    lat lon reverse
    """
    