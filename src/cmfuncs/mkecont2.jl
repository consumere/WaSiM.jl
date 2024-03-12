# 
function mkecont2(r::Raster;wasim=true,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            try
                replace!(z, missval=>missing)    
            catch
                @error "MethodError: Cannot `convert` an object of type Missing to an object of type Float64"
                @warn "replace missing failed!"
            end
            
        else
            z = r.data
            reverse!(z,dims=1)
            reverse!(z,dims=2)
            try
                replace!(z, missval=>missing)    
            catch
                @error "MethodError: Cannot `convert` an object of type Missing to an object of type Float64"
                @warn "replace missing failed!"
            end
        end
        ex = extrema(z|>skipmissing)
        lscale = ex[1]:10:ex[2] # Adjusted levels
        fig = Figure(size=(1200, 400), fontsize=22);
        axs = [Axis(fig[1, j], aspect=1, xlabel="x", ylabel=j == 1 ? "y" : "")
                for j in 1:3]
        p1 = heatmap!(axs[1], z, colormap=:plasma, levels=lscale)
            contour!(axs[2], z; color=:black) 
            heatmap!(axs[3], z; colormap=(:plasma, 0.85))
        contour!(axs[3], z; color=:white)
        Colorbar(fig[1, 4], p1, width=20, ticksize=20, tickalign=1)
        xmax = findmax(size(z))[1] #gets the index of the max value
        [limits!(axs[i], 1, xmax, 1, xmax) for i in 1:3]
        [hideydecorations!(axs[i], grid=false, ticks=false) for i in 2:3]
        return fig
    end

    