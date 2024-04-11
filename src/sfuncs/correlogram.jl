# 
function correlogram(df)
        rows = cols = size(df,2)
        plots = []
        for row = 1:rows, col = 1:cols
            if row == col
                push!(
                    plots,
                    histogram(df[:,row],bins=10, xtickfont = font(5), ytickfont = font(5), legend = false))
            else
                push!(
                    plots,
                    scatter(df[:,row], df[:,col], xtickfont = font(5), ytickfont = font(5), legend = false, markersize=1, alpha = 0.3, smooth = true,
                    linewidth=3, linecolor=:red),
                )
            end
        end
        Plots.plot(plots..., #size=(1200, 1000), 
            layout = (rows, cols))
    end

    """
    fillmissing(df,fillval=-9999)
    """
    