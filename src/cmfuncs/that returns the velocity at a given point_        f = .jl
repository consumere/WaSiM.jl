# 
function that returns the velocity at a given point
        f = (x, y) -> Point2(gx[round(Int, y), round(Int, x)], gy[round(Int, y), round(Int, x)])
        #f = (x, y) -> Point2([gx[y, x], gy[y, x]])
    
        # Define the x and y intervals
        xinterval = 1:size(z, 2)
        yinterval = 1:size(z, 1)       
        
        # #lscale = ex[1]:10:ex[2] # Adjusted levels
        # fig = Figure(
        #     #size=(800, 600), 
        #     fontsize=22);
        # axs = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")
        # xmax = findmax(size(z))[1] #gets the index of the max value
        # # Create the streamplot
        p1 = CairoMakie.streamplot(f, xinterval, yinterval)
        #Colorbar(fig[1, 2], p1, width=20, ticksize=20, tickalign=1)
        #limits!(axs, 1, xmax, 1, xmax)
        #return fig
        return p1
    end

end    ##end of module

println("cairomakie.jl module cmk loaded")

# #https://docs.makie.org/stable/explanations/fonts/
# Makie.to_font("Computer Modern")
# Makie.theme(:fonts)
# #https://docs.makie.org/stable/explanations/latex/
# using CairoMakie
# with_theme(theme_latexfonts()) do
#     fig = Figure()
#     Label(fig[1, 1], "A standard Label", tellwidth = false)
#     Label(fig[2, 1], L"A LaTeXString with a small formula $x^2$", tellwidth = false)
#     Axis(fig[3, 1], title = "An axis with matching font for the tick labels")
#     fig
# end
