using CairoMakie
using Pipe 

"""
taylor_plot(sdmod_rel::Vector{Float64}, correl::Vector{Float64}; bias=0.0, sdObs=1.0,   plotBias=false, plotVarianceErr=false)

taylor_plot([.3,.4,.1],[.66,.77,.11])
"""
function taylor_plot(sdmod_rel::Vector{Float64}, correl::Vector{Float64}; bias=0.0, sdObs=1.0,   plotBias=false, plotVarianceErr=false)
    sdmax = max(maximum(sdmod_rel), 1.0)
    
    ## Function to transform coordinates
    correl_sd2taylor(sd, correl) = (x=correl * sd, y=sqrt(sd^2-(correl*sd)^2))

    ## Create Taylorplot gridlines
    correls = reduce(vcat,[-1, -0.99, -0.95, -0.9:0.1:0.9, 0.95, 0.99, 1.0])
    sds = 0:0.1:sdmax |> collect

    fig, ax, sc = scatter(0,0, color=(:black, 0));
    for c in correls
        lines!(ax, c .* sds .* sdObs * 1.02, sqrt.(sds.^2-(c .* sds).^2).*sdObs * 1.02; linestyle= :dot, color=:grey) 
    end

    for s in sds
        lines!(ax, correls .* s .* sdObs, sqrt.(s.^2 .- (correls .* s).^2).*sdObs; linestyle= s==1 ? :dash : :dot, color=:grey) 
    end

    ## Creat isolines of root mean squared difference (RMSD) and labels
    RMSD = 0.2:0.2:2*sdmax+0.2 |> collect
    xvec =  -1.0:0.01:1 |> collect

    rmsd_lab_x=Vector{Float64}()
    rmsd_lab_y=Vector{Float64}()
    rmsd_lab_text=Vector{String}()

    for r in RMSD[1:end-1]
         
        x2 = @. (xvec * r + 1) * sdObs
        y2 = @. sqrt(1-xvec^2) * r * sdObs
        
       i_plot = [i for i in 1:length(x2) if sqrt(x2[i]^2+y2[i]^2)<=1.005*maximum(sds)*sdObs]

      # This also works, but not if the array sent back by the compr is empty (collect does not work)  
      # xp, yp = zip([(xi,yi) for (xi, yi) in zip(x2, y2) if sqrt(xi^2+yi^2)<=1.005*maximum(sds)*sdObs]...) |> collect
       
        lines!(ax, x2[i_plot], y2[i_plot], color=(:blue, 0.5))

        i_lab=[i for i in 1:length(x2) if (y2[i]>0.3^2*sdObs^2/abs(x2[i])) & (x2[i]>0.0) & (sqrt(x2[i]^2+y2[i]^2)<=1.005*maximum(sds)*sdObs)]
        if length(i_lab)>0
            push!(rmsd_lab_x, x2[i_lab[1]])
            push!(rmsd_lab_y, y2[i_lab[1]]) 
            push!(rmsd_lab_text, "$(round(r, digits=2))")
         end
       
    end
    
    #### Plotting data in Taylor space

    x,y = collect.(zip(correl_sd2taylor.(sdmod_rel, correl)...) |> collect) .* sdObs

    # First the arrows, so that the point is always on top
    if plotBias
        RMSD = @. sqrt(sdObs^2 + (sdmod_rel*sdObs)^2 - 2*sdObs^2*sdmod_rel*correl)
        RMSE = @. sqrt(bias^2+RMSD^2)
        RMSDangle = @. asin(y/RMSD)
        RMSEangle = @. asin(bias/RMSE)
        alpha = @. ifelse(x > sdObs, π - RMSDangle - RMSEangle, RMSDangle - RMSEangle)
        xb = @. sdObs - cos(alpha)* RMSE
        yb = @. sin(alpha) * RMSE
        arrows!(x, y, xb.-x, yb.-y; color=1:length(x), linewidth=3)
    end

    if plotVarianceErr
        x_no_var_err, y_no_var_err = collect.(zip(correl_sd2taylor.(1.0, correl)...) |> collect) .* sdObs
        arrows!(x,y, x_no_var_err.-x, y_no_var_err.-y; color=1:length(x), linewidth=3)
    end

    # if (plotNoVarErr) {
    #     modPointsNoVarErr <- correl_sd2Taylor(1, correl) * sdObs
    #     p  <- p + geom_point(data=modPointsNoVarErr, mapping=aes(x,y), color="red", size=2, alpha=0.5)

    ## Plot data points (original Taylorplot)
    scatter!(ax, x, y, color=1:length(x),  strokewidth=0.5, glowwidth = 5.0, glowcolor=(:black,0.5) )

    ## Plot point of "perfect model" 
    scatter!(ax, sdObs, 0.0, markersize=20, color=(:blue,0.5))
    
  
     
    ### Labels for coordinate system
    correlTicks = @pipe reduce(vcat, [-0.9:0.1:0.9, 0.95, 0.99]) |> setdiff(_, [0])
    text!(correlTicks*sdmax*sdObs*1.04, sqrt.(sdmax^2 .- (correlTicks*sdmax).^2)*sdObs*1.04, 
        text=["$(round(c, digits=2))" for c in correlTicks], align=(:center,:center), fontsize=10)
    #scatter!(correlTicks*sdmax*sdObs*1.04, sqrt.(sdmax^2 .- (correlTicks*sdmax).^2)*sdObs*1.04)
    scatter!(rmsd_lab_x,rmsd_lab_y, markersize=(30,20), color=:white, strokewidth=0.5)
    text!(rmsd_lab_x,rmsd_lab_y, text=rmsd_lab_text, align=(:center,:center), fontsize=10)

    ## Clipping plot
    #xlims!(ax, min(0, 1.1*minimum(x)), max(sdObs,1.1*maximum(x)))
    xlims!(ax, min(0, 1.5*minimum(x)), max(sdObs,1.5*maximum(x)))

    #(fig, ax)
    return fig



    # taylor_grid = @pipe Base.product(correls, sds) |> DataFrame |> rename(_, [:correls, :sds])

    # iso_rmsd = @pipe Base.product(0.2:0.2:2*sdmax+0.2 |> collect, -1.0:0.01:1 |> collect) |> 
    #     DataFrame |> rename(_, [:RMSD, :xvec])

    # @chain iso_rmsd begin
    #     @rtransform! begin
    #         :x2 = (:xvec * :RMSD + 1)* sdObs
    #         :y2 = sqrt(1-:xvec^2)*:RMSD*sdObs
    #     end
    #     @rsubset! sqrt(:x2^2+:y2^2)<=1.005*maximum(taylor_grid.sds)*sdObs
    # end
    
    # @rtransform! taylor_grid :xv = :correls*:sds
    # @rtransform! taylor_grid begin
    #     :yv = sqrt(:sds^2-:xv^2)*sdObs
    #     :xv=:xv*sdObs 
    # end
    # @transform! taylor_grid begin
    # :sdsCat = categorical(:sds)
    # :correlCat = categorical(:correls)
    # end

    # #plt = data(taylor_grid) * visual(Lines, linestyle=:dash) *  mapping(:xv, :yv, group=(:correlCat, :sdsCat))
    # #plt += data(taylor_grid) * visual(Lines, linestyle=:dash) *  mapping(:xv, :yv, group=:sdsCat)
    # #draw(plt; axis=(width = 1225, height = 825))
    # fig=Figure()


   

end