# 
function plot_grouped_metrics(dataframes::Vector{DataFrame};col=:ve,usethresold=true,threshold = -0.41,all=false,kw...)
        #nam = getnames(dataframes)
        # Create a new plot
        Plots.plot(xrotation = 35)

        # Iterate through each DataFrame in the input
        for df in dataframes
            # Extract year and metrics from the DataFrame
            year = df.year
            if all
                nam=collect(DataFrames.metadata(df))[1][2]|>basename
                nam=replace(nam,"-qoutjl"=>"",r"_" => " ",r"#" => "",r"qout" => "")
                # if usethresold
                #     kge = ifelse.(df.kge .<= threshold, -1, df.kge)
                #     nse = ifelse.(df.nse .<= threshold, -1, df.nse)
                #     ve = ifelse.(df.ve .<= threshold, -1, df.ve)
                #     yaxis!((threshold,1))
                # else
                #     kge = df.kge
                #     nse = df.nse
                #     ve = df.ve
                # end
                if usethresold
                    yaxis!((threshold,1))
                end
                kge = df.kge
                nse = df.nse
                ve = df.ve
                # Plot each metric with a different color
                plot!(year, kge, label="KGE $nam", marker=:circle)
                plot!(year, nse, label="NSE $nam", marker=:square)
                plot!(year, ve, label="VE $nam", 
                    marker=:diamond,
                    legend=:outertopright,
                    kw...)
                title!("Grouped Metrics by Year")
            else
                #dfilter = DataFrames.subset(df, col => ByRow(<=(threshold)))|>Matrix|>vec
                #DataFrames.transform(df, names(df)[2:end-1] .=> cor,renamecols=false)
                # if usethresold
                #     df[!, col] .= ifelse.(df[!, col] .<= threshold, -1, df[!, col])
                #     yaxis!((threshold,1))
                # end
                if usethresold
                    yaxis!((threshold,1))
                end
                ve = select(df,col)|>Matrix|>vec
                ti = select(df,col)|>names|>first|>uppercase
                nam=collect(DataFrames.metadata(df))[1][2]|>basename
                nam=replace(nam,r"-qoutjl"=>"",r"_" => " ",r"#" => "",r"qout" => "")
                plot!(year, ve, label=nam, 
                    title =ti, 
                    marker=:diamond,
                    legend=:outertopright,
                    kw...)
                Plots.xticks!(df[:, :year])
            end
        end
        # Customize the plot
        #xlabel!("Year")
        
        xx = uppercase(string(col))
        ylabel!("Score $xx")
    end

    """
    takes Vector{DataFrame} from:
    dataframes = map(byear,outd)
    annoated scatterplot of metrics
    pxm(dataframes::Vector{DataFrame}; col=:ve, usethresold=true, threshold=-0.410, all=false, kw...)
    """
    