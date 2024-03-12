# 
function pxm(dataframes::Vector{DataFrame}; col=:ve, usethresold=true, threshold=-0.410, all=false, kw...)
        p1=Plots.plot(xrotation = 35)
        # Iterate through each DataFrame in the input
        for df in dataframes
            # Extract year and metrics from the DataFrame
            year = df.year
            if all
                nam = collect(DataFrames.metadata(df))[1][2] |> basename
                nam = replace(nam, "-qoutjl" => "", r"_" => " ", r"#" => "", r"qout" => "")
                if usethresold
                    yaxis!((threshold, 1))
                end
                kge = df.kge
                nse = df.nse
                ve = df.ve
                # Plotting with annotations
                Plots.scatter!(year, kge, label="KGE $nam", marker=:circle)
                for (x, y) in zip(year, kge)
                    if y > threshold
                        Plots.annotate!([(x, y, 
                        Plots.text("$(round(y, digits=2))", 8, 
                        :black, :top))])
                    end
                end
                Plots.scatter!(year, nse, label="NSE $nam", marker=:square)
                for (x, y) in zip(year, nse)
                    if y > threshold
                        Plots.annotate!([(x, y, 
                        Plots.text("$(round(y, digits=2))", 8, 
                        :black, :top))])
                    end
                end
                Plots.scatter!(year, ve, label="VE $nam", marker=:diamond, legend=:outertopright, kw...)
                for (x, y) in zip(year, nse)
                    if y > threshold
                        Plots.annotate!([(x, y, 
                        Plots.text("$(round(y, digits=2))", 8, 
                        :black, :top))])
                    end
                end
                title!("Grouped Metrics by Year")
                Plots.xticks!(df[:, :year])
                ylabel!("Score")
            else
                if usethresold
                    yaxis!((threshold, 1))
                end
                
                ve = select(df, col) |> Matrix |> vec
                ti = select(df, col) |> names |> first |> uppercase
                nam = collect(DataFrames.metadata(df))[1][2] |> basename
                nam = replace(nam, r"-qoutjl" => "", r"_" => " ", r"#" => "", r"qout" => "")
                
                # Plotting with annotations
                Plots.scatter!(year, ve, 
                    #label=nam, 
                    label=false,
                    title=ti, 
                    marker=:diamond, legend=:outertopright, kw...)
                
                # Adding annotations for values above threshold
                for (x, y) in zip(year, ve)
                    if y > threshold
                        Plots.annotate!([(x, y, 
                        Plots.text("$(round(y, digits=2)) "*nam[1:3], 
                        8, :black, :top))])
                    end
                end
                
                Plots.xticks!(df[:, :year])
            end
            # xx = uppercase(string(col))
            # ylabel!("Score $xx")
        end
        return p1
    end
    


    """
    merges and saves to qoutjl
    outd = pout(infile;sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/",ofl = "route.txt")
    """
    