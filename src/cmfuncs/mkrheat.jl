# 
function mkrheat(x::Union{String,Regex,Raster};msk=true,missval=0.0,umask=10^6,mskval=0.001,layer=1,crs=EPSG(25832),mappedcrs=EPSG(25832),turn=true,kw...)    
        if isa(x,Regex)
            #fn = filter(z->endswith(z,"nc"),glob(x))[1]
            #ne = join([x,"nc","\$"],"+.")
            ds = try
                #fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]
                fn = first(filter(z -> endswith(z, "nc") && occursin(x, z), readdir()))
                Raster(fn;missingval=missval,kw...) #missingval,
            catch
                @error "no file found!"
                return
            end
            
        elseif isa(x,Raster)
            ds = x
            fn = string.(name(ds))
        else
            fn = x
            ds = Raster(fn;missingval=missval,kw...)
        end
        
        #A = ds[Ti=layer].data
        #A = ds[:, :, 1].data
        #mx = missingmask(ds)
            
        ds = ds[:, :, layer]
        if msk
            if umask==10^6
                umask = maximum(ds|>skipmissing)
            end
            #A = ds.data .> float(mskval)
            #ds = Rasters.rebuild(ds;missingval=0.0)
            #ds = replace_missing(ds,missingval)
            #bitmat = (dx .> msk) .& (dx .<= umask)
            zm = (ds .> float(mskval)) .& (ds .<= umask)
            #zm = ds .> float(mskval)
            #zm = boolmask(zm;missingval=0.0)
            zm = boolmask(zm;missingval=missval)
            
            A = Rasters.mask(ds; with=zm)
            min_val, max_val = extrema(skipmissing(A))
            @info join(["min: ", string(min_val), ", max: ", string(max_val)], " ")
            try
                @info A.metadata.val    #|>DataFame
            catch
                @warn "no metadata in raster!"
            end
            #A = Rasters.mask(ds; with=boolmask(ds;missingval=mskval))
            A = A.data
            #A = Rasters.rebuild(A;kw...)
            #A = ds.data
            #A = replace(A, missing => mskval)
            #A = A[A .> mskval]
        else
            A = ds.data
        end
        #see also Rasters.trim
        #extrema(A|>skipmissing)
        # A = replace(A, -9999.0 => missing)
        # A = replace(A, 0.0 => missing)
        #b = replace(A, [0 .. 70] => missing)
        #b = replace(A, float(-10..70) => missing)
        #b = map(x -> replace(x,float(-10..70) => missing),eachrow(A))
        

        mrows = map(x -> !all(ismissing, x),eachrow(A))
        mcols = map(x -> !all(ismissing, x),eachcol(A))
        A = A[mrows, mcols]
        
        if turn
            A = transpose(A)
            A = reverse(A, dims=2) #very important!
        else
            A = replace(A, missing => 0.0)
        end
        
        ti = replace(basename(fn),".nc"=>"","_"=>" ")
        #ti = replace(name(ds),".nc"=>"","_"=>" ")
        #https://docs.makie.org/stable/explanations/colors/
        fig = Figure(
            #color = Makie.wong_colors(),
            #color = :thermal, #see colormap below
            #size=(800, 600), 
            fontsize=22);
        axs = Axis(fig[1,1],
            title=ti, 
            aspect = DataAspect(), 
            xlabel="x", 
            ylabel="y")
        

        #p1 = heatmap!(A,colormap=(:turbo, 0.9))
        p1 = heatmap!(A,colormap=(:thermal, 0.9))
        #p1 = heatmap!(A,colormap=Makie.wong_colors())
        
            contour!(axs, A; color=:black) #, levels=lscale
            Colorbar(fig[1, 2], p1, 
                width=20, ticksize=20, 
                tickalign=1)
            #xmax,ymax = size(A)
            #limits!(axs, 1, xmax, 1, ymax)
        # hideydecorations!(axs, grid=true, ticks=false)
        # hidexdecorations!(axs, grid=true, ticks=false)
        return fig
    end

    """
    reads ascii grid and replaces -9999 \n returns Raster
    """
    