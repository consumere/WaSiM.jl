#rasterfuncs

#module rst
using Reexport
@reexport using WaSiM
#@reexport 
using DataFrames, CSV, Statistics, Dates, Distributions,StatsPlots, Plots.PlotMeasures
using DelimitedFiles, Grep, Printf, PrettyTables
using Rasters, ArchGDAL 

import ArchGDAL
import GeoDataFrames
import GeoInterface
import InteractiveUtils
import NCDatasets
import Shapefile


module rst

    using DataFrames, CSV, Statistics, Dates, Distributions,StatsPlots, Plots.PlotMeasures
    using DelimitedFiles, Grep, Printf, PrettyTables
    using Rasters, ArchGDAL 
    import NCDatasets, Shapefile
    import LibGEOS: centroid
    import GeoInterface: coordinates
    import GeoInterface
    import CFTime
    import GeoDataFrames
    import Conda #for gdal

    ##this imports to wa. see below
    # export addplot, agcont, agcont2, agheat, agjson, 
    # agmask, agread, agsurf, cntplt, cpal, cpl, cplt, 
    # crop_geodf, descr, ezplot, facets, facets_loop, 
    # filterplot, filterplot!, fp, fzplot, fzplot2, gdalop, 
    # jlzonal, jlzonal2, loadalloutput, mask_trim, 
    # maskplot, median_filter, ncdf, ncell, ncmean, 
    # process_rasters, project, project_o, rcont, 
    # readalloutput, readallras, readras, readras2, 
    # readrasrec, remocont, revcrds, reverse_coords, 
    # rhist, rlpro, rmask, rp, rp3, rpall, rplot, rpm, 
    # rpmcf, rpr, rst, rstats, rsubset, serplot, stackplot, 
    # stats, stp, stplot, stread, surf, tdifnc

    #begin

        # """
        # greps from current dir iRegex
        # """
        # function glob(x::AbstractString)
        #     filter(file -> occursin(Regex(x,"i"),file), readdir())
        # end

        # function glob(x::Regex)
        #     """
        #     greps from current dir Regex
        #     """
        #     filter(file -> occursin(x,file), readdir())
        # end

        
        """
        3D plot with geoarrays
        depricated!
        use rpr, surf or agsurf instead!
        """
        function rp3(x::String)
            @error "
            depricated!
            use rpr , surf or agsurf instead!"
            # ga = GeoArrays.read(x)
            # values = ga.A # a 3D array of raster values
            # #GeoArrays.coords(ga) # a tuple of x, y and band coordinates
            # #crs = ga.crs # a string of CRS definition
            # t = GeoArrays.coords(ga)|>size
            # coords = (1:t[1], 1:t[2]) # a Tuple{UnitRange{Int64}, UnitRange{Int64}}
            # ti=basename(x)     #title!("3D Raster Plot")
            # #p1=
            # Plots.surface(coords[1], coords[2], 
            # values[:, :, 1]    ,
            # xlabel="x",ylabel="y",zlabel="value",title=ti)
            
        end


        """
        reads non-recursively all ncdf files
        """
        function readallras(path::AbstractString)
            v = readdir(path);
            v = v[broadcast(x->endswith(x,"nc"),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        function readallras(path::AbstractString, ex::AbstractString)
            v = readdir(path);
            v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        function readallras(ex::Regex)
            v = readdir(".");
            v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        """
        reads timeseries and stores to Vector{DataFrame}
        reads NetCDFs and stores to Vector{Any}

        usage: 
        dfs,ncs = readalloutput()

        errors if not rmeq()

        """
        function readalloutput(;cwd = ".")    
            dfs=loadalldfs(cwd)
            ncs=readallras(cwd)
            return(dfs,ncs)
        end

        loadalloutput = readalloutput

        """
        rplot(x::Regex)
        reads first match of regex wasim ncs
        """
        function rplot(reg::Regex; lyr::Int=1)
            file = filter(x -> occursin(reg,x), readdir(pwd()))
            println("subsetting first nc of $file...")
            file = filter(x -> endswith(x,".nc"), file)|>first
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
        end

        """
        reads all into RAM
        readras(file::AbstractString,missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832);kw...)
        """
        function readras(file::AbstractString,missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832);kw...)
            x=read(Raster(file,missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832);kw...))
            #describe(x)
            return(x)
        end

        """
        reads first nc match
        """
        function readras(mm::Regex; missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832),kw...)
            v::Vector{String} = readdir();
            v = v[broadcast(x->endswith(x,"nc"),v)];
            file = v[(broadcast(x->occursin(mm,x),v))];
            println("reading first match of: ",file)
            f1 = first(file)
            x::Raster = read(Raster(f1,missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832);kw...))
            return(x)
        end
        

        """
        prefix::Regex; missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832),kw...)
        reads first match of regex raster
        """
        function readrasrec(prefix::Regex; missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832),kw...)
            rootdir="."
            results = []
            for (looproot, dirs, filenames) in walkdir(rootdir)
                for filename in filenames
                    #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                    if (occursin(prefix,filename)) && (endswith(filename,"nc"))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
            file = first(results)
            println("loading $file")
            x = read(Raster(file; missingval=0,crs=EPSG(25832),mappedcrs=EPSG(25832),kw...))
            describe(x)
            return(x)
        end

        function median_filter(ras::Raster)
            # Get the array and dimensions of the raster
            Z=Band(1)
            #arr = ras[:Z]
            arr = ras[Z]
            nx, ny = size(arr)
            # Create an output array with the same size and type
            out = similar(arr)
            # Loop over the pixels, excluding the borders
            for i in 2:nx-1, j in 2:ny-1
            # Get the values in the 3x3 window
            window = arr[i-1:i+1, j-1:j+1]
            # Calculate the median of the window
            out[i,j] = median(window)
            end
            # Return a new raster with the filtered array
            return rebuild(ras,out)
        end

        function cntplt(file::AbstractString)
            x=read(Raster(file,missingval=0))
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function cpl(file::Raster)
            #x=Rasters.rebuild(file;missingval=-9999)
            x=Rasters.rebuild(file;missingval=0)
            #x=x[t=1]
            Plots.contourf(x; 
            c=cgrad(:matter),
            xlabel="",
            ylabel="")
        end

        function cntplt(file::Union{Missing, Float64})
            x=file
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function cplt(file::AbstractString)
            x=read(Raster(file))
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function cpal(ext::AbstractString)
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext),file) && (!occursin("stack",file)) && endswith(file, ".nc")
            outname=replace(i,"nc"=>"jl.png");
                #println(outname)
                r=read(Raster(i,missingval=0));
                p=Plots.contourf(r;
                title=replace(basename(i),".nc"=>""), #split(outname,"/")[end], #basename(i)
                c=cgrad(:thermal),
                size=(1200, 800));
                savefig(p,outname)
                println(basename(outname)," saved!");
                end
        end
        end

        function stackplot(ext::AbstractString)
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext),file) && (occursin("stack",file)) && endswith(file, ".nc")
            outname=replace(i,"nc"=>"jl.png");
                #println(outname)
                r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
            #(i,missingval=-9999,mappedcrs=EPSG(25832))
            ee = Int(r.dims[3][end])
            rn = r[t=2:ee];    #subset till end
                p=Plots.plot(rn;
                #title=replace(basename(i),".nc"=>""), #no title cause problems
                c=cgrad(:thermal),
                size=(1200, 800));
                savefig(p,outname)
                println(basename(outname)," saved!");
                end
        end
        end

        # """
        # join([x1,y1],"+.*")
        # r"this+.*that"
        # """
        # function regand(v::Vector{String},x1::AbstractString,y1::AbstractString)
        #     needle=join([x1,y1],"+.*");
        #     z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))] 
        # return(z)
        # end
        
        # function regand(v::Vector{String},xv::Tuple{String, String})
        #     needle=join([xv[1],xv[2]],"+.*");
        #     z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
        # return(z)
        # end

        # #function regand(v::Vector{String},xv::Tuple{Symbol,Symbol})
        # function regand(v::Vector{String},xv::Vector{Symbol})
        #     needle=join(xv,"+.*");
        #     z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
        # return(z)
        # end

        # function regand(v::Vector{String},xv::Vector{String})
        #     needle=join(xv,"+.*");
        #     z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
        # return(z)
        # end

        # """
        # basic a + b regex
        # """
        # function regand(a::String, b::String)
        #     needle=join([a,b],"+.*");
        #     z=Regex(needle,"i")
        #     return(z)
        # end

        # """
        # here you can put any string to filter the Vector
        # like regand(getnames(dfs),"scn","ssr")
        # """
        # function regand(v::Vector{Any},a::String, b::String)
        #     needle=Regex(join([a,b],"+.*"),"i")
        #     z = v[(broadcast(x->occursin(needle,x),v))] 
        #     if length(z)==1
        #         return(only(z))
        #     else
        #         return(z)
        #     end
        # end


        # """
        # here you can put any regex to filter the Vector
        # like regand(getnames(dfs),r"tem")
        # """
        # function regand(v::Vector{Any},xv::Regex)
        #     z = v[(broadcast(x->occursin(xv,x),v))] 
        # return(z)
        # end


        #########tdiff 
        """
        path="."
        read non-recursivley and plots tdiff
        glob<->rglob
        stores tdiff.nc
        """
        function tdifnc(;writenetcdf=false)
            #join([x1,y1],"+.") r"Layer+.nc$"
            pot = filter(x->occursin(r"Layer.+.nc",x),glob(r"etp"))|>last|>readras
            real= filter(x->occursin(r"Layer.+.nc",x),glob(r"etr"))|>last|>readras
            td = pot-real
            nm = filter(x->occursin(r"Layer.+.nc",x),glob(r"etp"))|>last
            ti = split(nm,".")|>x->x[end-1]
            p = Plots.plot(td;
                title = "Tdiff[mm] of "*ti,
                c=cgrad(:matter))
            return p
            ##for overwriting force
            #write("tdiff.nc",td;force=true)
            if writenetcdf
                write("tdiffjl.nc",td;force=true)
                @info raw"saved to tdiffjl.nc ! "
            end
            #@info raw"for writing: write(tdiffjl.nc,RasterObject;force=true) "
            #return(td)
        end

        """
        ```
        mask_trim(raster, poly)
        Rasters.trim(Rasters.mask(raster; with=poly); pad=::Int64=10)
        ```
        """
        mask_trim(raster, poly;pad::Int64=10) = Rasters.trim(Rasters.mask(raster; with=poly); pad=pad)

        function rp(file::AbstractString; lyr=1)
            #xlyr = length(lyr)!=1 ? 1 : lyr
            ts = read(Raster(file,missingval=0))
            x = ts[t=lyr]
            #contourf(x; c=cgrad(:thermal),size=(1200, 800))
            contourf(x; 
            c=cgrad(:matter),
            xlabel="",ylabel="",
            title=name(x)
            )
        end

        function rplot(file::AbstractString;naval::Number=-9999)
            xr = read(Raster(file;crs=EPSG(25832),missingval=naval))
            Plots.plot(xr;c=cgrad(:thermal),
            xlabel="",
            ylabel="",
            size=(1200*.66, 800*.66))
        end

        function rplot(file::AbstractString, lyr::Int)
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
        end   
        
        function rplot(x::Raster;ex=1)
            """
            subset raster by last dimension.
            excludelayer and plot the rest.
            rplot(tras;ex=2)
            rplot(tras,2)
            keyword arguments by using a semicolon (;) in the parameter list. 
            """
            if x.dims|>length > 2
                @warn "subsetting to layer $ex !"
                xr = x[Dim{Rasters.name(x.dims)[end]}(Rasters.Where(x -> x >= ex))]
            else
                xr = x
            end
            Plots.plot(xr;
                c=cgrad(:thermal),
                size=(1200*.66, 800*.66),
                xlabel="",
                ylabel="",
                title=Rasters.name(xr))
        end

        """
        function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
        """
        function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
            #check it rootdir is a dir
            rootdir = isdir(rootdir) ? rootdir : throw(ArgumentError("rootdir is not a directory!"))    
            results = []
            # checks if all items in rootdir are directories
            if (any(x->isdir(joinpath(rootdir, x)), readdir(rootdir)))
                for (looproot, dirs, filenames) in walkdir(rootdir)
                    for filename in filenames
                        if pattern isa AbstractString
                            if (occursin(Regex("$(pattern)", "i"), filename))
                                push!(results, joinpath(looproot, filename)) 
                            end
                        else
                            if (occursin(pattern, filename))
                                push!(results, joinpath(looproot, filename)) 
                            end
                        end
                    end
                end
            else
                printstyled("no dirs in $rootdir !\n", color=:light_red)
                for filename in (filter(x->isfile(x), readdir(rootdir; join=true)))
                    if pattern isa AbstractString
                        if (occursin(Regex("$(pattern)", "i"), filename))
                            push!(results, filename) 
                        end
                    else
                        if (occursin(pattern, filename))
                            push!(results, filename) 
                        end
                    end
                end
            end
            return results
        end


        """
        subset raster by mask and last dimension.
        excludelayer and plot the rest.
        keyword arguments by using a semicolon (;) in the parameter list. 
        """
        function rp(x::Regex;msk::Float64=0.0001,gt=false)
            #msk=0 #msk::Float64
            if x isa Regex
                v = rglob(x)
                x = v[broadcast(x->endswith(x,"nc"),v)]|>first
                x = read(Raster(x,missingval=0;lazy=true))
            end
            
            #lastdim = x.data|>size|>last
            #lastdim = Int(x.dims[3][end])
            #rn = r[t=2:ee];    #subs
            #dimname = (Rasters.name(x.dims)[end]) #gehtnet
            #x = x[t=1]
            #x = x[Rasters.Where(k)=lastdim]
            #xs = x[lastdim=Rasters.Where(Symbol.(Rasters.name(x.dims)[end]))]
            
            xs = x[t=Int(x.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            rmsk = Rasters.mask(xs; with=zm); 
            #size=(1200, 800),
            #Plots.contourf(
            Plots.plot(
                rmsk;
                #c=cgrad(:thermal),
                c=cgrad(:matter),
                xlabel="",
                ylabel="",
                title = name(rmsk)
                )
        end

        function rp(x::Raster;msk::Float64=0.0001,gt=false)
            """
            subset raster by mask only on one dimension.
            """
            #msk=0 #msk::Float64          
            
            #xs = x[t=Int(x.dims[3][end])+1]
            
            if msk !== Float64
                msk = tryparse(Float64,string(msk))
            end
            xs = x
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            Plots.contourf(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                #size=(1200, 800),
                xlabel="",
                ylabel=""
                )
        end

        """
        rpall(file::AbstractString)
        read(Raster(file;crs=EPSG(25832),missingval=0)
        """
        function rpall(file::AbstractString)
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            #xr = read(Raster(file;crs=EPSG(25832),missingval=-9999))
            Plots.plot(xr;c=cgrad(:thermal),
            xlabel="",
            ylabel="",
            size=(1200*.8, 800*.8))
        end

        """
        path="."
        read non-recursivley
        glob<->rglob
        stores to mean.nc

        """
        function ncmean(x::Regex)
            #x=r"win"
            file = filter(file -> occursin(x,file) && endswith(file,".nc"), 
                readdir())|>first
            td::Raster = read(Raster(file,missingval=0))
            ##for overwriting force
            td = td/365
            outname = basename(file)
            m = match(r".*[.]",basename(file))
            outfile = contains(basename(file),".") ? string(m.match,"nc") : basename(file)*".html"
            # replace!(outname,".nc","mean.nc")
            # replace!(outname,"sum","")    
            write("mean.nc",td;force=true)
        end


        function facets_loop(ext::AbstractString)
            """
            like stackplot, but for interactive view
            """
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext,"i"),file) && (occursin("stack",file)) && endswith(file, ".nc")
                    @warn("subsetting first layer...")
                    r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
                    ee = Int(r.dims[3][end])
                    rn = r[t=2:ee];    #subset till end
                    p=Plots.plot(rn;
                    xlabel="",
                    ylabel="",
                    c=cgrad(:thermal),
                    size=(1200, 800));
                    return p
                end
        end
        end

        function facets(ext::AbstractString)
            """
            like stackplot, but for interactive view
            """
            grids = filter(x -> occursin(Regex(ext,"i"),x), readdir())
            grids = filter(x->endswith(x,"nc"),grids)
            if (length(grids)>0)
                file = first(grids)
                @warn("taking first match of $grids\n -> $file")
            else
                dir=pwd()
                @warn "No netcdf match found in $dir !"
                return
            end
                if isfile(file)
                    r=read(Raster(file,missingval=0,mappedcrs=EPSG(25832)));
                    if (r.dims[end]|>length == 1)
                        @warn("only one layer available...")
                        p=Plots.plot(r;
                        c=cgrad(:matter));
                        #size=(1200, 800));
                        return p
                    else
                        @info("subsetting first layer...")
                        ee = Int(r.dims[3][end])
                        rn = r[t=2:ee];    #subset till end
                        p = Plots.plot(rn;
                        xlabel="",
                        ylabel="",
                        #title=replace(basename(i),".nc"=>""), #no title cause problems
                        c=cgrad(:thermal));
                        #size=(1200, 800));
                        return p
                #savefig(p,outname)
                #println(basename(outname)," saved!");
                end
            end
        end

        """
        like stackplot, but for interactive view
        ext::Union{Regex,Raster};kw... -> Plot
        """
        function facets(ext::Union{Regex,Raster};kw...)

            if ext isa Raster
                r = ext
                if (r.dims[end]|>length == 1)
                    @warn("only one layer available...")
                    p=Plots.plot(r;
                    xlabel="",
                    ylabel="",
                    #c=cgrad(:thermal),
                    c=cgrad(:matter),
                    kw...);       
                    #size=(1200, 800));
                    return p
                else
                    @warn("subsetting first layer...")
                    ee = Int(r.dims[3][end])
                    rn = r[t=2:ee];    #subset till end
                    p=Plots.plot(rn;
                        xlabel="",
                        ylabel="",
                        c=cgrad(:thermal),
                        kw...);                       #title=replace(basename(i),".nc"=>""), #no title cause problems
                        #size=(1200, 800));
                    return p
                end
            else
                grids = filter(
                    (x -> occursin(ext,x) && endswith(x,"nc")),readdir())
                
                if (length(grids)>0)
                    file = first(grids)
                    @info("taking first match of $grids\n -> $file")
                else
                    dir=pwd()
                    @error "No netcdf match found in $dir !"
                    return
                end
                    if isfile(file)
                        r = read(Raster(file,missingval=0,mappedcrs=EPSG(25832)));
                        if (r.dims[end]|>length == 1)
                            @warn("only one layer available...")
                            p=Plots.plot(r;
                            xlabel="",
                            ylabel="",
                            #c=cgrad(:thermal),
                            c=cgrad(:matter),
                            kw...);       
                            #size=(1200, 800));
                            return p
                        else
                            @warn("subsetting first layer...")
                            ee = Int(r.dims[3][end])
                            rn = r[t=2:ee];    #subset till end
                            p=Plots.plot(rn;
                                xlabel="",
                                ylabel="",
                                c=cgrad(:thermal),
                                kw...);                      
                                #title=replace(basename(i),".nc"=>""), #no title cause problems
                                #size=(1200, 800));
                            return p
                    #savefig(p,outname)
                    #println(basename(outname)," saved!");
                    end
                end
            end
        end

        """
        Rasterstats
        """
        function descr(r::Union{Raster,String};missval=Float32(-9999),kw...)
            if r isa String
                rs = Raster(r;kw...)
                #missingval=missval
                #-9999.0f0|>typeof #Float32
                #r = replace_missing(rs, missval)
                r = replace_missing(rs, 0)
            end
            nm = name(r)
            if length(dims(r))>2
                for i in 1:last(size(r))
                    printstyled("$nm t $i \n",color=:green)
                    describe(r[:,:,i])
                end
            else
                describe(r)
            end
        end

        """
        reads with rglob
        rstats(x::Union{String,Regex})
        """
        function rstats(x::Union{String,Regex})
            v = rglob(x)
            a = []
            for i in v
                printstyled("read $i ...\n",color=:green)
                r = Raster(i)
                r = replace_missing(r, 0)
                m = mean(r)
                n = minimum(r)
                x = maximum(r)
                d = median(r)
                s = std(r)
                arr = [m, n, x, d, s]'
                df = DataFrame(arr, :auto)
                nm = ["mean", "min", "max", "median", "sd"]
                rename!(df, nm)
                df.nm .= i
                push!(a, df)
            end
            a = reduce(vcat, a)
            sort!(a, :mean, rev=true)
            return a
        end

        # function nconly(x1::Union{AbstractString,Regex})
        #     v::Vector{String} = readdir();
        #     v = v[broadcast(x->endswith(x,"nc"),v)];
        #     if x1 isa Regex
        #         z = v[(broadcast(x->occursin(x1,x),v))] 
        #     else
        #         z = v[(broadcast(x->occursin(Regex(x1,"i"),x),v))] 
        #     end
        #     return(z)
        # end

        function stats(r::Union{Raster,String,Regex};missval::Number=0,kw...)
            if r isa String
                r = Raster(r;missingval=NaN,kw...)
            end

            if r isa Regex
                v = filter(x->endswith(x,"nc"),readdir());
                r = first(v[(broadcast(x->occursin(r,x),v))])
                #r = filter(file -> occursin(r,"i"),v)|>first
                r = Raster(r;missingval=NaN,kw...)
            end
            # get the number of missing values for each band
            c = missval
            # @info "missval is set to $c"
            tbl=r.metadata|>collect|>Dict|>DataFrame
            @info "metadata:"
            pretty_table(tbl,header=uppercasefirst.(names(tbl)), backend = Val(:text))
            
            r = replace_missing(r, c)
            m = mean(r|>skipmissing) # get the mean for each band
            n = minimum(r|>skipmissing) # get the minimum for each band
            x = maximum(r|>skipmissing) # get the maximum for each band
            d = median(r|>skipmissing) # get the median for each band
            s = std(r|>skipmissing) # get the standard deviation for each band
            ct = r|>skipmissing|>collect|>length
            rs = size(r)
            c = reduce(*, rs) - ct #count(ismissing, r) # get the number of missing values for each band
            rname = string(name(r))
            
            arr=[m,n,x,d,s,ct,c]'
            #println("$nm\n",arr)
            df = DataFrame(arr,:auto)
            nm=["mean", "min", "max", "median", "sd", 
                "count",
                "missval"]
            rename!(df,nm)
            #ternary if else
            #df.name .= isnothing(rname) ? tbl.filepath[1] : rname
            
            df.name .= isempty(rname) ? normpath(tbl.filepath[1]) : rname

            return(df)
        end

        function rpr(x::String)
            """
            3D plot with Rasters
            """
            #x="D:/temp/saale/output/thulba/v2/All_HydrologicResponseUnits.nc"
            ga = Rasters.read(Rasters.Raster(x))
            values = ga.data[:,:,1]
            #coords = (1:length(lookup(ga, X)), 1:length(lookup(ga, Y)))
            s=size(ga)
            coords = (1:s[2], 1:s[1]) #swiched on netdf!
            ti=basename(x)
            #p1=
            Plots.surface(        coords[1], coords[2], 
                values    ,   xlabel="x",
                ylabel="y",zlabel="rastervalue",title=ti)
        end

        function rpr(x::Raster)
            """
            3D plot with Rasters
            x .> 0 
            """
            msk = float(0.001)
            xs = (length(x.dims) > 2) ? x[t=Int(x.dims[3][end])+1]  : x
            zm = (xs .> msk)
            # #zm = (gt) ? (xs .> msk) : (xs .< msk)
            # x = Rasters.mask(xs; with=zm); 
            # values = x.data[:,:,1] .> 0
            # s=size(x)
            # coords = (1:s[2], 1:s[1]) #swiched on netdf!
            ti=Rasters.name(x)
            M=Rasters.mask(xs; with=zm);
            #M=Rasters.rebuild(M;missingval=msk)
            #M=Rasters.rebuild(M;missingval=msk)|>skipmissing
            Plots.surface(
                M,
                #legend=:outerbottomleft,
                #legend=:bottom,
                #legend= false,
                xlabel=" ",
                ylabel=" ",
                zlabel=" ",
                title=ti,
                camera = (-20, 75)        
                )
            # Plots.surface(coords[1], coords[2], 
            #     values    ,   
            #     xlabel="x",
            #     ylabel="y",zlabel=ti,
            #     #title=ti,
            #     camera = (-20, 75))
        end
        
        function filterplot(regex::AbstractString,ncs::Vector{Raster})
            "selects first match and plots..."
            r = ncs[map(n->occursin(Regex(regex,"i"),n),
            map(x->string.(name(x)),ncs)
            )] |> first
            plot(trim(r),xlab="",ylab="",title=name(r))
            cz=r.metadata.val|>collect|>last|>last
            Plots.annotate!(:bottomright,
            text("cellsize: $cz", 8, :black, :right))
        end

        function filterplot!(regex::AbstractString,ncs::Vector{Raster})
            "selects first match and add to plot..."
            r = ncs[map(n->occursin(Regex(regex,"i"),n),
            map(x->string.(name(x)),ncs)
            )] |> first
            plot!(trim(r),xlab="",ylab="",title=name(r))
        end

        function rhist(r::Raster)
            M=r.data[:,:,1]
            Plots.histogram(skipmissing(M),title=r.name,legend=false)
        end

        function rp3(x::Raster)
            """
            3D plot with Rasters from Raster.
            """
            #t = x|>size
            #coords = (1:t[1], 1:t[2])
            ti= try 
                basename(x);
            catch 
                @warn "no basename available, trying to parse name..."; 
                ti=Rasters.name(x); 
            end
            x.dims    
            Plots.surface(x[:, :, 1],
            xlabel="x",ylabel="y",zlabel="value",title=ti)
        end

        """
        subset raster by mask and last dimension.
        excludelayer and plot the rest.
        keyword arguments by using a semicolon (;) in the parameter list. 
        rpm(r"rad";msk=0.0,gt=true) ->plots
        """
        function rpmcf(x::Regex;msk::Float64,gt=false)
            #msk=0 #msk::Float64
            v = glob(x)
            r = v[broadcast(x->endswith(x,"nc"),v)]|>first
            #r = nconly(x)|>first

            #first("∀ϵ≠0: ϵ²>0", 1)

            r = read(Raster(r,missingval=0;lazy=true))
            xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.contourf(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end

        """
        subset raster by mask and last dimension.
        excludelayer and plot the rest.
        keyword arguments by using a semicolon (;) in the parameter list. 
        rpm(r"rad";msk=0.0,gt=true) ->plots
        """
        function rpm(x::Regex;msk::Float64,gt=false)
            #msk=0 #msk::Float64
            v = glob(x)
            r = v[broadcast(x->endswith(x,"nc"),v)]|>first
            #r = nconly(x)|>first

            #first("∀ϵ≠0: ϵ²>0", 1)

            r = read(Raster(r,missingval=0;lazy=true))
            xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.plot(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end


        """
        subset raster by mask and plot.
        no subset on last dim!
        keyword arguments by using a semicolon (;) in the parameter list. 
        maskplot(r;msk=0.0) ->plots
        """
        function maskplot(xs::Raster;msk::Float64,gt=true)
            #msk=0 #msk::Float64
            #xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.plot(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:matter),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end
        
        function readras2(file::AbstractString)
            x=read(Raster(file,missingval=-9999.000000)) #read all in RAM
            return(x)
        end

        function fzplot(dirpath=".";kwargs...)
                
            files = readdir(dirpath;join=true)
            rs = try
                filter(file -> occursin(Regex(".fzs\$"),file),files)|>first

            catch
                @error "No .fzs file found in $dirpath !"
                return nothing  # 
            end 
            r = rs
            r = Raster(r;missingval=-9999)
            r = r ./ 3600|>Rasters.trim
            # msk = float(0.001)
            # #msk = Float32(0)
            # zm = (r .> msk)
            # r = Rasters.mask(r; with=zm); 
            r = Rasters.rebuild(r,missingval=minimum(r));
            plot(r,c=:blues,
            #title="flowtime in [h]",
            title="",
            colorbar_title = "Fliesszeitsummen in [h]",
            xaxis="",yaxis="";kwargs...)
        end

        """
        AG only.
        fzplot2(x::AbstractString)
        """
        function fzplot2(x::AbstractString)
            r = ArchGDAL.readraster(x)
            r = ArchGDAL.getband(r,1)
            dx = r ./ 3600
            dx = permutedims(dx, (2, 1))
            dx = reverse(dx, dims=1)
            Plots.contour(dx,
            colorbar_title = "flowtime in [h]",
            title="",xaxis="",yaxis="")
            #colorbar_title = "",
            # xm = .!all(ismissing, dx, dims=2)
            # ym = .!all(ismissing, dx, dims=1)
            # A = dx[xm[:], ym[:]]
            # Plots.contour(A,title="flowtime in [h]",xaxis="",yaxis="")
            
        end

        """
        reads ascii grid and replaces -9999 \n returns Raster
        """
        function agread(x::String)
            dataset = ArchGDAL.read(x)
            # Get the first band
            band = ArchGDAL.getband(dataset, 1)
            # A = permutedims(band, (2, 1))
            # A = A.a
            # A = reverse(A, dims=1)
            A = reverse(band, dims=2)
            B = replace(A, -9999.0 => missing)
            xm = .!all(ismissing, B, dims=2)
            ym = .!all(ismissing, B, dims=1)
            A = B[xm[:], ym[:]]
            k = Raster(A,(X,Y))
            println(extrema(k))
            return k
        end

        """
        `s::AbstractString;
        step=30,lyr=1,msk=0.001,
                umask::Number=10^6,
                anns=true,
                cmap::Symbol=:greys, #:thermal,:lightrainbow
                roundto=false,
                noturn=false,kw...`\n
        kw passed to heatmap
        plots heatmap from raster with ArchGDAL\n
        roundto can also be an integer...
        """
        function agheat(s::Union{AbstractString,Raster};
                step=30,
                lyr=1,
                msk=0.001,
                anns=true,
                umask::Number=10^6,
                cmap::Symbol=:greys, #:thermal,:lightrainbow
                roundto=false,
                noturn=false,kw...)
            if s isa Raster
                println(s.metadata)
                dx = try 
                    convert(Matrix{Float32},s)
                catch
                    @error "Raster could not be converted to Matrix{Float32}!"
                end
                csx = try 
                        Rasters.cellsize(dx)|>first
                    catch
                        @error "could not be retrieve cellsize!"
                        200 #return 200
                    end
            else
                if !isfile(s)
                    error("file not found!")
                end
                        
                r = ArchGDAL.readraster(s)
                dx = ArchGDAL.getband(r,lyr)
                x = ArchGDAL.getgeotransform(r)[2]
                y = ArchGDAL.getgeotransform(r)[6]
                @info "gdal cellsize info: x: $x y: $y"
                #cellsize (x-direction) / 8
                csx = ArchGDAL.getgeotransform(r)[2] ./ 8
            end
                step = step
                if step < csx
                    (round(csx))
                end
                @info "annotation step is $step"

                if umask==10^6
                    umask = maximum(dx)
                end
            
            println("MIN:",minimum(r),"\nMAX:",maximum(r))
            
            bitmat = (dx .> msk) .& (dx .<= umask)
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            #dx_output .= NaN #fill!
            fill!(dx_output,NaN32)
            # dx_output = Matrix{Int8}(undef, size(dx, 1),size(dx,2))
            # dx_output = typeof(dx)(undef, size(dx, 1),size(dx,2))
            dx_output[bitmat] .= dx_filtered[bitmat]
            dx = dx_output
            println("masking value range to $msk and $umask ...")
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF!")
                dx = permutedims(dx, (2, 1))
                dx = reverse(dx, dims=1)
            elseif noturn
                @info("processing $s without flip...")
                #dx = Rasters.Raster(s)
            else
                @info("processing $s ...")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end   ##:lightrainbow
            heatmap_plot = heatmap(dx, 
                            c=cmap,
                            title=basename(s), 
                            xaxis="", yaxis="",kw...);
            
            if anns
                for i in 1:step:size(dx, 1)
                    for j in 1:step:size(dx, 2)
                        if roundto==false
                            value = dx[i, j]
                            if !isnan(value)
                                value = round(Int64, value)
                            end
                        else
                            value = round(dx[i, j]; digits=roundto)
                        end
                        color = isnan(value) ? :transparent : :black
                        annotate!(j, i, Plots.text(
                            #replace(string(value),r".0$"=>""), 
                            string(value),
                            7, color, :center, 
                            halign=:center, rotation=-35.0))
                    end
                end
            end
                # Show the plot
            #display(heatmap_plot)
            return heatmap_plot
        end

        """
        AG only.
        import ArchGDAL as AG
        plots raster with ArchGDAL
        """
        function agcont(s::AbstractString;lyr=1,msk=0.001)
            if !isfile(s)
                error("file not found!")
            end
            r = ArchGDAL.readraster(s)
            dx = ArchGDAL.getband(r,lyr)
            # Filter the dx matrix using bitmat
            bitmat = dx .> msk
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]
            dx = dx_output

            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF!")
                dx = permutedims(dx, (2, 1))
                dx = reverse(dx, dims=1)
            else
                @info("processing $s ...")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end   
            #println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            println("MIN:",minimum(r),"\nMAX:",maximum(r))
            #println(extrema(dx))
            #dx = reverse(dx, dims=1)|>x->reverse(x, dims=2)
            #dx = reverse(dx, dims=:) #alldims
            #dx = reverse(dx, dims=2) #alldims

            #Plots.contour(dx,c=:matter,title=basename(s),xaxis="",yaxis="")
            Plots.contourf(dx, fill=true, levels=10,
                #c=:matter,
                title=basename(s),xaxis="",yaxis="")

        end

        """
        AG only.
        import ArchGDAL as AG
        agcont2(f;lyr=5,msk=-1)
        agcont2(f;lyr=2,msk=-1,step=50)
        s::AbstractString;step=30,lyr=1,msk=0.001,roundto=false
        """
        function agcont2(s::AbstractString;step=30,lyr=1,msk=0.001,roundto=false)
            if !isfile(s)
                error("file not found!")
            end
            r = ArchGDAL.readraster(s)
            dx = ArchGDAL.getband(r,lyr)
            x = ArchGDAL.getgeotransform(r)[2]
            y = ArchGDAL.getgeotransform(r)[6]
            @info "gdal cellsize info: x: $x y: $y"
            #cellsize (x-direction) / 2
            step = step
            csx = ArchGDAL.getgeotransform(r)[2] ./ 8
            if step < csx
                step = Int(round(csx))
            end
            @info "annotation step is $step"
            println("MIN:",minimum(r),"\nMAX:",maximum(r))       
            # Filter the dx matrix using bitmat
            bitmat = dx .> msk
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]
            dx = dx_output

            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF!")
                dx = permutedims(dx, (2, 1))
                dx = reverse(dx, dims=1)
            else
                @info("processing $s ...")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end   
            #println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            
            
            c_plot = Plots.contour(dx,
                    #c=:matter,
                        title=basename(s),xaxis="",yaxis="");
            #step = step
            for i in 1:step:size(dx, 1)
                for j in 1:step:size(dx, 2)
                    if roundto==false
                        value = dx[i, j]
                        if !isnan(value)
                            value = round(Int64, value)
                        end
                    else
                        value = round(dx[i, j]; digits=roundto)
                    end
                    color = isnan(value) ? :transparent : :black
                    annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                        halign=:center, rotation=35.0))
                end
            end

            # for i in 1:step:size(dx, 1)
            #     for j in 1:step:size(dx, 2)
            #         value = round(dx[i, j]; digits=2)
            #         color = isnan(value) ? :transparent : :black
            #         rotation_char = isnan(value) ? "" : "↺"
            #         annotate!(j, i, text(string(value) * rotation_char, 5, color, :center))
            #     end
            # end

            # Show the plot
            display(c_plot)
        end

        """
        s::AbstractString;
        camera=(60, 30, 30),
        lyr=1,msk=0.1,umask=10^6
        """
        function agsurf(s::AbstractString;camera=(60, 30, 30),lyr=1,msk=0.1,umask=10^6)
            if !isfile(s)
                error("file not found!")
            end
            r = ArchGDAL.readraster(s)
            dx = ArchGDAL.getband(r,lyr)
            if umask==10^6
                umask = maximum(dx)
            end
            println("MIN:",minimum(r),"\nMAX:",maximum(r))
            println("masking value range to $msk and $umask ...")
            bitmat = (dx .> msk) .& (dx .<= umask)
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]
            dx = dx_output
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF!")
                dx = permutedims(dx, (2, 1))
                dx = reverse(dx, dims=1)
            else
                @info("processing $s ...")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end 
            
            c_plot = Plots.surface(dx,
                    c=:matter,
                    camera = camera,
                        title=basename(s),
                        xaxis="",yaxis="");
            return (c_plot)
        end

        """
        AG only.
        plots heatmap from raster with ArchGDAL
        upper and lower bound
        """
        function agmask(s::AbstractString;step=50,lyr=1,lower=0,upper=1000)
            if !isfile(s)
                error("file not found!")
            end
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF !")
            end
            r = ArchGDAL.readraster(s)
            dx = ArchGDAL.getband(r,lyr)
            #ArchGDAL.getnodatavalue(dx)
            
            if endswith(s,".nc")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end
            #findmin(dx)
            if (minimum(dx) <= -9999.0) # || (maximum(dx) >= 9999.0)
                println("extrema: ",join(extrema(dx)," <-> "))
                dmin=minimum(dx)
                @info "
                $dmin set to NaN!
                lower bound: $lower
                upper bound: $upper
                "
                dx = replace(dx, minimum(dx)=>NaN)
                #replace!(dx, minimum(dx)=>missing) #err
                #replace!(dx, minimum(dx)=>NaN)
            end
                    
            #println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            
            bitmat = (dx .> lower) .& (dx .< upper)
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]
        
            if allequal(dx_output)
                @error("all values are equal!")
                return
            end
        
            dx = dx_output
            heatmap_plot = heatmap(dx, c=:matter,
                            title=basename(s), 
                            xaxis="", yaxis="");
            step = step
            for i in 1:step:size(dx, 1)
                for j in 1:step:size(dx, 2)
                    value = round(dx[i, j]; digits=2)
                    color = isnan(value) ? :transparent : :black
                    annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                        halign=:center, rotation=-35.0))
                end
            end
            # Show the plot
            display(heatmap_plot)
        end

        """
        ezplot(dirpath::String=".";kwargs...)
        """
        function ezplot(dirpath::String=".";kwargs...)
            #cd(dirpath)
            
            files = readdir(dirpath;join=true)
            rs = try
                filter(file -> occursin(Regex(".ezg\$"),file),files)|>first
            catch
                @error "No .ezg file found in $dirpath !"
                return nothing  # 
            end 
            
            lstr = filter(file -> occursin(Regex(".lnk\$"),file),files)|>first
            ezg = Raster(rs;missingval=-9999)   |>Rasters.trim
            lnk = Raster(lstr;missingval=-9999) #|>Rasters.trim
            p = Plots.plot(ezg,
                #c=:lightrainbow,
                c = :matter,
                #alpha = 0.7,
                xaxis="",yaxis="");
            Plots.plot!(lnk,c=:white,
                title=dirname(dirpath),
                xaxis="",yaxis="";kwargs...)
        
            return p
            
            # r = ArchGDAL.readraster(rs)
            # dx = ArchGDAL.getband(r,1)
            # dp = AG.polygonize(dx)
        end

        function ncdf(str::String)
            rr = Raster(str)
            @info "this extracts values from the midpoint of the raster.
            rename!(df,1=>Rasters._maybename(rr),2=>...)
            last dimension is time..."

            #nm = Rasters.name(rr)
            nm = Rasters._maybename(rr)
        
            if dims(rr)[1]|>name != :X
                new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                # Create a new raster with the new dimensions and the same data
                rr = try
                    Raster(rr.data, dims=new_dims)
                catch
                    @error "var conversion failed!"
                    return
                end 
            end
        
            ti = Rasters.lookup(rr,3)
            x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
            df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
            df = hcat(df, parent(ti), makeunique=true)
            #rename!(df,1=>Rasters._maybename(rr),2=>"date")
            rename!(df,1=>Symbol(nm),2=>:date)
            if typeof(df.date) == Vector{CFTime.DateTimeNoLeap}
                df.date = map(x->tryparse(DateTime,string(x)),df.date)
            end
            DataFrames.metadata!(df, "filename", str, style=:note);        
            return df
        end

        function ncdf(rr::Raster)
            @info "this extracts values from the midpoint of the raster.
            rename!(df,1=>Rasters._maybename(rr),2=>...)
            last dimension is time..."

            #nm = Rasters.name(rr)
            nm = Rasters._maybename(rr)
        
            if dims(rr)[1]|>name != :X
                new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                # Create a new raster with the new dimensions and the same data
                rr = try
                    Raster(rr.data, dims=new_dims)
                catch
                    @error "var conversion failed!"
                    return
                end 
            end
        
            ti = Rasters.lookup(rr,3)
            x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
            df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
            df = hcat(df, parent(ti), makeunique=true)
            #rename!(df,1=>Rasters._maybename(rr),2=>"date")
            rename!(df,1=>Symbol(nm),2=>:date)
            DataFrames.metadata!(df, "filename", nm, style=:note);
            
            #if df.date[1]|>typeof == CFTime.DateTimeNoLeap
            #parse(DateTime,string(df.date[1]))
            if df.date[1] != Date
                df.date = map(x->tryparse(DateTime,string(x)),df.date)
            end
            

            return df
        end

        """
        3D plot with Rasters
        rst.surf(rr[Ti=12];msk=3,cam=(12,45),c=:lightrainbow)
        
        """
        function surf(x::Raster;msk=0.001,c=:jet,cam=(-20, 40))
            #0.001 isa Float64
            #ti=Rasters._maybename(x)
            #ti=Rasters.name(x)

            if msk !== Float64
                msk = tryparse(Float64,string(msk))
            end       
            
            #xs = (length(x.dims) > 2) ? x[t=Int(x.dims[3][end])+1]  : x
            xs = x
            zm = (xs .>= msk)
                
            M = Rasters.mask(xs; with=zm);
            
            Plots.surface(
                M,
                c = c,
                #title = ti,
                legend= false,
                xlabel=" ",
                ylabel=" ",
                #zlabel=" ",
                camera = cam        
                )
            #     camera = (-20, 75))
        end

        """
        use of Rasters.lookup
        @info "this extracts values from the midpoint of the raster.
        rename!(df,1=>Rasters._maybename(rr),2=>...)
        last dimension is time..."

        ser = RasterSeries(nconly("thaw")[2:end], :t;missingval=0.0)
        ok = ncdf(ser)
        ok = ok[!,Not(:date)]
        plot(Matrix(ok)'|>reverse,names(ok)|>reverse,xrotation = 45,legend=false)

        """
        function ncdf(xs::RasterSeries)
            #for rr in xs;println(name(rr));end
            dfs = []
            for rr in xs;
                @info "generating" (name(rr))
                nm = Rasters._maybename(rr)
                if dims(rr)[1]|>name != :X
                    new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                    # Create a new raster with the new dimensions and the same data
                    rr = try
                        Raster(rr.data, dims=new_dims)
                    catch
                        @error "var conversion failed!"
                        return
                    end 
                end
            
                ti = Rasters.lookup(rr,3)
                x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
                df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
                df = hcat(df, parent(ti), makeunique=true)
                #rename!(df,1=>Rasters._maybename(rr),2=>"date")
                rename!(df,1=>Symbol(nm),2=>:date)
                #DataFrames.metadata!(df, "filename", str, style=:note); 
                push!(dfs,df)
            
                end
                
            odf = reduce((left, right) -> 
                innerjoin(left, right, on = :date,makeunique=true), 
                dfs)

            return hcat(odf[!,Not(Cols(r"date"))],odf[:,Cols(r"date")])
        end

        """
        filterplot for rasters regex
        """
        function fp(x::Regex;)
            x = glob(x)
            filter!(z->endswith(z,".nc"),x)
            x = first(x)
            r = read(Raster(x,missingval=0;lazy=true))
            println(descr(r))
            Plots.plot(
                r;
                #c=cgrad(:thermal),
                c=cgrad(:matter),
                xlabel="",
                ylabel="",
                title = name(r)
                )
        end

        """
        reprojects a wasim NetCDF 
        using Rasters and ArchGDAL 
        to a new coordinate system
        rebuilds the raster using misval argument
        """
        function project_o(x::Union{String,Raster,String};
            misval::Number=0,
            #layer = 1,
            s_srs="EPSG:25832",
            t_srs="EPSG:4326",method="bilinear")
            if x isa Regex
                #x = nconly(x)|>first
                try
                    v::Vector{String} = readdir();
                    v = v[broadcast(f->endswith(f,"nc"),v)];
                    z = v[(broadcast(f->occursin(x,f),v))];
                    file = z|>first    
                    r = Raster(file,missingval=-9999;crs=s_srs)
                catch
                    throw("no file found!")
                end
            elseif x isa String
                try
                r = Raster(x,missingval=-9999;crs=s_srs)
                catch
                    throw("no file found!")
                end           
            elseif x isa Raster
                file = x
                #works, because self
                r = Raster(file,missingval=-9999;crs=s_srs)
            else
                throw(x)
                @error("x must be a Raster or String")
            end
            
            #r = r[t=layer]
            flags = Dict(
                "s_srs"=>s_srs,
                "t_srs"=>t_srs,
                "r"=>method)
            rs = Rasters.warp(r,flags)
            rout = replace_missing(rs, missval)
            #rs = Rasters.trim(Rasters.rebuild(rs,missingval=misval))
            return rout
        end

        """
        project(x::Union{Regex,Raster,String}
        reprojects a wasim NetCDF 
        using Rasters and ArchGDAL to a new coordinate system
        rebuilds the raster using misval argument
        ```
        args:
        misval::Number=0,
        src="EPSG:25832"
        dst="EPSG:4326",
        imethod="bilinear"
        ```
        kwargs... passed to Rasters.rebuild
        note that ArchGDAL.wkbPolygon is reporjected inplace!
        see also @doc Rasters.resample
        """
        function project(x::Union{Regex,Raster,String,ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}};
            misval::Number=0,
            src="EPSG:25832",
            dst="EPSG:4326",
            imethod="bilinear",kwargs...)
            if x isa Regex
                #x = nconly(x)|>first
                try
                    z = filter(f -> occursin(x,f), readdir())
                    file = z|>first    
                    A = Raster(file)
                catch
                    throw("no file found!")
                end
            elseif x isa String
                try
                    A = Raster(x)
                catch
                    throw("no file found!")
                end           
            elseif x isa Raster
                file = x
                #works, because self
                A = Raster(file) #;crs=s_srs
            elseif x isa ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}
                #z = deepcopy(x) #does not work
                reprojected_polygon = ArchGDAL.reproject(
                        x,
                        EPSG(src),
                        EPSG(dst)
                    )
                return reprojected_polygon
            else
                throw(x)
                @error("x must be a Raster or String")
            end       

            flags = Dict(
                "s_srs"=>src,
                "t_srs"=>dst,
                "r"=>imethod)
            rs = Rasters.warp(A,flags)
            rs = Rasters.trim(Rasters.rebuild(rs,missingval=misval;kwargs...))
            return rs
        end

        """
        reprojects a wasim NetCDF from EPSG:25832 to EPSG:4326
        masks 0 
        file = nconly(x)|>last
        r = Raster(file;kwargs...)
        """
        function rlpro(x::String;kwargs...)
            flags = Dict(
                "s_srs"=>"EPSG:25832",
                "t_srs"=>"EPSG:4326",
                "r"=>"bilinear")
            file = nconly(x)|>last
            r = Raster(file;kwargs...)
            #return Rasters.warp(r,flags)
            rs = Rasters.warp(r,flags)
            return Rasters.trim(Rasters.rebuild(rs,missingval=0))
        end

        """
        typeof(dx) Matrix{Float32}
        
        first(t) * last(t)
        """
        function ncell(dx::Union{Raster,Matrix{Any}})
            t = size(dx)
            return (first(t) * last(t))
        end
        
        #ncell(t::Tuple{Int64, Int64}) = return (first(t) * last(t))
        #ncell(t::Tuple{Int64, Int64}) = return (first(t) * last(t))
        #size(xr)

        """
        zonal stats for a raster
        jlzonal(shapefilepath, rasterpath;agg=mean,kwargs...)
        
        """
        function jlzonal(shapefilepath, rasterpath;agg=mean,kwargs...)
            shp = Shapefile.Table(shapefilepath)
            #shp = GeoDataFrames.read(shapefilepath)
            raster = Raster(rasterpath)
            #boundary: for polygons, include pixels where the :center is inside the polygon, where the line :touches the
            #pixel, or that are completely :inside inside the polygon. The default is :center.
            df = Float64.(zonal(agg, raster; of=shp,kwargs...))        
            df = DataFrame(geometry=shp.geometry,val=df)
            rename!(df,:val => string(agg))
            #df = Float64.(zonal(agg, raster; of=shp.geometry, boundary=:center))
            return df
        end

        """
        crops outline of polygon from raster
        r::Raster,geom::Vector{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}}
        same would be: Rasters.mask(r;with=gd.geometry|>first)
        """
        function rmask(r::Raster,geom::Vector{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}})
            ks = geom|>first
            x = [first(pt) for pt in GeoInterface.coordinates(ks)|>first]
            y = [last(pt) for pt in GeoInterface.coordinates(ks)|>first]
            poly_points = [(x, y) for (x, y) in zip(x, y)]
            # Create a single linear ring
            linear_ring = [poly_points]
            # Create a polygon
            polygon = GeoInterface.Polygon(linear_ring)
            # crop it
            xrm = Rasters.mask(r;with=polygon)
            return xrm
        end

        """
        rcont(x::Union{String,Regex,Raster};t=1,lw=0.25,kw...)
        """
        function rcont(x::Union{String,Regex,Raster};t=1,lw=0.25,kw...)
            #Plots.colorbar_style(orientation=:h)
            if isa(x,Regex)
                fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]
                r = Raster(fn;kw...) #missingval,
            elseif isa(x,Raster)
                r = x
                fn = string.(name(ds))
            else
                fn = x
                r = Raster(fn,missingval=0)
            end
        
        
            z = r[t=t]
            contourf(z;c=:matter,cbar=false)
            contour!(z,color=:black,legend=false,
                title="",xlab="",ylab="",linewidth=lw,
                clims=extrema(z.data|>skipmissing), 
                cbar=true,
                colorbar_title = "",
                margins=5mm,
                kw...)
            
        end

        """
        change rastervalues to newval
        """
        function process_rasters(input_folder::AbstractString, inras::AbstractString, outras::AbstractString, newval::Real)
            cd(input_folder)
            
            if isfile(outras)
                @info "File $outras exists, overwriting ..."
            end
            
            begin 
                raster = Raster(inras)
                raster.data .= round(newval, digits=4)
                write("tmp.asc", raster; force=true)
                run(`wsl ./rwperl.sh tmp.asc $outras`)
                rm("tmp.asc")
            end
            
            
        end

        
        """
        facetsplot of RasterSeries by timekey= :t
        """
        function serplot(ser::Union{RasterSeries,Regex,AbstractString}; 
            skipfirst::Bool=false,
            misval::Number=-9999,   #0.0,
            timekey::Symbol=:t,
            c=:matter, title="", xlabel="", ylabel="")
            if !(ser isa RasterSeries) 
                matches = try 
                        nconly(ser)
                    catch
                        @error "no matches of $ser found!"
                        return
                    end
                if skipfirst
                    matches = matches[2:end]
                end
                ser = RasterSeries(matches,
                    timekey;
                    missingval=misval)
            end
            @info "missingval is $misval ..."
            Plots.plot(ser; 
            c=cgrad(c), title=title, xlabel=xlabel, ylabel=ylabel)
        end


        """
        r = rst.readras(file)
        xmin = dims(r, (X, Y))[1]|>first
        xmax = dims(r, (X, Y))[1]|>last
        ymin = dims(r, (X, Y))[2]|>first
        ymax = dims(r, (X, Y))[2]|>last
        crop_poly = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
        crop_poly = ArchGDAL.createpolygon(crop_poly)
        xdf = crop_geodf(sm, crop_poly)
        edf = reduce(GeoDataFrames.union,xdf.geometry)
        """
        function crop_geodf(sm::DataFrame, crop_poly::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon})
            # Perform intersection
            np = [gdf.intersection(geom, crop_poly) for geom in sm.geometry]
            
            # Filter out empty geometries
            npw = filter(x -> !ArchGDAL.isempty(x), np)
            
            # Get the indices of the non-empty geometries
            indices = findall(x -> !ArchGDAL.isempty(x), np)
            
            # Create a new GeoDataFrame with the non-empty intersections
            new_gdf = DataFrame(geometry = npw)
            
            # Copy the attributes from the original GeoDataFrame
            for col in names(sm)
                if col != :geometry
                    new_gdf[!, col] = sm[!, col][indices]
                end
            end
            
            return new_gdf
        end   

        """
        Rasters
        """
        function remocont(s::AbstractString;lyr=1,msk=0.001)
            if !isfile(s)
                error("file not found!")
            end
            # r = ArchGDAL.readraster(s)
            # dx = ArchGDAL.getband(r,lyr)
            # # Filter the dx matrix using bitmat
            # bitmat = dx .> msk
            # dx_filtered = dx .* bitmat
            # dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            # dx_output .= NaN
            # dx_output[bitmat] .= dx_filtered[bitmat]
            # dx = dx_output        
            dx = Raster(s)
            @info("processing $s ...")
            dx = reverse(dx, dims=1)
            dx = reverse(dx, dims=2)
            println("MIN:",minimum(r),"\nMAX:",maximum(r))
            Plots.contourf(dx, 
                fill=true, 
                levels=10,
                title=basename(s),xaxis="",yaxis="")
        end

        """
        reads stationdata
        """
        function stp(fn::String;proj=false,dst=ProjString("+proj=longlat +datum=WGS84 +no_defs"))
            fl = CSV.read(fn,DataFrame;limit=4)
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            if first(yc) isa String31
                yc = tryparse.(Float64,yc)
            end
            if first(xc) isa String31
                xc = tryparse.(Float64,xc)
            end

            pts = ArchGDAL.IGeometry[]
            if proj #projects to 4326
                for i in 1:length(xc)
                    pt = ArchGDAL.createpoint([xc[i],yc[i]])
                    pt = ArchGDAL.reproject(pt,EPSG(25832),
                    dst)
                    push!(pts,pt)
                end
            else
                for i in 1:length(xc)
                    pt = ArchGDAL.createpoint([xc[i],yc[i]])
                    #pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                    push!(pts,pt)
                end
            end
            nd = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            return nd
        end

        stread=stp

        """
        reads from stationdata, reprojects and plots
        """
        function stplot(fn::String)
            fl = CSV.read(fn,DataFrame;limit=4)
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            p = plot(od.geometry);
            for (i, pt) in enumerate(od.geometry)
                #x = od.xc[i]
                x = ArchGDAL.getx(od.geometry[i], 0)
                #y = od.yc[i]
                y = ArchGDAL.gety(od.geometry[i], 0)
                name = od.name[i]
                annotate!(x, y, text(name, 8, :black, :bottom, :left))
            end
            return plot!(p)
        end

        """
        adds points to an existing plot from a specfile
        """
        function addplot(specfile::Union{String,DataFrame};reproject=false)
            if specfile isa String
                fl = CSV.read(specfile,DataFrame;limit=4)
                xc = fl[2,5:end]|>collect
                yc = fl[3,5:end]|>collect
                pts = ArchGDAL.IGeometry[]
                for i in 1:length(xc)
                    pt = ArchGDAL.createpoint([xc[i],yc[i]])
                    if reproject
                        pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                    end
                    push!(pts,pt)
                end
                od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            else
                od = specfile
            end
            
            Plots.plot!(od.geometry;c=:black,shape = :x) #:utriangle)
            for (i,pt) in enumerate(od.geometry)
                x = ArchGDAL.getx(od.geometry[i], 0)
                y = ArchGDAL.gety(od.geometry[i], 0)
                name = od.name[i]
                #Plots.text(name,:consolas, 8, :black, :bottom, :left)
                #Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
                Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
            end
            Plots.plot!()
        end

        function agjson(jsonfile::AbstractString)
            """
            reads json and transforms points...
            but AG transformation is wrong :(


            ERROR: type UnionAll has no field read
            works only in WSL
            """
            
            geojson_file=jsonfile
            jsonbytes = read(geojson_file) # read the geojson file as bytes
            fc = GeoJSON.read(jsonbytes)
            pts=[]
            for geom in fc.geometry
                xc = [(x) for x in geom]|>first|>first
                yc = [(x) for x in geom]|>first|>last
                pt = ArchGDAL.createpoint(xc,yc)
                pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            # df = [] ##geht auch
            # for x in pts
            #     x1 = AG.getpoint(x,0)
            #     tmp = DataFrame(x=[x1[1]], y=[x1[2]])
            #     push!(df,tmp)
            # end
            # df = reduce(vcat, df) 
            #Plots.plot(df.x, df.y, seriestype=:scatter)
            df = DataFrame( 
                "x" => [AG.getpoint(x,0)[1] for x in pts],
                "y" => [AG.getpoint(x,0)[2] for x in pts] )
            return(df)
        end

        """
        wrapper for ArchGDAL.unsafe_gdaldem
            operation:
            one of "hillshade", "slope", 
            "aspect", "color-relief", "TRI", "TPI",
            "Roughness"
        """
        function gdalop(input_file, output_file, operation::String="TPI")
            dataset = ArchGDAL.read(input_file)
            ArchGDAL.unsafe_gdaldem(dataset, operation, ["-of", "AAIGrid"], dest=output_file)
        end

        """
        function to subset a raster by values
        """
        function rsubset(r::Union{String,Raster},lb::Real,ub::Real;wasim=true,lyr::Int64=1)
            if isa(r,String)
                r = Raster(r;missingval=-9999,crs=EPSG(25832))
                r = r[:,:,lyr] #first layer
                #r = replace_missing(r, -9999)
            end
            A = Matrix(r)    
            if wasim
                #A = reverse(A, dims=2)  #if wasim netcdf
                A = reverse(A, dims=1)
                # A = replace_missing(A, -9999)
            end
            A = replace(A, missing => -9999)    
            A = map(x -> x <= lb ? -9999 : x, A)
            A = map(x -> x >= ub ? -9999 : x, A)
            B = replace(A, -9999 => missing)
            # Find rows and columns where all values are missing
            xm = .!all(ismissing, B, dims=2)
            ym = .!all(ismissing, B, dims=1)
            # Remove rows and columns where all values are missing
            A = B[xm[:], ym[:]]
            return A
        end

        """
        build a GeoDataFrame from a shapefile
        """
        function jlzonal2(shapefilepath, rasterpath;agg=mean,kwargs...)
            shp = Shapefile.Table(shapefilepath)
            raster = Raster(rasterpath)
            df = Float64.(zonal(agg, raster; of=shp,kwargs...))        
            gdf = DataFrame(geometry=shp.geometry, val=df)
            rename!(gdf,:val => Symbol(agg))
            centerpoints = [centroid(gdf.geometry[i]) for i in 1:length(gdf.geometry)]
            cr = coordinates.(centerpoints)
            p1 = Plots.plot(gdf.geometry;fillcolor=:cividis,fillalpha=0.5)
            #anns = string.(select(gdf, 2))
            anns = getproperty(gdf,Symbol(agg))
            Plots.annotate!(first.(cr),last.(cr),
                Plots.text.(anns, 7, color, :center, 
                    halign=:center, rotation=35.0;family="Computer Modern"))
        

            display(p1)
            return gdf
        end
            
        """
        lat lon reverse
        """
        function reverse_coords(polygon)
            # Extract the coordinates from the polygon
            # crds = GeoInterface.coordinates.(polygon)
            # # Reverse each pair of coordinates
            # reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
            if polygon isa DataFrame
                out = []
                for i in 1:length(polygon.geometry)
                    # Get the coordinates of the polygon
                    coords = GeoInterface.coordinates(polygon.geometry[i])
                    # Flatten the list of points
                    ptc = coords[1]
                    ptc_tuples = [tuple(pt...) for pt in ptc]
                    prev = [(Y,X) for (X,Y) in ptc_tuples]
                    reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
                    push!(out,reversed_polygon)       
                end
                pout = vcat(out...)   
                return pout
            end
            # coords = GeoInterface.coordinates.(gd.geometry)
            coords = GeoInterface.coordinates(polygon)
            # Flatten the list of points
            #points = reduce(vcat, coords[1])
            ptc = coords[1]
            #ptc = GeoInterface.coordinates.(points)
            ptc_tuples = [tuple(pt...) for pt in ptc]
            prev = [(Y,X) for (X,Y) in ptc_tuples]
            #
            #rpt = ArchGDAL.createpoint.(prev)
            #reversed_polygon = ArchGDAL.createpolygon(coords)
            reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
            
            # Create a new polygon with the reversed coordinates
            #reversed_polygon = ArchGDAL.createpolygon.(reversed_crds)
            #reversed_polygon = ArchGDAL.createpolygon.(prev)

            return reversed_polygon
        end

        """
        flip x and y coords.
        df.geometry .= revcrds(df)

        npt = ArchGDAL.createpoint.(
            ArchGDAL.gety.(stsub.geometry,0),
            ArchGDAL.getx.(stsub.geometry,0))
        """
        function revcrds(stsub::DataFrame)
            npt = ArchGDAL.createpoint.(
                ArchGDAL.gety.(stsub.geometry,0),
                ArchGDAL.getx.(stsub.geometry,0))
            return npt
        end

        
        """
        function to add time_bounds(filename::String)
            on a NetCDF file inplace
        """
        function add_time_bounds(filename::String)
            # Open your NetCDF file
            ds = NCDatasets.Dataset(filename, "a")
            
            # Check if the TIME variable exists
            if haskey(ds, "time")
                # Get the time variable
                time = ds["time"]
                
                #renameVar(ds, "bnds", "time_bnds")
                #renameDim(ds, "time_bnds","bnds")


                # Create a new dimension called 'bnds'
                # NCDatasets.defDim(ds, "bnds", 2)

                # Create a new variable called 'time_bnds'
                time_bnds = try
                    #NCDatasets.defVar(ds, "time_bnds", time, ("time", "bnds"))
                    NCDatasets.defVar(ds, "time_bnds", time, ("time", "x", "y"))
                    #lon × lat × time
                    #NCDatasets.defVar(ds, "time_bnds", time, ("time", "lon", "lat"))
                catch e
                    @error "Error creating the time_bnds variable: $e"
                end

                # Calculate the bounds for each time value
                # # This is just an example, you might need to adjust this for your specific case
                # for (i, t) in enumerate(time)
                #     time_bnds[i, 1] = t - 0.5
                #     time_bnds[i, 2] = t + 0.5
                # end
            else
                println("The time variable does not exist in the NetCDF file.")
            end
            
            # Close the NetCDF file
            NCDatasets.close(ds)
        end

        """
        subset, replace missings, calc mean and plot
        add geodataframe
        rmeanplot(ras::Raster,yr=2000,ge::DataFrame;contour=true,kwargs...)
        """
        function rmeanplot(ras::Raster,yr=2000;
            df::DataFrame=DataFrame(),
            contour=true,kwargs...)
            kras = ras[Ti(Where(x->year.(x)==yr))]
            kras = replace_missing(kras,0)
            kras = mean(kras,dims=(Ti))
            if contour
                Plots.contourf(kras,
                    margins=5mm;
                    #colorbar_title = "mean values of the year $yr",
                    #colorbar_title_padding = 5,
                    kwargs...)
            else
                Plots.plot(kras;
                c = :matter,
                margins = 5mm,
                #colorbar_title = "mean values of the year $yr",
                kwargs...)
            end
            if !isempty(df)
                Plots.plot!(df.geometry;
                fillcolor=:transparent,
                linecolor=:lightgrey,
                linesize=1.9,linestyle=:dashdot)
            else
                title!("Mean values of the year $yr")
            end

        end

        """
        takes mean of 3D Array and plots at time axis
        """
        function ratime(filepath::String, key::Symbol;kwargs...)
            raster = Raster(filepath, key=key)

            Plots.plot(raster[X(Near(mean(extrema(lookup(raster, :X))))),
                        Y(Near(mean(extrema(lookup(raster, :Y)))))],
                margins=5mm;kwargs...)
        end

        function ratime(raster::Raster;kwargs...)
            Plots.plot(raster[
                X(Near(mean(extrema(lookup(raster, :X))))),
                Y(Near(mean(extrema(lookup(raster, :Y)))))],
                margins=5mm;kwargs...)
        end

        """
        extract data from a raster at a point
        use `rst.stp` to return a GeoDataFrame from Station Data
        `ras = A[ X(Near(xc)),Y(Near(yc)) ]`
        `df = ncd(ras::Raster)`
        """
        function ncd(ras::Raster)
            tv = parent(Rasters.lookup(ras,:Ti))
            vl = DataFrame(ras',:auto)|>permutedims
            vl[!,1] .= tryparse.(Float64,string.(vl[!,1])) 
            # would return nothing if not possible
            vl[!,1] .= replace(vl[!,1],nothing=>missing)
            #vl[!,1] .= parse.(Float64,string.(vl[!,1]))            
            df = hcat(vl, tv, makeunique=true)
            rename!(df,1=>name(ras),2=>:date)
            df.date = Date.(df.date)
            #df.pre = parse.(Float64,string.(df[!,Not(:date)]))
            return df
        end

        """
        returns x,y cell size of a raster
        ```	
        get_cell_size(raster_path::AbstractString)
        ```
        """
        function get_cell_size(raster_path::AbstractString)
    
            if endswith(raster_path, ".nc")
                  # Open the NetCDF file
                nc_dataset = NCDatasets.Dataset(nc_file, "r")
                
                # Read the x and y variables
                x_var = nc_dataset[x_var_name]
                y_var = nc_dataset[y_var_name]
                
                # Get the cell size from the first two elements of each dimension
                cell_size_x = x_var[2] - x_var[1]
                cell_size_y = y_var[2] - y_var[1]
                
                # Close the NetCDF file
                close(nc_dataset)
                
                return cell_size_x, cell_size_y
            else    
                # Open the raster dataset
                dataset = ArchGDAL.readraster(raster_path)
                
                # Get the geotransform parameters
                geotransform = ArchGDAL.getgeotransform(dataset)
                
                # Extract the cell size from the geotransform parameters
                cell_size_x = geotransform[2]
                cell_size_y = geotransform[end]  # In case pixel is not a square, use geotransform[5]
                
                # Close the dataset
                #ArchGDAL.close(dataset)
                
                return cell_size_x, cell_size_y
            end
        end

        """
        subset Raster by value
        """
        function maskbyvalue(r::Rasters.Raster, value::Union{Float64, Vector{Float64}})
            if value isa Float64
                zm = (r .== value)
            else
                zm = in.(r, Ref(value))
            end
            
            zm = boolmask(zm;)
            A = Rasters.mask(r; with=zm)
            return A
        end

        
    """
    Plots outputraster with annotations from EZG shapefile
    ```
    ezout(r::Raster, ez::DataFrame; kw...)
    ```
    """
    function ezout(r::Raster, ez::DataFrame;leg=true, kw...)
        plot(r, c=cgrad(:matter, rev=true),
            xlabel="",ylabel=""; kw...)
        if leg
            mids = [ GeoDataFrames.centroid(b) for b in ez.geometry ]
            crds = GeoInterface.coordinates.(mids)
            for (N,i) in enumerate(crds)
                x = first(i)
                y = last(i)
                name = ez[:,end][N]
                annotate!(x, y, 
                    text("TG $name", 10, :black,
                    :bottom, :left, :sans;) )
            end
            plot!(mids,color=:black)
        end
        plot!(ez.geometry,fillcolor=:transparent,linecolor=:black,linewidth=1.5)
    end

    """
    ```
    polygonize_gdal(gdal_polygonize_path::String, input_raster_path::String, output_file_path::String, output_format::String="GeoJSON")
    "ESRI Shapefile"
    ```
    see also 
    ```
    gdal-config  --version 
    ```
    """
    function polygonize_gdal(input_raster_path::String, output_file_path::String, output_format::String="GeoJSON")
        #gdal_polygonize_path = Conda.PREFIX*"/Scripts/gdal_polygonize.py"
        gdal_polygonize_path = Conda.script_dir(Conda.ROOTENV)*"/gdal_polygonize.py"
        cmd = pipeline(`python $(gdal_polygonize_path) $(input_raster_path) -f $(output_format) $(output_file_path)`)
        #cmd = pipeline(`mamba run -n ds python $(gdal_polygonize_path) $(input_raster_path) -f $(output_format) $(output_file_path)`)
        run(cmd)
    end
    

    #end #endof funcs

end #endof module rst

#end #end of endof module rst


# fnames = names(rst, all=true)
# for submodule in fnames
#     @eval import Main.rst.$submodule
# end

# #this is necessary to use the modules in the REPL
# using DataFrames, CSV, Statistics, Dates, Distributions,StatsPlots, Plots.PlotMeasures
# using DelimitedFiles, Grep, Printf, PrettyTables
# using Rasters, ArchGDAL, Grep
# import NCDatasets




# """
# pycall function to polygonize a raster
# polygonize_raster(input_raster_path::String, output_shapefile_path::String;epsg=25832)
# """
# function polygonize_raster(input_raster_path::String, output_shapefile_path::String;epsg=25832)
#     gdal = pyimport("osgeo.gdal")
#     ogr = pyimport("osgeo.ogr")
#     osr = pyimport("osgeo.osr")

#     # Open the raster dataset
#     dataset = gdal.Open(input_raster_path)

#     # Get the first band
#     band = dataset.GetRasterBand(1)

#     # Get the "ESRI Shapefile" driver
#     driver = gdal.GetDriverByName("ESRI Shapefile")

#     # Create a new shapefile dataset
#     out_ds = driver.Create(output_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)

#     # Create a spatial reference object
#     srs = osr.SpatialReference()
#     srs.ImportFromEPSG(epsg)

#     # Create a new layer
#     layer = out_ds.CreateLayer("polygonized", srs, ogr.wkbPolygon)

#     # Polygonize the raster
#     try
#         gdal.Polygonize(band, py"None", layer, -1, [], callback=py"None")
#         @info "new shapefile created at: $output_shapefile_path ..."
#     catch
#         @error "gdal.Polygonize failed ..."
#     end
    
#     # Close the dataset to write it to the disk
#     out_ds = py"None"
# end
