module pyjl
    using PyCall
    using PyPlot
    include("smallfuncs.jl")
    #@info "dont forget using PyCall; pygui(true) in vscode!"
    @pyimport hydroeval as he


    function looks_like_number(str::AbstractString)
        try
            parse(Float64, str)
            return true
        catch
            return false
        end
    end

    """
    smaller looks_like_number
    """
    function llnp(str::Any) 
        #ternary operator makes no sense
        #flag = false ? !isnothing(tryparse(Float64, str)) : true
        #return
        !isnothing(tryparse(Float64, str))
    end

    """
    smaller looks_like_number with regex
    """
    function lln(str::Any)
        is_number = occursin(r"^[-+]?\d+(?:\.\d+)?$", str) || occursin(r"^[-+]?\d+(?:\.\d+)?e[-+]?\d+$", str)
        return is_number
    end

    function pyread(x)
        """
        pyreader, reads all as stings, conversion later.
        """
        pd = pyimport("pandas")
        df = pd.read_table(x, 
            engine="c",
            #verbose=true, #depricated
            low_memory=false,
            header=0,
            skipinitialspace=true,
            dtype="str",              #new!
            na_values=[-9999]
            )
        col_names = df.columns  # Get the column names from the Python DataFrame
        col_names = convert(Array, col_names)
        col_arrays = [convert(Array, df[col]) for col in col_names]
        filtered_rows = broadcast(x->looks_like_number(x),col_arrays[1])
        
        df = DataFrame(col_arrays, :auto)
        
        df = df[filtered_rows, :]
        
        rename!(df, col_names)
        if "YY" ∉ names(df)
            println("Column 'YY' not found in the CSV file.")
            return nothing
        end

        function parse_date(row)
            try
                year = parse(Int, row[1])
                month = parse(Int, row[2])
                day = parse(Int, row[3])
                return Date(year, month, day)
            catch
                return missing  # Return missing if parsing fails
            end
        end
        
        df.date = map(parse_date, eachrow(df))
        
        df = select(df, Not(1:4))
        
        for x in names(df)
            if looks_like_number(x)
                newname=replace(x,"$x"=>"C$x", count=1)
                rename!(df,Dict(x=>newname))
            end
        end

        #map(y->typeof(y),eachcol(df))

        # Iterate over column names
        for colname in names(df)
            # Check if the column type is not Date
            if eltype(df[!, colname]) != Date
                df[!, colname] .= tryparse.(Float64, df[!, colname])
            end
        end

        DataFrames.metadata!(df, "filename", x, style=:note)
        return df 
    end

    function pydf_to_julia(py_df::PyObject)
        """
        no transposing
        """
        # Convert each column of the Python DataFrame to a Julia array
        col_names = py_df.columns  # Get the column names from the Python DataFrame
        col_arrays = [convert(Array, py_df[col]) for col in col_names]
        # Create a Julia DataFrame using the converted arrays and column names
        julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, col_arrays))    
        return julia_df
    end 
    
    function pydf(py_df::PyObject)
        """
        Convert each column of the Python DataFrame to a Julia array
        """
        # fn = filename
        # pyo = py"""waread3($fn).reset_index(drop=False)"""
        # pdf = wa.pydf(pyo)
        # names(pdf)
        # dfp(pdf)
        col_names = py_df.columns  # Get the column names from the Python DataFrame
        col_names = convert(Array, col_names)
        col_arrays = [convert(Array, py_df[col]) for col in col_names]
        jdf = DataFrame(col_arrays, :auto)
        #size(jdf)
        fn = try
            py_df.filename
        catch
            @info "no filename present"
        end

        metadata!(jdf, "filename", fn, style=:note);
        rename!(jdf, col_names);
        return jdf
    end

    function pyread_meteo(s::AbstractString;hdr=0)
        pd = pyimport("pandas")
        ddf = pd.read_csv(s, delim_whitespace=true, 
            header=hdr,
            na_values=-9999,
            low_memory=false,
            verbose=true)
        #ddf.filename=basename(s)
        ddf = pydf(ddf)
        
        #dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
        dt = select(ddf,1:3)
        Z = []
        for i in eachcol(dt)
            col = tryparse.(Int, i)
            push!(Z,col)
        end
        dt = DataFrame(Z,:auto)
        
        msk = tryparse.(Int, ddf.YY)
        ddf = ddf[findall(!isnothing,msk),:]
        dt = dt[findall(!isnothing, dt.x1),:]
        dt = Date.(map(k -> join(k, "-"), eachrow(dt[:, 1:3])))
        nd = hcat(ddf[!, 5:end], dt)
    
        for col in names(nd)[(eltype.(eachcol(nd)) .<: String)]
            nd[!, col] .= tryparse.(Float64, col)
        end
        rename!(nd, ncol(nd) =>"date")
        metadata!(nd, "filename", s, style=:note);
        return nd
    end

    function pyread_old(s::AbstractString;hdr=0)
        pd = pyimport("pandas")
        ddf = pd.read_csv(s, delim_whitespace=true, 
            header=hdr,
            na_values=-9999,
            low_memory=false,
            verbose=true)
        #ddf.filename=basename(s)
        ddf = pydf(ddf)
        
        dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
        
        nd = hcat(ddf[!, 5:end], dt)

        for col in names(nd)[(eltype.(eachcol(nd)) .<: String)]
            nd[!, col] .= tryparse.(Float64, col)
        end
        rename!(nd, ncol(nd) =>"date")
        metadata!(nd, "filename", s, style=:note);
        return nd
    end

    function waread3_py(x, flag=true)
        # use py"""...""" syntax to access the Python function
        py"""
        import pandas as pd
        import datetime

        def waread3(x, flag=True):
            if flag:
                df = pd.read_csv(x, delim_whitespace=True, header=0,
                                na_values=-9999, verbose=True,engine='c')
                if 'YY' not in df.columns:
                    print("Column 'YY' not found in the CSV file.")
                    return None
                if df.iloc[:, 0].dtype == 'object':
                    print('First column contains strings, subsetting to Int...')
                    df = df[~df.iloc[:, 0].str.contains("[A-z]|-", na=False)]
                source_col_loc = df.columns.get_loc('YY')        
                df['date'] = df.iloc[:, source_col_loc:source_col_loc +
                                    3].apply(lambda x: "-".join(x.astype(str)), axis=1)
                df = df.iloc[:, 4:]
                df['date'] = pd.to_datetime(df['date'])
                #df.set_index('date', inplace=True)
                df.iloc[:,0:-2] = df.iloc[:,0:-2].apply(lambda x: x.astype(float), axis=1)
                df.filename = x
                print(df.filename,"done!")
                return df
            else:
                print('Date parse failed, try reading without date transformation...')
                df = pd.read_csv(x, delim_whitespace=True, comment="Y", skip_blank_lines=True).dropna()
                df.filename = x
                print(df.filename,"done!")
                return df
        """
        # call the Python function with the arguments and return the result
        return py"waread3($x, $flag)"
    end

    function xrfacets(a::AbstractString;maskval=0)
        """
        uses xarray to plot a 4xn grid of wasim stacks.
        """
        xr  = pyimport("xarray")
        plt = pyimport("matplotlib.pyplot")
        ad = xr.open_dataset(a)
        m = ad.keys()|>collect|>last    
        ad[m].where(ad[m]>maskval).transpose().plot(
            col="t",
            col_wrap=4,
            robust=true,
            cmap=plt.cm.RdYlBu_r);
        #return p1.fig #now we can see the plot inside vscode, if PyCall gui works
        # pygui_start()
        #p1.fig.show()
        #display(p1)
        plt.show()     #thats the right way.
    end

    function xrplot(x::AbstractString;maskval=0,lyr=0)
        x = nconly(x)|>last
        py"""
        from xarray import open_dataset
        from matplotlib.pyplot import show
        maskval = float($maskval)
        dx = open_dataset($x,mask_and_scale=True).isel(t=$lyr).transpose().to_array()
        dx.where(dx.values>maskval).plot(cmap="cividis")
        show()
        """
    end

    #ENV["LD_PRELOAD"] = .so to the missing c++ library...

    
    #filter(z -> endswith(z,".nc"),readdir()[map(x->occursin(rx,x),readdir())])


    """
    pyjl.xrp(r"sb.+06") 
    xrp(x::Union{String,Regex}; maskval=0, lyr=0)
            if isa(x,Regex)
                x = nconly(x)|>last
            end
    """
    function xrp(x::Union{String,Regex}; maskval=0, lyr=0)
        if isa(x,Regex)
            #x = nconly(x)|>last
            v = readdir();
            z = v[map(k->occursin(x,k),v)]
            println(z)
            z = z[map(k->endswith(k,"nc"),z)][end]
            x = z
        end
        xr = pyimport("xarray")
        plt = pyimport("matplotlib.pyplot")
        dx = xr.open_dataset(x,mask_and_scale=true)
        m = dx.keys()|>collect|>last    
        dx[m].where(dx[m]>maskval).isel(t=lyr).transpose().plot(cmap="turbo")
        #.isel(t=lyr).transpose().to_array()
        #p1 = dx.where(dx.values>maskval).plot(cmap="turbo")
        #@pyimport matplotlib.pyplot as plt #plt as const
        ti=basename(x)
        plt.title(ti)
        plt.show()
    end

    function pyplot_df(df::DataFrame;log=false)
        x = df.date
        ln = (filter(x -> !occursin(r"date|month|year", x), names(df)))
        for col in ln
            y = df[!, Symbol(col)]
            PyPlot.plot(x, y, label=col)
        end

        if log
            PyPlot.yscale("log")
        end
    
        PyPlot.xlabel("Date")
        PyPlot.ylabel("")
        PyPlot.legend()
        ti = only(values(DataFrames.metadata(df)))
        PyPlot.title(ti)
        PyPlot.grid(true)
    end

    function reorder_df(df::DataFrame)
        """
        date to last position
        """
        df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        return(df)
    end

    function pybar(pt::Union{String,Regex,DataFrame};
        fun=mean, box=false,
        log=false)

        if pt isa String
            df = waread(pt)
        elseif pt isa Regex
            df = waread(pt)
        else
            df = pt
        end

        #x = df.date
        #ln = (filter(x -> !occursin(r"date|month|year", x), names(df)))
        # Add a Year column
        df.Year = year.(df.date)
        select!(df,Not(:date))
        
        ln = (filter(x -> !occursin(r"date|month|year"i, 
            x), names(df)))
        # Aggregate data by year and column
        # df_agg = DataFrames.combine(
        #     groupby(df, :Year), 
        #     ln .=> fun .=> ln);
    
        df_agg = combine(groupby(df, :Year), 
            (ln .=> fun .=> ln)...)

        if box
            # Get unique years
            years = unique(df_agg.Year)

            # Create a new figure
            figure()

            # Iterate over columns
            for (i, col) in enumerate(ln)
                # Create a new subplot for each column
                subplot(length(ln), 1, i)
                title(col)

                # Create a boxplot for each year
                boxplot([df_agg[df_agg.Year .== year, 
                    col] for year in years], 
                    labels=years)
            end
            # for (i, col) in enumerate(ln)
            #     PyPlot.boxplot(
            #         df_agg[:,col]
            #     )
            # end
            #PyPlot.xlabel("Year")
        else           
            # Create grouped bar plot
            for (i, col) in enumerate(ln)
                PyPlot.bar(
                    df_agg.Year .+ (i-1)*0.2, 
                df_agg[:, col], 
                width=0.2, 
                label=col, 
                align="center")
            end
            # for col in ln
            #     y = df_agg[!, Symbol(col)]
            #     #PyPlot.plot(x, y, label=col)
            #     PyPlot.bar(
            #         df_agg.Year, # .- 0.2, 
            #         y, 
            #         width=0.4, 
            #         label=col, 
            #         align="center")
            # end
        end

        if log
            PyPlot.yscale("log")
        end
        
        PyPlot.ylabel("")
        PyPlot.legend()
        ti = only(values(DataFrames.metadata(df)))
        PyPlot.title(ti)
        PyPlot.grid(true)
    end

    
    """
    monthly hydrgraph plot using PyPlot
    selection of first column by default
    x::Union{Regex,String,DataFrame}; col = 1,
        leg = "best", logy = false)
    """
    function pyhydro(x::Union{Regex,String,DataFrame}; 
        col = 1,
        leg::String = "best", logy = false)
        if isa(x,DataFrame)
        df = (x)
        else
        #df = waread(x)
        df = pyread(x)
        end

        ti = try
        DataFrames.metadata(df)|>values|>only
        catch
        @warn "No basename in metadata!"
        raw""
        end

        df = select(df,Cols(col,:date))
        @info "selected column:",names(df)

        colname = names(df[!,Not(:date)])|>only

        s = filter(x -> !occursin(r"date|year|month"i, string(x)), names(df))
        years = unique(year.(df.date))
        years_str = string.(years)

        ylog = logy ? :log : :identity

        mn = [ monthabbr(x) for x in unique(month.(df.date)) ]
        PyPlot.rc("font", family="serif", serif=["cmr10"])
        PyPlot.figure()
        PyPlot.set_cmap("cividis")
        #PyPlot.rc("font",**{"family":"serif","serif":["cmr10"]}) #thats the one
        ## PyPlot.axes.formatter.use_mathtext(true) #<-not available
        for yr in years
            su = filter(row -> year(row.date) == yr, df)
            PyPlot.plot(vec(Matrix(select(su, Not(:date)))), label = yr)
        end
        #PyPlot.xlabel("Time")
        PyPlot.ylabel(colname)
        PyPlot.title(ti)
        # locs = ["best", "upper right", "upper left", "lower left", 
        # "lower right", "right", "center left", "center right", 
        # "lower center", "upper center", "center"]
        PyPlot.legend(loc=leg)
        if logy
            PyPlot.yscale("log")
        end
        PyPlot.xticks(rotation=45)
        PyPlot.grid(true)
        # Add x-axis labels with month abbreviations
        # x_values = 1:length(mn)
        # x_labels = mn
        # PyPlot.xticks(x_values, x_labels)
        # Add x-axis labels with month abbreviations
        #x_values = 1:366/12:366
        x_values = range(1, stop=360, length=12)
        #x_labels = repeat(mn, inner=round(Int, 366/12))
        #x_labels = repeat(mn, inner=31)
        x_labels = mn
        PyPlot.xticks(x_values, x_labels)
        PyPlot.show()
    end

    """
    monthly boxplot plot using PyPlot
    selection of first column by default
    x::Union{Regex,String,DataFrame}; col = 1,
        leg = "best", logy = false, fun=mean)
    with annotations of mean values.
    """
    function pybox(x::Union{Regex,String,DataFrame}; 
        col = 1,
        leg::String = "best", logy = false, fun=mean)
        if isa(x,DataFrame)
        df = (x)
        else
        #df = waread(x)
        df = pyread(x)
        end

        ti = try
        DataFrames.metadata(df)|>values|>only|>basename
        catch
        @warn "No basename in metadata!"
        raw""
        end

        df = select(df,Cols(col,:date))
        @info "selected column:",names(df)

        colname = names(df[!,Not(:date)])|>only

        df.Month = month.(df.date)
        str = [ @sprintf("%02i", x) for x in (df.Month) ]
        month_abbr = ["Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]

        ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))

        grouped_df = groupby(df, :Month)
        
        PyPlot.figure()
        PyPlot.rc("font", family="serif", serif=["cmr10"])
        for (i, group) in enumerate(grouped_df)
            PyPlot.boxplot(group[!, ln[1]], positions=[i], 
            autorange=true,
            #fliers = false,
            notch = true,
            widths=0.7)
        end
        PyPlot.xticks(1:12, month_abbr)
        colname = replace(colname,r"_" => " ")
        ti = replace(ti,r"_" => " ")
        PyPlot.ylabel(colname)
        PyPlot.title(ti)

        if fun != mean
            #funval = DataFrames.combine(grouped_df, ln[1] => fun)[!,end]
            #means = DataFrames.combine(grouped_df, ln[1] => mean)
            #means.val .= funval
            funval = DataFrames.combine(grouped_df, ln[1] => fun)
            #years = unique(year.(df.date))
            #funval[!,end] = funval[!,end] ./ length(years)
            for i in eachrow(funval)
                m = 0 #i[2] #yposition, could also be meansval
                val = i[2]  #i[3]
                PyPlot.annotate(round(val; digits=2), (i.Month, m), 
                    textcoords="offset points", 
                    #xytext=(0,25), 
                    xytext=(0,-10), 
                    ha="center")
            end

        else
            means = DataFrames.combine(grouped_df, ln[1] => mean)
            for i in eachrow(means)
                m = 0 #i[2] #yposition, could also be meansval
                val = i[2]
                PyPlot.annotate(round(val; digits=2), (i.Month, m), 
                    textcoords="offset points", 
                    #xytext=(0,10), 
                    xytext=(0,-10), 
                    ha="center")
            end
        end
        
        return gcf()
    end

    """
    doyplot(simh, simp, obsh)
    xr  = pyimport("xarray")
    sh="d:/remo/qm/tas/simh.nc"
    simh=xr.open_dataset(sh)
    sp="d:/remo/qm/corgrids/tas/tas_cor_raw.nc"
    simp=xr.open_dataset(sp)
    tl="D:/remo/cordex/eobs/v28/tas/tas_obs.nc"
    obsh=xr.open_dataset(tl)
    pyjl.doyplot(simh,simp,obsh)
    keys have to be the same!

    ad = xr.open_dataset(a)
    m = ad.keys()|>collect|>last    
    """
    function doyplot(simh, simp, obsh; tosum::Bool=false)
        plt = pyimport("matplotlib.pyplot")
        k = simh.keys()|>collect|>last
        
        if k != simp.keys()|>collect|>last || k != obsh.keys()|>collect|>last
            @error "varibale keys have to be the same!"
            @info "check with ds.variables|>collect "
            return
        end

        fns = [simh, simp, obsh] #|>map(x->x.filename|>basename)

        # for z in enumerate(fns)
        #     println(z)
        # end
        
        for z in fns
            #idx = obsh.variables|>collect
            idx = z.variables|>collect
            #if (!occursin("longitude", idx)) || (!occursin("latitude", idx))
            #map(x->occursin(x,"time"),z.coords|>collect)
            if (!any(x->occursin(r"longitude|latitude",x),idx))
                lng = ([i for i in idx if startswith(i, "lon") || endswith(i, "x")])|>first
                #lng = join([i for i in idx if startswith(i, "lon") || endswith(i, "x")])
                ltt = ([i for i in idx if startswith(i, "lat") || endswith(i, "y")])|>first
                #tim = ([i for i in idx if startswith(i, "tim")])|>first
                #dimmns = Dict(tim=>"t", lng=>"longitude", ltt=>"latitude")
                dimmns = Dict(lng=>"longitude", ltt=>"latitude")
                z = z.swap_dims(dimmns)
                z = z.rename(dimmns)
                if !haskey(z.coords, "latitude")
                    #z = z.set_index(z, Dict("latitude" => "latitude"))
                    z = z.set_coords("latitude")
                end
                if !haskey(z.coords, "longitude")
                    z = z.set_coords("longitude")
                end

                #z = z.rename_dims(dimmns) #.set_index(dimmns)
                #z = z.swap_dims(dimmns).set_index(dimmns).rename(dimmns)
                #z = z = z.rename(dimmns).set_index(dimmns)
                @info "renamed $idx !"
            end
        end
        

        if tosum
            grouped_simh = simh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
            grouped_simp = simp[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
            grouped_obsh = obsh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
        else
            # Group by dayofyear and calculate mean for each DataFrame
            grouped_simh = simh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").mean()
            grouped_simp = simp[k].mean("longitude").mean("latitude").groupby("time.dayofyear").mean()
            grouped_obsh = obsh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").mean()
        end        
        # Create the figure
        #plt.figure(figsize=(10, 5), dpi=216)
        #PyPlot.rc("font", family="serif", serif=["cmr10"])
        vr = uppercase(first(k))
        plt.rc("font", family="serif", serif=["cmr10"])
        plt.figure()
        plt.plot(grouped_simh, label="\$$vr _{sim,h}\$")
        plt.plot(grouped_simp, label="\$$vr _{sim,p}\$")
        plt.plot(grouped_obsh, label="\$$vr _{obs,h}\$")
        plt.title("Historical modeled and predicted versus observed $(uppercasefirst(k))")
        plt.xlim(0, 365)
        plt.gca().grid(alpha=0.3)
        plt.legend()
        plt.show()
    end

    """
    pycall function to polygonize a raster
    uses osgeo.gdal and osgeo.ogr    
    """
    function polygonize_raster(input_raster_path::String, output_shapefile_path::String)
        gdal = pyimport("osgeo.gdal")
        ogr = pyimport("osgeo.ogr")
        osr = pyimport("osgeo.osr")

        # Open the raster dataset
        dataset = gdal.Open(input_raster_path)

        # Get the first band
        band = dataset.GetRasterBand(1)

        # Get the "ESRI Shapefile" driver
        driver = gdal.GetDriverByName("ESRI Shapefile")

        # Create a new shapefile dataset
        out_ds = driver.Create(output_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)

        # Create a spatial reference object
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(25832)

        # Create a new layer
        layer = out_ds.CreateLayer("polygonized", srs, ogr.wkbPolygon)

        # Polygonize the raster
        gdal.Polygonize(band, py"None", layer, -1, [], callback=py"None")

        # Close the dataset to write it to the disk
        out_ds = py"None"
    end

    """
    dfvec(df::DataFrame, col::Int64)
    """
    function dfvec(df::DataFrame, col::Int64)
        return vec(df[!,col])
    end

    """
    KGE, r, α, β = he.kge(dfvec(ts,1),dfvec(ts,2))
    hekge(x::Union{String,Regex,DataFrame};c1::Int64=1,c2::Int64=2)
    c1 = colindex for simulations
    c2 = colindex for observations\n
    vst = glob(r"qoutjl")
    dfs = map(pyjl.hekge,vst)
    hcat(vcat(dfs...),basename.(vst))
    \nall to one line:\n
    see allkge or do:
    odf = rename(hcat(vcat(map(pyjl.hekge,
        filter( x->endswith(x,"qoutjl"),readdir() ))...),
        basename.(filter(x->endswith(x,"qoutjl"),readdir()))),Dict(5=>"fn"))
    """
    function hekge(x::Union{String,Regex,DataFrame};c1::Int64=1,c2::Int64=2)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            ts = pyread(x)
            #dropmissing!(ts)    #
            printstyled("colorder: $(names(ts))\n",color=:light_green)
        elseif x isa Regex
            ts = pyread(x)
            fn = try 
                collect(DataFrames.metadata(ts))[1][2]|>basename
                catch
                @info "no filename"
                end
            printstyled("$fn loaded\n",color=:light_red)
            #dropmissing!(ts)    #
            printstyled("colorder: $(names(ts))\n",color=:light_green)
        else
            ts = x
            printstyled("colorder: $(names(ts))\n",color=:light_green)
        end    
        KGE, r, α, β = he.kge(dfvec(ts,c1),dfvec(ts,c2))
        kge_matrix = [KGE, r, α, β]
        # Column names see doc he.kge
        col_names = ["KGE", "r", "α", "β"]
        df = DataFrame(Dict(col_names.=>kge_matrix))
        return df
    end

    function kgevec(a::Vector{Float64},b::Vector{Float64})
        return DataFrame(Dict(["KGE", "r", "α", "β"].=>he.kge(a,b)))
    end

    """
    allkge(;dirpath::String=pwd(),suffix::String="qoutjl")
    KGE, r, α, β = he.kge(dfvec(ts,1),dfvec(ts,2))
    """
    function allkge(;dirpath::String=pwd(),suffix::String="qoutjl")
        odf = rename(hcat(vcat(map(hekge,
            filter( x->endswith(x,suffix),
            readdir(dirpath) ))...),
            basename.(filter(x->endswith(x,suffix),
            readdir(dirpath)))),
            Dict(5=>"fn"))
        return odf
    end

    """
    day of year plot using xarray and matplotlib
    """
    function doy(simh, simp, obsh; tosum::Bool=false)
        plt = pyimport("matplotlib.pyplot")
        xr = pyimport("xarray")
        calendar = pyimport("calendar")
        if isa(simh,PyObject)
            simh = simh
        else
            simh = xr.open_dataset(simh)
        end
        if isa(simp,PyObject)
            simp = simp
        else
            simp = xr.open_dataset(simp)
        end
        if isa(obsh,PyObject)
            obsh = obsh
        else
            obsh = xr.open_dataset(obsh)
        end

        
        k = simh.keys()|>collect|>last
        
        if k != simp.keys()|>collect|>last || k != obsh.keys()|>collect|>last
            @error "varibale keys have to be the same!"
            @info "check with ds.variables|>collect "
            return
        end
    
        fns = [simh, simp, obsh] #|>map(x->x.filename|>basename)
        outds = []
        for z in fns
            idx = z.variables|>collect
            if (!any(x->occursin("time",x),z.coords|>collect))
                #ti = [i for i in idx if startswith(i, "lon")]
                z = z.rename(Dict("t"=>"time"))
                #z = z.swap_dims(Dict("t"=>"time"))
                z = z.set_coords("time")
                @info "renamed time!"
            end
            
            #if (!any(x->occursin(r"longitude|latitude",x),idx))
                lng = ([i for i in idx if startswith(i, "lon") || endswith(i, "x")])|>first
                ltt = ([i for i in idx if startswith(i, "lat") || endswith(i, "y")])|>first
                if tosum
                    grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").sum()
                else
                    grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").mean()
                end        
            push!(outds,grouped_ds)
        end
        vr = uppercase(first(k))
        plt.rc("font", family="serif", serif=["cmr10"])
        plt.figure()
        plt.plot(outds[1], label="\$$vr _{sim,h}\$")
        plt.plot(outds[2], label="\$$vr _{sim,p}\$")
        plt.plot(outds[3], label="\$$vr _{obs,h}\$")
        plt.title("Historical modeled and predicted versus observed $(uppercasefirst(k))")
        plt.xlim(0, 365)
        plt.gca().grid(alpha=0.3)
        plt.legend()
        month_labels = ["Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
       # Set the x-axis tick locations to be the middle day of each month
        #tick_locations = [calendar.monthrange(2022, i)[1] // 2 + sum(calendar.monthrange(2022, j)[1] for j in range(1, i)) for i in range(1, 13)]
        tick_locations = [15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
       
        plt.xticks(tick_locations)

        # Update the tick labels to display the month names
        #month_labels = [calendar.month_abbr[i] for i in range(1, 13)]
        plt.gca().set_xticklabels(month_labels) 

        plt.show()
    end

    """
    reads and plots a netcdf file using xarray and matplotlib
    """
    function fdoy(ds::Union{String,PyObject}; txy::Bool=false, tosum::Bool=false, ti::String="")
        plt = pyimport("matplotlib.pyplot")
        xr = pyimport("xarray")
        if ds isa String
            z = xr.open_dataset(ds)
        else
            z = ds
        end
        
        k = z.keys()|>collect|>last
        @info "key is $k"

        idx = z.variables|>collect
        if (!any(x->occursin("time",x),z.coords|>collect))
            #ti = [i for i in idx if startswith(i, "lon")]
            z = z.rename(Dict("t"=>"time"))
            #z = z.swap_dims(Dict("t"=>"time"))
            z = z.set_coords("time")
            @info "set_coords from t to time!"
        end
        if txy
            lng = "x"
            ltt = "y"
            z = z.assign_coords(t=z.time)
        else
            lng = ([i for i in idx if endswith(i, "x") || startswith(i, "lon")] )|>first
            ltt = ([i for i in idx if endswith(i, "y") || startswith(i, "lat")] )|>first
        end
        
        if tosum
            grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").sum()
        else
            grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").mean()
        end        
        vr = uppercase(first(k))
        nm = try 
            basename(ds)
        catch
            ""
        end
        plt.rc("font", family="serif", serif=["cmr10"])
        plt.figure()
        plt.plot(grouped_ds, label="\$$nm _{$vr}\$")
        
        #filter(x -> occursin(r"name", string(x)), collect(keys(z.attrs)))
        #filter(x -> occursin(r"on", string(x)), collect(keys(z.attrs)))

        plt.title("$(uppercasefirst(nm))"*ti)
        plt.xlim(0, 365)
        plt.gca().grid(alpha=0.3)
        plt.legend()
        plt.show()
    end

    """
    export to Main
    exportall(Main.pyjl) 
    """
    function exportall(e::Module)
        #for n in names(@__MODULE__; all=true)
        for n in names(e; all=true)
            if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                @eval export $n
            end
        end
    end

    function importall(e::Module)
        for n in names(e; all=true)
            if Base.isidentifier(n) && n ∉ (Symbol(e), :eval, :include)
                @eval import $(e).$n
            end
        end
    end
    

    # """
    # toMain(Main.pyjl)
    # """
    # function toMain(e::Module)
    #     fnames = names(e, all=true)
    #     for submodule in fnames
    #         @eval import e.$submodule
    #     end
    # end
          
    # fnames = names(Main.pyjl, all=true)
    # for submodule in fnames
    #     @eval import Main.pyjl.$submodule
    # end
    function toMain()
        for submodule in names(pyjl, all=true)
            @eval import Main.pyjl.$submodule
        end
    end
    

    """
    selects from dataframe date an column
    selt(x::DataFrame, col::Any;dtcol=:date)
    """
    function selt(x::DataFrame, col::Any; dtcol=:date)
        df = select(x,Cols(col,dtcol))
        println(names(df))
        return df 
    end

    """
    from a vector of dfs.
    fkge(dfs::Vector{DataFrame})

    see also pxm(map(byear,outd))
    """
    function fkge(dfs::Vector{DataFrame})
        nms = []
        for i in dfs
            x=collect(DataFrames.metadata(i))[1][2]|>basename
            push!(nms, x)
        end
        odf = map(hekge,dfs)
        return rename(hcat(vcat(odf...),nms),Dict(5=>"fn"))
    end

    function qall(;recursive=false)
        if recursive
            files = rglob("qgko")
        else
            files = glob("qgko")
        end
        outdf::Vector{DataFrame} = []
        for file in files
            # Load the file into a DataFrame
            #x = "qgkofab.m6.2010"
            x = file
            try
                df = DataFrame(CSV.File(x, header=1, 
                                    delim="\t",
                                    skipto=366,
                                    ntasks=1,
                                    types = String,
                                    silencewarnings=true,
                                    ignorerepeated=true,
                                    ignoreemptyrows=true,
                                    stripwhitespace=true))
                                    
                #df = CSV.read(x,DataFrame;ntasks=1)
                println(x)
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                
                #mapcols!(x -> parse(Float64, x), dx)
                #mapcols(x -> x.^2, dx)

                #basins = copy(df[1,5:end])
                #AsTable(basins)
                basins = []
                #for i in copy(df[1,5:end])
                for i in names(df)[5:end]
                    push!(basins,i)
                end
                #size(basins)
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins

                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                for col in names(dout)[1:end-1]
                    dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
                end          
                #mapcols!(x -> parse(Float64, x), dout)
                #dout = hcat(dx[!,Cols(r"bas")],dx[:,Not(Cols(r"bas"))])
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                dout.basin = parse.(Int64,dout[!, :basin])
                dout.nm .= replace(file,".\\"=>"")
                push!(outdf,dout)
            catch
                @warn("error! file $x can not be loaded as a DataFrame! ")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
        return(outdf)
    end



    """
    find LOG. R-SQUARE > .4 recursivley
    """
    function findlog(;lb=.4)
        v = qall(;recursive=true)
        #map(x->DataFrames.subset(x,3 => ByRow(<(1))),v)
        k = try 
            map(x->DataFrames.subset(x,3 => ByRow(>(lb))),v)
            catch
                @error "no df on lowerbound! "
                return
        end
        df = try 
            reduce(vcat,k) 
            catch
                @error "vcat failed! "
                return
        end
        
        df = try 
            DataFrames.subset(df,3 => ByRow(<(1.0)))
            catch
                println(first(k))
                @error "no df on upperbound! "
                return
        end
        
        
        return df
    end


end #end of module #endof


@info "running using PyCall; pygui(true) now..."
using PyCall; pygui(true) 

#@doc PyPlot.boxplot
#cd("/mnt/d/Wasim/sinn/out/b1")
#pyjl.pybar(r"sb05")
#pyjl.pybar(r"sb05";fun=sum)
#import Statistics
#pyjl.pybar(r"wind";fun=Statistics.median,box=true)
#pyjl.pybar(r"wind";box=true)
#pyjl.pybar(r"qges";fun=Statistics.median)
#pyjl.pybar(r"qges";fun=z->Statistics.quantile(z,0.9))
# ap = pyjl.pyplot_df
# r = pyjl.pyread
# ap(r"Sch"|>r;log=true)
# r"Sch"|>r|>p->ap(p;log=true)
# r"Unter"|>r|>p->ap(p;log=true)

#pyjl.pyhydro(r"qbas")