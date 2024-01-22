#would be redefinition of macro pj
# platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
        
# if platform == "windows"
#     src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
#     macro pj() pt=raw"C:\Users\Public\Documents\Python_Scripts\julia\pyjl.jl";include(pt);end
# else
#     src_path = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
#     macro pj() pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/pyjl.jl";include(pt);end
# end   

module pyjl
    using PyCall
    using PyPlot
    include("win/smallfuncs.jl")
    #@info "dont forget using PyCall; pygui(true) in vscode!"

    function looks_like_number(str::AbstractString)
        try
            parse(Float64, str)
            return true
        catch
            return false
        end
    end

    function pyread(x)
        """
        pyreader, reads all as stings, conversion later.
        """
        pd = pyimport("pandas")
        df = pd.read_table(x, 
            engine="c",
            verbose=true,
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
        #PyPlot.rc("font",**{"family":"serif","serif":["cmr10"]}) #thats the one
        PyPlot.rc("font", family="serif", serif=["cmr10"])
        PyPlot.figure()
        PyPlot.set_cmap("cividis")
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
    simp=xr.open_dataset("tas_cor_raw.nc")
    tl="D:/remo/cordex/eobs/v28/tas/tas_obs.nc"
    obsh=xr.open_dataset(tl)
    pyjl.doyplot(simh,simp,obsh)
    keys have to be the same!

    ad = xr.open_dataset(a)
    m = ad.keys()|>collect|>last    
    """
    function doyplot(simh, simp, obsh;tosum::Bool=false)
        plt = pyimport("matplotlib.pyplot")
        k = simh.keys()|>collect|>last

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
        PyPlot.rc("font", family="serif", serif=["cmr10"])
        plt.figure()
        # plt.plot(grouped_simh, label=join("\$"*k*"_{sim,h}\$"))
        # plt.plot(grouped_simp, label=join("\$"*k*_"{sim,p}\$"))
        # plt.plot(grouped_obsh, label=join("\$"*k*"_{obs,h}\$"))
        # Set plot title and limits
        # Plot the mean temperature for simh
        plt.plot(grouped_simh, label="\$T_{sim,h}\$")

        # Plot the mean temperature for simp
        plt.plot(grouped_simp, label="\$T_{sim,p}\$")

        # Plot the mean temperature for obsh
        plt.plot(grouped_obsh, label="\$T_{obs,h}\$")
        plt.title("Historical modeled and observed and predicted $k")
        plt.xlim(0, 365)
        # Add grid
        plt.gca().grid(alpha=0.3)
        # Add legend
        plt.legend()
        # Show the plot
        plt.show()
    end

end #end of module

@info "running using PyCall; pygui(true) now..."
using PyCall; pygui(true) 

#@doc PyPlot.boxplot
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