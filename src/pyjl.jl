platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
        
if platform == "windows"
    src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
    macro pj() pt=raw"C:\Users\Public\Documents\Python_Scripts\julia\pyjl.jl";include(pt);end
else
    src_path = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
    macro pj() pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/pyjl.jl";include(pt);end
end   

module pyjl
    using PyCall
    using PyPlot
#    include("win/smallfuncs.jl") only in locdir not here
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
        ddf = wa.pydf(ddf)
        
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

    function xrp(x::AbstractString; maskval=0, lyr=0)
        xr = pyimport("xarray")
        plt = pyimport("matplotlib.pyplot")
        x = nconly(x)|>last
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





end

@info "running using PyCall; pygui(true) now..."
using PyCall; pygui(true) 

# ap = pyjl.pyplot_df
# r = pyjl.pyread
# ap(r"Sch"|>r;log=true)
# r"Sch"|>r|>p->ap(p;log=true)
# r"Unter"|>r|>p->ap(p;log=true)