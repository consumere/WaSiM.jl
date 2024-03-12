# 
function pyread(x::Union{String,Regex})
    """
    pyreader, reads all as stings, conversion later.
    """
    if x isa Regex
        x = dfonly(x)|>first
    end
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

    