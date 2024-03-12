# 
function pydf(py_df::PyObject)
        #py_df.reset_index(inplace=true)
        py_df = py_df.reset_index(inplace=false)
        col_names = py_df.columns  # Get the column names from the Python DataFrame
        col_names = convert(Array, col_names)
        col_arrays = try 
                    [convert(Array, py_df[col]) for col in col_names]
                catch
                    @error "error in converting!"
                    return 
                end
        #col_arrays = [convert(Array, df[col]) for col in col_names]
        cas = [convert(Vector{Float64}, z) for z in col_arrays]
        #julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, cas))
        jdf = DataFrame(cas, :auto)
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

    """
    Convert DataFrame Column to a Vector
    vec(Matrix(select(df,col)))
    """
    