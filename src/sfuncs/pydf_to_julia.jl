# 
function pydf_to_julia(py_df::PyObject)
        # Convert each column of the Python DataFrame to a Julia array
        py_df = py_df.reset_index(inplace=false)
        col_names = py_df.columns  # Get the column names from the Python DataFrame
        col_arrays = try 
            [convert(Array, py_df[col]) for col in col_names]
        catch
            @error "error in converting!"
            return 
        end
        cas = [convert(Vector{Float64}, z) for z in col_arrays]
        # Create a Julia DataFrame using the converted arrays and column names
        julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, cas))
        #julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, col_arrays))    
        return julia_df
    end 
    
    """
    Convert each column of the Python DataFrame to a Julia array
    """
    