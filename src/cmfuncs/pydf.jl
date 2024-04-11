# waread3($fn).reset_index(drop=False)
function pydf(py_df::PyObject)
    #     """
    #     Convert each column of the Python DataFrame to a Julia array
    #     """
    #     # fn = filename
    #     # pyo = py"""waread3($fn).reset_index(drop=False)"""
    #     # pdf = wa.pydf(pyo)
    #     # names(pdf)
    #     # dfp(pdf)
    #     col_names = py_df.columns  # Get the column names from the Python DataFrame
    #     col_names = convert(Array, col_names)
    #     col_arrays = [convert(Array, py_df[col]) for col in col_names]
    #     jdf = DataFrame(col_arrays, :auto)
    #     #size(jdf)
    #     fn = try
    #         py_df.filename
    #     catch
    #         @info "no filename present"
    #     end

    #     metadata!(jdf, "filename", fn, style=:note);
    #     rename!(jdf, col_names);
    #     return jdf
    # end

    