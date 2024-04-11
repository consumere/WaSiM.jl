# 
function rename_duplicates(df::DataFrame)
        # loop over the column names
        for (i, name) in enumerate(names(df))
            # check if the name contains any repeated words
            if occursin(r"(\w+)[ _]\1", name)
                # replace the duplicates with a single word
                new_name = replace(name, r"(\w+)[ _]\1" => s"\1")
                # rename the column
                rename!(df, name => new_name)
            end
        end
        # return the modified dataframe
        return df
    end


    