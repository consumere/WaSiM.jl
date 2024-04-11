# 
function subset_dataframe_by_mask(df::DataFrame, msk::DataFrame)
        """
        map(typeof, eachcol(df)) #check types of cols
        msk = broadcast(x->typeof(x)==Vector{Float64},df)
        """
        # Get column names that satisfy the condition
        columns_to_keep = names(df)[collect(msk[1, :])]
        # Subset DataFrame using the mask
        subset_df = select(df, columns_to_keep)
        return subset_df
    end

    