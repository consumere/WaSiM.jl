# 
function fsize()
    """
    ALL names and size in MB via readdir
    """
    files = filter(x -> isfile(x), readdir())
    fzs = [(file, filesize(file) / 2^20) for file in files]
    tot = round(sum(map(x->x[2],fzs));digits=3)
    fzs = sort(fzs, by = x -> x[2], rev = true)
    odf = rename(DataFrame(fzs),["file","size"])
    DataFrames.metadata!(odf, "Total Size", tot, style=:note)
    printstyled("Total Size: $tot MB\n",color=:green)
    return(odf)
end


"""
Read the text file, preserve line 1 as header column
Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
This avoids reading the entire file into memory at once, 
which can be more memory-efficient for large datasets.
"""
