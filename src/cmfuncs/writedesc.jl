# 
function writedesc(file, table)
        CSV.write(file, describe(table), transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end

    #wc -l in julia:
    