# 
function writedf(file, table)
        CSV.write(file, table, transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end

    