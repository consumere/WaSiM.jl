# 
function readf(x::Regex)
        """
        Read the text file, preserve line 1 as header column
        """
        x=glob(x)|>first
        df = CSV.read(x, DataFrame, delim="\t",header=1,
        silencewarnings=true,types=Float64)
        dropmissing!(df,1)
        for i in 5:size(df,2)
            df[!,i]=replace(df[!,i],-9999.0 => missing)
        end 
        # map to int for dates
        for i in 1:3
            df[!,i]=map(x ->Int(x),df[!,i])
        end
        #and parse dates...
        df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end

    """
    waread2 on regex
    kwargs... passed to CSV.File
    Read the text file, preserve line 1 as header column
    Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
    This avoids reading the entire file into memory at once, 
    which can be more memory-efficient for large datasets.
    """
    