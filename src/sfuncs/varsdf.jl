# 
function varsdf()
        # Get variable info as a string
        varinfo_str = InteractiveUtils.varinfo(;sortby=:size, minsize=1)

        # Split the string into lines
        lines = split(string(varinfo_str), '\n')

        # Remove the header and footer lines
        lines = lines[3:end-1]

        # Join the lines with newline characters
        csv_str = join(lines, '\n')

        # Parse the CSV string into a DataFrame
        df = CSV.read(IOBuffer(csv_str), DataFrame, header=["name", "size", "summary"], delim='|', ignorerepeated=true)
        sort!(df, :size, rev=true)

        return df
    end



    #if (occursin(Regex(prefix,"i"),filename))
    """
    join([x1,y1],"+.*")
    r"this+.*that"
    """
    