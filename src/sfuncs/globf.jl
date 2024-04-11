# 
function globf(x::AbstractString)
        """
        greps first match current dir Regex
        """
        first(filter(file -> occursin(Regex(x,"i"),file), readdir()))
    end
    
    #readall = loadalldfs
    dfread = readf
    readmeteo = waread
    loaddf = readdf


    