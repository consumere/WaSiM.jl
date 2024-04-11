# 
function globdf(x::AbstractString)
        """
        greps ts from current dir Regex
        """
        # Use a regex to match x and exclude the unwanted extensions
        #regex = Regex(x * "(?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg))","i")
        #regex = Regex(x * "(?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg))")
        filter(file -> (occursin(
            Regex(x,"i"),file) & 
        (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
        ), readdir())
        #regex = Regex("\\Q$x*?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg)","i")
        # Use a list comprehension to filter the files
        #files = [file for file in readdir() if occursin(regex, file)]
    end

   
    