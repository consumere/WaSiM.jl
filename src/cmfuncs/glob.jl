# 
function glob(x::Regex)
        """
        greps from current dir Regex
        """
        filter(file -> occursin(x,file), readdir())
    end

    