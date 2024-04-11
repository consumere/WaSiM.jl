# 
function ftp(z::Regex)
        """
        first match of regex and qoutjl
        """
        theplot(first(
        filter(x->occursin(z,x),
        filter(x->endswith(x,"qoutjl"),readdir()))
        )
        )
    end

    