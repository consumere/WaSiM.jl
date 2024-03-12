# 
function pw()
        """
        cb(pwd())
        """
        pt = replace(pwd(),"\\"=> "/")
        println("$pt in clipoard...")
        pt |> clipboard

    end

    pww = pw

    