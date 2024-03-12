# 
function vars(pt::AbstractString)
        #varinfo(Core,r".*field.*")
        #varinfo(Main,r".*load*")
        varinfo(Main,Regex(".*pt*"))
    end

    #if (occursin(Regex(prefix,"i"),filename))

    