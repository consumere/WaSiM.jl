# 
function homrc()
    if platform == "unix"
        cd("/mnt/d/Wasim/regio/out/rc200");
    else
        cd("D:/Wasim/regio/out/rc200");
    end
    println("you are here: ",pwd())
    fdi()
end


