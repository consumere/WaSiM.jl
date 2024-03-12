# 
function mkestream(r::Raster;wasim=true,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            # replace!(z, missval=>missing)
        else
            z = r.data
            # reverse!(z,dims=1)
            # reverse!(z,dims=2)
            # replace!(z, missval=>missing)
        end
        #ex = extrema(z|>skipmissing)
        
        # Define a 