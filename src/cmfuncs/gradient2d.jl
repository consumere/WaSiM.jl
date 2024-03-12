# 
function gradient2d(z)
            dy, dx = size(z)
            gx = [z[i, j+1] - z[i, j] for i in 1:dy, j in 1:dx-1]
            gy = [z[i+1, j] - z[i, j] for i in 1:dy-1, j in 1:dx]
            return gx, gy
        end
    
        # Calculate the gradient of z
        gx, gy = gradient2d(z)
    
        # Create a 