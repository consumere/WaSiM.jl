# 
function dfonly(x::Vector{Any})
        
        z = try filter(f->!occursin(r"^wq|xml|nc|png|svg|jpg",f),x)
        catch
            @warn "vector should be of type string.."
            return
        end
        # z = filter(file -> occursin(x,file), 
        # readdir()[broadcast(x->!occursin(r"^wq|xml|nc|png|svg|jpg",x),readdir())]);
        return(z)
    end

    """
    removes empty TS recursively; 
    use with caution!
    """
    