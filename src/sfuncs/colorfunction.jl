# 
function colorfunction(v::Vector)
        max_val = maximum(v)  # calculate the maximum once, outside the loop
        f = 1/length(v)
        return [x < f ? :red :
                x < 2*f ? :teal :
                x < 3*f ? :yellow :
                x < 4*f ? :green :
                x < 5*f ? :cyan :
                x < 6*f ? :purple :
                x < 7*f ? :orange :
                x < 8*f ? :magenta :
                x < 9*f ? :pink :
                :blue for x in v ./ max_val]  # normalize v here
    end
