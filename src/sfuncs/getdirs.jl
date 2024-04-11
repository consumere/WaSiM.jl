# 
function getdirs(kw...)
    printstyled(filter(x->isdir(x),readdir(kw...)),color=:blue)
end

"""
dfroute(;ofl="route.txt")
reads from routeg(infile, ofl) and returns a DataFrame with the following columns:
    - sim: simulated flow
    - obs: observed flow
    - name: name of the station
"""
