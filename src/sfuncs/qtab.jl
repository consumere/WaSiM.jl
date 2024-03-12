# 
function qtab(;todf=true)
    dfs = getq();
    # z = ctlx()
    # if isempty(z)
    #     printstyled("no control file found!\n",color=:light_red)
    #     pretty_table(dfs,header=uppercasefirst.(names(dfs));)
    #     return
    # end
    ofl = "route.txt"
    if isempty(ofl)
        printstyled("no $ofl found!\n",color=:light_red)
        pretty_table(dfs,header=uppercasefirst.(names(dfs));)
        return
    end   
    dx = dfroute(;ofl=ofl);
    rename!(dx,"sim"=>"Basin");
    dfs.Basin = parse.(Int64,dfs.Basin)
    kd  = innerjoin(dfs, dx, on=:Basin)
    if todf
        return kd
    else
        pretty_table(kd,header=uppercasefirst.(names(kd));)
    end
end

"""
returns DataFame with qgk Vals recursively
"""
