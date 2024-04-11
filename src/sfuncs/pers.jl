# 
function pers(x::Union{String,Regex}=r"qoutjl$";c1=5,c2=4,f2=nothing,tofile=false)
        bsp = "/mnt/c/Users/Public/Documents/Python_Scripts/bashscripts/pearson_wsl.pl"
        dirpt = last(split(dirname(pwd()),"\\"))
        res = "pers_"*dirpt*".txt"
        if f2 != nothing
            f1 = isa(x,String) ? x : first(x)
            if tofile
                #res = readchomp(pipeline(cmd))
                #run(`wsl sh -c "(echo bla) | tee z.txt" `)
                run(`wsl sh -c "$bsp $f1 $f2 $c1 $c2 | tee $res" `)
            else
                run(`wsl $bsp $f1 $f2 $c1 $c2`)
            end
        else   
        files = glob(x)
            for x in files
                if tofile
                    res = "pers_"*dirpt*".txt"
                    res = replace(res,"pers"=>x)
                    #run(`wsl sh -c "($cmd) | tee $res" `)
                    run(`wsl sh -c "$bsp $x $x $c1 $c2 | tee $res" `)
                else
                    run(`wsl $bsp $x $x $c1 $c2`)
                end
                
            end
        end
    end

    """
    takes Vector{DataFrame} from:
    dataframes = map(byear,outd)
    plot_grouped_metrics(dataframes::Vector{DataFrame};col=:ve,all=false,kw...)
    """
    