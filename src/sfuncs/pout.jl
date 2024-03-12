# 
function pout(infile;sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/",ofl = "route.txt")
        routeg(infile, ofl)
        sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
        specfile=joinpath(sfpt,sfn)
        obs = readdf(specfile)
        df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
        rename!(df,1=>"sim",2=>"obs",3=>"name")
        df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
        df.name=map(x->replace(x,r"_>.*" => ""),df.name)
        sort!(df, :sim)
        sim = r"qges"|>glob|>first|>waread
        @info "innerjoin has to be sim -> obs!"
        outd = []
        for i in eachrow(df)
            println(i[1],"->",i[3])
            try
                dm =innerjoin(
                    sim[!,Cols("C"*string(i[1]),end)],
                    obs[!,Cols(Regex(i[3],"i"),end)],    
                    on=:date)
                onam = i[3]*"-qoutjl"
                wawrite(dm,onam)
                println("$onam saved!")
                DataFrames.metadata!(dm, "filename", onam, style=:note);
                push!(outd,dm)
                println(names(dm)," on $onam pushed!")
            catch
                onam = i[3]*"-qoutjl"
                @warn "merge is empty for $onam ! ..."
                continue
            end
        end
        vz = filter(xx->xx[2]==5,cntcolv("outjl"))
        if length(vz)>0
            @warn "some files are not properly merged..."
        end
        return outd
    end

    """
    dfroute(;ofl="route.txt")
    reads from routeg(infile, ofl) and returns a DataFrame with the following columns:
        - sim: simulated flow
        - obs: observed flow
        - name: name of the station
    """
    