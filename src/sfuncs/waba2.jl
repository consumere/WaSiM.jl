# 
function waba2()
        begin
            af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
            if length(af) <= 2 || any(ismissing.(af))
                error("match failed \n ... abort.\nno special output files present!")
                return #exit(86)
            end
        
            println("calculating yearly water balance of WaSiM special output data..\n
            bwvr = rain + snow + uprs - perc - qb - qd - qi - etr_ - ei_ - etrs_\n")
        
            re = Regex("preci|snow_storage_tota|Capi|Perc|baseflo|directflo|interflo|real_evap|real_tran|interception_evaporatio|snow_evaporatio","i")
            my = filter(x -> occursin(re, x), af)
            printstyled("loading...\n $my\n",color=:green)
        
            if (length(my) .!= 11)==true
                lng=length(my)
                printstyled("found only $lng files...\n $my\n",color=:yellow)
                error("\nfiles are missing!\ncheck so files...")
            end
        
            rain = filter(x -> occursin("precip", x), af)|>only|>so_read
            snow = filter(x -> occursin("snow_storage_total", x), af)|>only|>so_read
            uprs = filter(x -> occursin("Capil", x), af)|>only|>so_read
            perc = filter(x -> occursin("Perco", x), af)|>only|>so_read
            qb = filter(x -> occursin("baseflow", x), af)|>only|>so_read
            qd = filter(x -> occursin("directflow", x), af)|>only|>so_read
            qifl = filter(x -> occursin("interflow", x), af)|>only|>so_read
            etr = filter(x -> occursin("real_evapo", x), af)|>only|>so_read
            etrans = filter(x -> occursin("real_trans", x), af)|>only|>so_read
            ei = filter(x -> occursin("interception_evaporation", x), af)|>only|>so_read
            etrs = filter(x -> occursin("snow_evaporation",x), af)|>only|>so_read
        
            l = [rain, snow, uprs, perc, qb, qd, qifl, etr, etrans, ei, etrs]
            #typeof(l)
            nm = [names(l[i])[1] for i in 1:size(l, 1)]
            #same:
            #map(x->names(x)[1],l)
            println("loaded dataframes:\n$nm")
        
            if (length(l) .!= 11)==true
                error("files are missing!\ncheck so files...")
            end
        
            d = mall(l)
            xd=copy(d)
            writewa("waba-input.wa",xd)
        
            #dyr = yrsum(d)
        
            pos = d[!,Cols(r"date|^(prec)|^(snow_stora)|^(Capi)")]
            pos = yrsum(pos)
            # calculate the sum of each row
            psum = DataFrame(
                possums = [sum(eachrow(pos)[i]) for i in 1:size(pos, 1)],
                year=pos[!,:year]
            )
        
        
            neg = d[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)"))]
            neg = sum.(yrsum(neg))
        
        
            nsum = DataFrame(
                negsums = [sum(eachrow(neg)[i]) for i in 1:size(neg, 1)],
                year=neg[!,:year]
            )
        
        
            bw = innerjoin(psum, nsum, on=:year)
            bw.bw = bw[!,:possums] .- bw[!,:negsums]
        end
        ti="water-balance of "*basename(pwd())
        #theme(:vibrant)
        #theme(:wong)
        #Plots.theme(:dao)       #latex fonts.
        #theme(:mute)
        #theme(:sand)
        ann = map(x->string.(round(x;sigdigits=3))*" [mm]",bw.bw)
        fact=.60
        plotsize = (1600*fact,800*fact)
        p1 = @df bw Plots.plot(
            :year,:bw,
            #annotations =(bw.year, bw.bw, ann, :top),
            #annotations = (bw.year,bw.bw,(ann,10,:center,:top,:black)),
            #annotations = (bw.year, bw.bw, ann, 8, :left, :top, :black),
            legend = false, 
            seriestype=:bar,
            xticks = bw.year,
            xtickfont = 12,
            xlabel = "",
            ylabel = "[mm]",
            ytickfont = 12,
            title = ti,
            fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"),
            size=plotsize,
            #xrotation = 60);
            left_margin = 10mm,
            bottom_margin = 10mm, 
            #bottom_margin = 10px, 
            xrotation = 45)
            
        #Plots.annotate!(bw.year,bw.bw,(ann,10,:left,:top,:black))
                    #text("hey", 14, :left, :top, :green)
        for i in 1:length(bw.year)
            Plots.annotate!(bw.year[i],bw.bw[i],(ann[i],10,:center,:top,:black))
            #println(ann[i]*" added")
        end
        #display(p1)
        #return(bw)
        return p1
    end

    # ##very nice meta programming... move to smallfuncs.jl
    # macro vv(s) vgjl(s);end
    # #@vv "unter"
    # macro vpy(s) vgpy(s);end
    # #@vpy "climate"
    # macro vr(s) vgr(s);end
    # #@vr "climate"
    # macro vct(s) vgctl(s);end
    # #@vct "das ist"
    # macro rg(s) rglob(s);end
    # macro glb(s) glob(s);end
    # macro gl(s) glob(s)|>first;end
    # #fastplot
    # macro fp(s) dfp(Regex(s));end
    # macro flog(s) dfl(Regex(s));end
    # macro ncrm() ncrem=src_path*"/ncremover.jl";include(ncrem);end
    # macro rasrm() remer=src_path*"/raster_remover.jl";include(remer);end
    # macro nco(s) first(nconly(s));end
    # macro dfo(s) first(dfonly(s));end
    # macro wajs() pt=src_path*"/wajs.jl";include(pt);end
    # macro cmk() pt=src_path*"/cairomakie.jl";include(pt);end
    # macro rcall() pt=src_path*"/rcall.jl";include(pt);end
    # macro hd(df) df[1:4,:];end
    # macro pyjl() pt=src_path*"/pyjl.jl";include(pt);end

    