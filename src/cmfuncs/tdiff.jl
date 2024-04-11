# 
function tdiff()
        npot = glob("so_pot_trans")|>first|>dfread
        nreal =  glob("so_real_trans*")|>first|>dfread
        td = innerjoin(npot,nreal,on=:date)
        td = reorder_df(td)
        #Tdiff = hcat(td,td[!,1] .- td[!,2])
        td.Tdiff = td[!,1] .- td[!,2]
        te(td)
        return(reorder_df(td))
        #wawrite(td,"tdiff-jl.txt")
    end


    