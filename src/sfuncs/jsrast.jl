# 
function jsrast(R::String)
        rpath="/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions"
        run(`wsl -d Ubuntu-18.04 -e Rscript "$rpath/plotly_rast.R" $R`)
    end
    
    