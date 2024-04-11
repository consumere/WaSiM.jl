# 
function odfr(x::Regex)
        """
        --- reader with fewer constrains ---
        with |> dropmissing on 2nd col
        df[!,Cols(r"^Col|date")] |>dfp  
        """
        x=globdf(x)|>first
        println("reading $x ...")
        ms=["-9999","lin","log","--","A-z"]
        df = CSV.read(x,    
        DataFrame,    
        missingstring=ms,
        #delim="\t",    
        types = Float64,
        normalizenames=true,
        drop=(i, nm) -> i == 4) 
        dropmissing!(df , 2) #2nd column
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,4:end]
        df=df[!,Not(Cols(r"^Column"))] #drops names starting with Column, usually all missing
        #renamer
        for x in 1:size(df,2)-2
            rename!(df,x=>"C"*names(df)[x])
        end
        DataFrames.metadata!(df, "filename", x, style=:note);
    end

    