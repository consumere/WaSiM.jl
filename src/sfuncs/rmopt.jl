# 
function rmopt()
        # Get the list of files
        files = rglob("")

        #!occursin(r"pl|sh|csv|html|xml|fzt|ftz|log|ini|^wq|yrly|nc|tif|jpeg|png|svg|txt", file)
        # py|R|ftz_0|tex
        # Define the file extensions to keep
        keep_exts = [".ipynb", ".py", ".R", ".Rmd", ".log",
            ".tif", ".jpeg", ".png", ".svg",
            ".cpg", ".shx", ".dbf", ".prj", ".shp", ".tex", 
            ".csv", 
            ".html", ".ftz", ".ftz_0", ".txt", 
            ".list", ".nc", ".xml", ".sh", ".grd", ".yrly"]
            
        files = filter(file -> 
            !occursin(r"(wq_|pl$|sh$|fzt|ftz|log$|ini|otherdata|intern)", file)
            , files)
        
            # Process each file
        for file in files
            # Check if the file extension is in the list of extensions to keep
            if !(splitext(file)[2] in keep_exts)
                try
                # Read the file as a DataFrame
                println("Processing $file ...")
                ms = ["-9999", "lin", "log", "--"]
                #df = CSV.read(file, DataFrame)
                df = CSV.File(file; 
                    delim="\t", 
                    header=1,
                    silencewarnings=true, 
                    normalizenames=false, 
                    missingstring=ms, 
                    types=Float64) |> DataFrame
                dropmissing!(df,ncol(df))
                # dropmissing!(df,4)
                # Calculate the sums of the 5th column and the last column
                y = sum(df[:, 5])
                x = sum(df[:, end])

                # Check the conditions
                if x <= (nrow(df) - 5) * -9999 && y <= (nrow(df) - 5) * -9999 || x == 0 && y == 0
                    # Remove the file
                    rm(file)
                    printstyled("$file removed ...\n", color=:red)
                end
                catch e
                    println("Error processing $file: $e")
                end
            end
        end
    end

    