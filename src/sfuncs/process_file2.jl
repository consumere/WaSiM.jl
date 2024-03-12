# 
function process_file2(input_file::String)
        """
        removes " inside textfile
        """
        lines = readlines(input_file)
        output_file = input_file * ".tmp"
        
        open(output_file, "w") do file
            for line in lines
                line = replace(line, "\"" => "")
                if length(line) > 0
                    println(file, line)
                end
            end
        end
        
        # Rename the temporary file to replace the original file
        mv(output_file, input_file; force=true)
    end

    