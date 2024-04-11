# 
function getmoduleorder(file1::AbstractString)
        # Run the awk command and capture the STDOUT
        result = read(`awk '$0 ~ /^[[]/{c=1} c&&c--' $file1`, String)
    
        # Split the result into lines and strip \r and \t characters
        lines = replace(split(result, "\n"), r"[\r\t]" => "")
    
        # Extract the names inside square brackets using regular expressions
        pattern = r"\[(.*?)\]"
        names_inside_brackets = [match(pattern, line) === nothing ? "" : match(pattern, line).captures[1] for line in lines]
    
        # Create a DataFrame with the extracted names inside the square brackets
        df = DataFrame(Line = names_inside_brackets)
    
        return df
    end

    