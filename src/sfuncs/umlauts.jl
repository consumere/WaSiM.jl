#  which python 
function umlauts(input_file::AbstractString, output_file::AbstractString)
    #     # Read input file
    #     data = readdlm(input_file, String)

    #     # Apply string replacements
    #     for i in 1:size(data, 1)
    #         # data[i] = replace(data[i], r"ß" => "ss")
    #         # data[i] = replace(data[i], r"\/" => "_")
    #         # data[i] = replace(data[i], r"_" => "-")
    #         # data[i] = replace(data[i], r",," => "")
    #         # data[i] = replace(data[i], r"\xc4" => "Ae")
    #         # data[i] = replace(data[i], r"\xd6" => "Oe")
    #         # data[i] = replace(data[i], r"\xdc" => "Ue")
    #         # data[i] = replace(data[i], r"\xe4" => "ae")
    #         # data[i] = replace(data[i], r"\xf6" => "oe")
    #         # data[i] = replace(data[i], r"\xfc" => "ue")
    #         # data[i] = replace(data[i], r"\xdf" => "ss")
    #     end

    #     # Write modified data back to the file
    #     output_file = open(output_file, "w")
    #     for i in 1:size(data, 1)
    #         println(output_file, data[i])
    #     end
    #     close(output_file)
    # end
    
    #umlauts(input_file, output_file)

    macro bash_str(s) open(`bash`,"w",stdout) do io; print(io, s); end;end
    #bash""" which python """

    macro pwrs_str(s) open(`powershell -noprofile`,"w",stdout) do io; print(io, s); end;end
    #pwrs""" which python """
    #pwrs""" pwd """
    