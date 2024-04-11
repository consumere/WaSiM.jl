# 
function process_folders_and_subfolders(root_dir::AbstractString, log_file::AbstractString = "output.log")
        open(log_file, "w") do io
            for (looproot, dirs, filenames) in walkdir(root_dir)
                for filename in filenames
                    if occursin(r"xml", filename)
                        entry = joinpath(looproot, filename)
                        result = extract_duration_from_xml(entry)
    
                        if result !== nothing
                            folder_name = dirname(entry)
                            println(io, "XML File: $entry")
                            hours, minutes = result
                            println(io, "Duration: $hours hours and $minutes minutes")
                            println(io, "="^50)
                        end
                    end
                end
            end
        end
    end

    