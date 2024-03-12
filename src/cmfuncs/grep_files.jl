# 
function grep_files(pattern,file_paths, context)
        for file_path in file_paths
            grep_with_context(pattern, file_path, context)
        end
    end

    