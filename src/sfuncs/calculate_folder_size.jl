# 
function calculate_folder_size(directory)
        size = 0
        count = 0
        for (root, dirs, files) in walkdir(directory)
            for file in files
                size += stat(joinpath(root, file)).size
                count += 1
            end
        end
        return size, count
    end

    