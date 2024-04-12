# if length(ARGS) == 0
# 	println("need args! <file>...")
#     exit()
# end
#using Printf

this=pwd()
println("performing rmdub on $this...")

using SHA

"""
recursively removes duplicates 
uses SHA
"""
function rmdub(;directory=pwd())
    
    # Create a dictionary to store file hashes as keys and file paths as values
    #hash_dict = Dict{String, String}()
    hash_dict = Dict{Vector{UInt8}, String}()

    # Get a list of all files in the directory
    #files = readdir(directory)
    for (root, dirs, files) in walkdir(directory)
        for file in files
            filepath = joinpath(directory, file)
            if isfile(filepath)
                # Calculate the SHA-256 hash of the file's contents
                println("Hashing file: $filepath")
                io = open(filepath, "r")
                filehash = sha256(io)
                close(io)
                # If the hash is not already in the dictionary, add it
                if !haskey(hash_dict, filehash)
                    hash_dict[filehash] = filepath
                else
                    # If a file with the same hash is found, delete it
                    printstyled("Deleting duplicate: $filepath \n",color=:red)
                    rm(filepath)
                end
            end
        end
    end
end

rmdub()

#in wsl:
#julia --startup-file=no $(wslpath "C:\Users\Public\Documents\Python_Scripts\julia\rmdub.jl")
#julia --startup-file=no --optimize=0 --threads 8  "C:\Users\Public\Documents\Python_Scripts\julia\rmdub.jl"

