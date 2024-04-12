#pdf grep in wsl julia cleaned up from pdfgrep.jl

if length(ARGS) <= 0
	println("need pattern! <pdfgrep>...")
    exit()
end

pattern=ARGS[1];
#using Glob
using PDFIO
using Distributed
# Set the number of workers (threads) for parallel processing
# The number of workers should be equal to or less than the number of CPU cores available.
addprocs(4)  # Set the number of workers

function rglob(prefix::Regex)
    rootdir::String = "."
    results::Vector{String} = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (occursin(prefix,filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end


function pdgrep(pattern)
    pdf_files = rglob(r"pdf$"i)
    for src in pdf_files
        # handle that can be used for subsequence operations on the document.
        doc = pdDocOpen(src)
        # Returns number of pages in the document       
        npage = pdDocGetPageCount(doc)
        println("reading $src ...")
        printstyled("check $src 
        ($npage pages) ->\n",color=:yellow)

        for i in 1:npage
            # handle to the specific page given the number index. 
            #println("check page ")
            #printstyled("$i ",color=:green)
            #printstyled("check $src page $i ...\n",color=:yellow)
            page = pdDocGetPage(doc, i)
            printstyled("$i ",color=:yellow)
            # Use a StringIO object as an IO stream to capture the extracted text.
            io = IOBuffer()
            pdPageExtractText(io, page)
            # Get the extracted text from the StringIO object and convert it to a string.
            extracted_text = String(take!(io))
            # Check if the pattern exists in the extracted text
            if occursin(Regex(pattern,"i"), extracted_text)
                printstyled("Pattern '$pattern' found on page $i of PDF file: $src\n",color=:green)
                #println("Extracted Text:")
                println(extracted_text)
                println("="^50)
            end
        end
        # Close the document handle. 
        pdDocClose(doc)
    end
end

pdgrep(pattern)

#"/mnt/c/Users/chs72fw/Documents/Rhön/Endnote/cs_endnote_lib_2020_v_93.Data/PDF/4283062392"|>cd
#pdgrep("extra")
#pdgrep("entweder")