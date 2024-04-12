#pdf grep in wsl julia 

if length(ARGS) <= 0
	#println("need args! <timeseries file regex!>...")
	println("need pattern! <pdfgrep>...")
    exit()
end

#file=ARGS[1];
pattern=ARGS[1];

#using Glob
using PDFIO

# """
# ​```
#     getPDFText(src, out) -> Dict 
# ​```
# - src - Input PDF file path from where text is to be extracted
# - out - Output TXT file path where the output will be written
# return - A dictionary containing metadata of the document
# """
# function getPDFText(src, out)
#     # handle that can be used for subsequence operations on the document.
#     doc = pdDocOpen(src)
    
#     # Metadata extracted from the PDF document. 
#     # This value is retained and returned as the return from the function. 
#     docinfo = pdDocGetInfo(doc) 
#     open(out, "w") do io
    
#         # Returns number of pages in the document       
#         npage = pdDocGetPageCount(doc)

#         for i=1:npage
        
#             # handle to the specific page given the number index. 
#             page = pdDocGetPage(doc, i)
            
#             # Extract text from the page and write it to the output file.
#             pdPageExtractText(io, page)

#         end
#     end
#     # Close the document handle. 
#     # The doc handle should not be used after this call
#     pdDocClose(doc)
#     return docinfo
# end

# #getPDFText(pdf_file,di)

# function rglob(prefix::AbstractString)
#     rootdir::String = "."
#     results::Vector{String} = []
#     #results = [] #Vector{Any}
#     for (looproot, dirs, filenames) in walkdir(rootdir)
#         for filename in filenames
#             #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#             if (occursin(Regex(prefix,"i"),filename))
#                 push!(results, joinpath(looproot, filename)) 
#             end
#         end
#     end
#     return results
# end

# function search_and_extract_text(root_dir::AbstractString, pattern::AbstractString)
#     # Find all PDF files matching the pattern in the given directory
#     #pdf_files = Glob.glob("$(root_dir)/**/*.pdf", case_insensitive=true)
#     #"/mnt/c/Users/chs72fw/Documents/Rhön/Endnote/cs_endnote_lib_2020_v_93.Data/PDF/4283062392"|>cd
#     pattern = "Lei"
#     pdf_files = rglob("pdf")

#     for pdf_file in pdf_files
#         # Extract text from the PDF file
#         pdf_file = first(pdf_files)
#         #text = read(pdf_file, PDFText)
#         doc = PDFIO.pdDocOpen(pdf_file)
#         npage = pdDocGetPageCount(doc)
#         page = pdDocGetPage(doc, 1)
        
#         pdPageExtractText(page)

#         # Check if the pattern exists in the extracted text
#         if occursin(pattern, )
#             println("Pattern '$pattern' found in PDF file: $pdf_file")
#             # Optionally, you can process or display the extracted text here.
#             # For example, print the first 100 characters of the text:
#             println("Extracted Text:")
#             println(text[1:min(end, 100)])
#             println("="^50)
#         end
#         PDFIO.pdDocClose(pdf_file)
#     end
# end

# # Example usage:
# root_dir = pwd()
# search_and_extract_text(root_dir, "geodata")


# using PDFIO

# function searchAndPrintPDFText(src, pattern)
#     # handle that can be used for subsequence operations on the document.
#     doc = pdDocOpen(src)

#     # Returns number of pages in the document       
#     npage = pdDocGetPageCount(doc)

#     for i in 1:npage
#         # handle to the specific page given the number index. 
#         page = pdDocGetPage(doc, i)
        
#         # Use a StringIO object as an IO stream to capture the extracted text.
#         io = IOBuffer()
#         pdPageExtractText(io, page)
        
#         # Get the extracted text from the StringIO object and convert it to a string.
#         extracted_text = String(take!(io))
        
#         # Check if the pattern exists in the extracted text
#         if occursin(pattern, extracted_text)
#             println("Pattern '$pattern' found on page $i of PDF file: $src")
#             println("Extracted Text:")
#             println(extracted_text)
#             println("="^10)
#         end
#     end

#     # Close the document handle. 
#     # The doc handle should not be used after this call
#     pdDocClose(doc)
# end

# # Example usage:
# #pdf_file #= "path/to/your/pdf_file.pdf"
# search_pattern = "Leit"
# searchAndPrintPDFText(pdf_file, search_pattern)

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
        pdDocClose(doc)
    end

    # Close the document handle. 
    # The doc handle should not be used after this call
    #pdDocClose(doc)
end

pdgrep(pattern)



#"/mnt/c/Users/chs72fw/Documents/Rhön/Endnote/cs_endnote_lib_2020_v_93.Data/PDF/4283062392"|>cd
#pdgrep("extra")