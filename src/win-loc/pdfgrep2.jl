#pdf grep in wsl julia 

if length(ARGS) <= 0
	#println("need args! <timeseries file regex!>...")
	println("need pattern! <pdfgrep>...")
    exit()
end

#file=ARGS[1];
pattern=ARGS[1];

using Distributed
#using PDFIO

# Set the number of workers (threads) for parallel processing
# The number of workers should be equal to or less than the number of CPU cores available.
addprocs(4)  # Set the number of workers

# Load the required functions and packages on all workers
@everywhere using PDFIO

function rglob(prefix::Regex)
    rootdir::String = "."
    results::Vector{String} = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (occursin(prefix, filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

# Define the function on all workers using @everywhere
@everywhere function extract_text_from_pdf(pdf_file::AbstractString)
    # handle that can be used for subsequence operations on the document.
    doc = pdDocOpen(pdf_file)

    # Returns number of pages in the document       
    npage = pdDocGetPageCount(doc)

    extracted_texts = Vector{String}(undef, npage)

    for i in 1:npage
        # handle to the specific page given the number index. 
        page = pdDocGetPage(doc, i)
        
        # Use a StringIO object as an IO stream to capture the extracted text.
        io = IOBuffer()
        pdPageExtractText(io, page)
        
        # Get the extracted text from the StringIO object and convert it to a string.
        extracted_text = String(take!(io))
        extracted_texts[i] = extracted_text
    end

    # Close the document handle. 
    # The doc handle should not be used after this call
    pdDocClose(doc)

    return extracted_texts
end

function parallel_pdf_text_extraction(root_dir::AbstractString, pattern::AbstractString)
    # Find all PDF files matching the pattern in the given directory
    pdf_files = rglob(r"pdf$"i)

    # Use pmap to perform parallel processing of PDF text extraction
    results = pmap(pdf_file -> extract_text_from_pdf(pdf_file), pdf_files)

    # Print the results
    for (i, (pdf_file, extracted_texts)) in enumerate(zip(pdf_files, results))
        for (j, extracted_text) in enumerate(extracted_texts)
            if occursin(Regex(pattern, "i"), extracted_text)
                println("Pattern '$pattern' found on page $j of PDF file: $pdf_file")
                println("Extracted Text:")
                println(extracted_text)
                println("="^50)
            end
        end
    end
end


root_dir = "."

parallel_pdf_text_extraction(root_dir, pattern)

