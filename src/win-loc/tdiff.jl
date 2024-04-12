#!/usr/bin/env julia
# tdiff sheet LWF 66/2008 Falk

using Printf
using Statistics

function looks_like_number(x)
    tryparse(Float64, x) !== nothing
end

function prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
    if looks_like_number(prev)
        ref = (sum / cnt) * fact
        @printf("count: %d | Tdiff in %.2f per %.2f days | sums: %.2f %.2f %.2f |", cnt, ref, sum3, sum2, sum)
        sum, sum2, sum3, cnt = 0, 0, 0, 0
        mean3, cnt = 0, 0
        if ref <= 0
            println("\n")
        elseif ref <= 5
            println("conditions are very moist")
        elseif ref <= 10
            println("conditions are moist")
        elseif ref <= 15
            println("conditions are rather moist")
        elseif ref <= 20
            println("conditions are quite moist")
        elseif ref <= 30
            println("conditions are quite dry")
        elseif ref <= 40
            println("conditions are rather dry")
        elseif ref <= 50
            println("conditions are dry")
        elseif ref <= 70
            println("conditions are very dry")
        elseif ref <= 7e10
            println("conditions are exceptionally dry")
        end
    end
end

function main()
    if length(ARGS) < 3
        println("you don't have enough parameters!\nusage: pmm <infile> <vegetation period> <factor>!")
    else
        f1, fact, fact2 = ARGS[1:3]
        f1 = open(f1)
        if looks_like_number(fact)
            fact = parse(Float64, fact)
            println("\nvegetation period: $fact days on file: $(ARGS[1])")
        else
            fact = 100
            println("fact set to $fact!\n")
        end
        if looks_like_number(fact2)
            fact2 = parse(Float64, fact2)
            println("factor set to $fact2 ...\n")
        else
            fact2 = fact
        end

        sum, sum2, sum3, cnt = 0, 0, 0, 0
        prev = ""
        for line in eachline(f1)
            line = chomp(line)
            Fld = split(line, '\t')
            date = 0
            if looks_like_number(Fld[1])
                date = parse(Int64, Fld[1])
            end
            if date != prev
                prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
                prev = date
            end
            sum += parse(Float64, Fld[end])
            sum2 += parse(Float64, Fld[end - 1])
            sum3 += parse(Float64, Fld[end - 2])
            cnt += 1
        end
        prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
        close(f1)
    end
end

main()


# #tidiff
# using Printf

# function looks_like_number(str)
#     return tryparse(Float64, str) !== nothing
# end

# function prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
#     if looks_like_number(prev)
#         ref = (sum / cnt) * fact
#         @printf("count:%d|Tdiff in %.2f per %d days:%.2f| sums:%.2f, %.2f, %.2f|\n", cnt, prev, fact2, ref, sum3, sum2, sum)
#         sum, sum2, sum3, cnt = 0, 0, 0, 0
#         mean3, cnt = 0, 0
#         if ref <= 0
#             println()
#         elseif ref <= 5
#             println("conditions are very moist")
#         elseif ref <= 10
#             println("conditions are moist")
#         elseif ref <= 15
#             println("conditions are rather moist")
#         elseif ref <= 20
#             println("conditions are quite moist")
#         elseif ref <= 30
#             println("conditions are quite dry")
#         elseif ref <= 40
#             println("conditions are rather dry")
#         elseif ref <= 50
#             println("conditions are dry")
#         elseif ref <= 70
#             println("conditions are very dry")
#         elseif ref <= 7e10
#             println("conditions are exceptionally dry")
#         end
#     end
#     return sum, sum2, sum3, cnt
# end

# function process_file(f1, fact, fact2)
#     sum, sum2, sum3, cnt = 0, 0, 0, 0
#     mean3, cnt = 0, 0
#     prev = ""
#     date = 0

#     open(f1) do file
#         for line in eachline(file)
#             line = strip(line)
#             if !startswith(line, '#') || !startswith(line, 'Y') && line != ""
#                 fields = split(line, '\t')
#                 if looks_like_number(fields[1])
#                     date = parse(Float64, fields[1])
#                 end
#                 if date != prev
#                     sum, sum2, sum3, cnt = prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
#                     prev = date
#                 end
#                 sum += parse(Float64, fields[end])
#                 sum2 += parse(Float64, fields[end-1])
#                 sum3 += parse(Float64, fields[end-2])
#                 cnt += 1
#             end
#         end
#         sum, sum2, sum3, cnt = prtMean(prev, sum, sum2, sum3, cnt, fact, fact2)
#     end
# end

# function main()
#     args = ARGS

#     if length(args) < 3
#         println("you don't have enough parameters!\nusage: julia script.jl <infile> <vegetation period> <factor>!")
#     else
#         f1 = args[1]
#         fact = parse(Float64, args[2])
#         fact2 = parse(Float64, args[3])

#         if looks_like_number(args[2])
#             println("vegetation period: ", fact, " days on file: ", f1)
#         else
#             fact = 100
#             println("fact set to ", fact, "!")
#         end

#         if looks_like_number(args[3])
#             println("factor set to ", fact2, " ...")
#         else
#             fact2 = fact
#         end

#         process_file(f1, fact, fact2)
#     end
# end

# main()

# #pt=/mnt/c/Users/Public/Documents/Python_Scripts/julia/tdiff.jl
# #"C:\Users\Public\Documents\Python_Scripts\julia\tdiff.jl"
