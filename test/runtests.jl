using WaSiM
using Test

"""
    ftest(x)

Helper function to test if input equals one.
"""
ftest(x) = @test isone(x)

@testset verbose=true "Basic Arithmetic Tests" begin
    θ = 2/3 * π
    @test sin(-θ) ≈ -sin(θ)
    ftest(1)
end

@testset verbose=true "File System Operations" begin
    # Test directory listing and disk usage
    ls()
    du()
end

# Uncommented section for future implementation
#=
@testset verbose=true "Additional Tests" begin
    # Placeholder for WaSiM specific tests
    # @test WaSiM.llf()
    # sf("func")
end

# Static analysis
# using JET
# report_package("WaSiM")
=#
