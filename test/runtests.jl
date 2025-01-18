# %%
using WaSiM
using Test
# %%
ftest(x) = @test isone(x)

@testset verbose = true "arithmetrics" begin
    # Write your tests here.
    θ = 2/3*π
    @test sin(-θ) ≈ -sin(θ)
    ftest(1)
    #@test WaSiM.llf()
end

# @testset verbose = true "my ls test" begin
#    ls()
# end

# @testset verbose = true "my du test" begin
#    du() # also works
# end

# @testset verbose = true "my ls test" begin
   # ls();
   # du();
   # #sf("func")
# end

#using JET
#report_package("WaSiM")