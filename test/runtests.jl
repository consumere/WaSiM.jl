# %%
using WaSiM
using Test
# %%
ftest(x) = @test isone(x)
@testset ftest(1)

@testset verbose = true "arithmetrics" begin
    # Write your tests here.
    θ = 2/3*π
    @test sin(-θ) ≈ -sin(θ)
    #@test ct()
    #@test WaSiM.llf()
end


@testset verbose = true "my ls test" begin
   ls()
end

@testset verbose = true "my du test" begin
   du() # also works
end
