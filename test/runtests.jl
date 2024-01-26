using WaSiM
using Test

f(x) = @test isone(x)
@testset f(1)

@testset verbose = true "WaSiM.jl" begin
    # Write your tests here.
    θ = 2/3*π
    @test sin(-θ) ≈ -sin(θ)
    #@test ct()
    #@test WaSiM.llf()
end
