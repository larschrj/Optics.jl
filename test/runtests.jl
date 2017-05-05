using Base.Test
import Optics

# write your own tests here
@testset "SphericalLens Tests" begin
    R1 = 1
    R2 = -1
    d = 1
    t = 1
    nl = 1
    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    l = Optics.SphericalLens(R1, R2, d, nl)
    @test l.R1 ≈ R1
    @test l.R2 ≈ R2
    @test l.t ≈ Optics.sphericalLensMinThick(R1, R2, d)
    @test l.d ≈ d
    @test l.nl ≈ nl
    @test l.xV1 ≈ xV1
    @test l.xV2 ≈ xV2
    @test l.xR1 ≈ xR1
    @test l.xR2 ≈ xR2

    #l = Optics.SphericalLens(R1, R2, d, t, nl)
end
