using Base.Test
import Optics

# write your own tests here
@testset "SphericalLens Tests" begin
    println("Biconcave spherical lens")
    @show R1 = 1
    @show R2 = -1
    @show d = 1
    @show t = 1
    @show nl = 1.5

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    @test xR1 ≈ R1
    @test xR2 ≈ -0.7320508075688772
    @test xV1 ≈ 0
    @test xV2 ≈ 0.2679491924311228

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexThick(R1, R2, t)
    @test xR1 ≈ R1
    @test xR2 ≈ 0
    @test xV1 ≈ 0
    @test xV2 ≈ t

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d)
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
    @test l.xe1 ≈ xe1
    @test l.xe2 ≈ xe2
    @test l.xmin ≈ xV1
    @test l.xmax ≈ xV2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2

    #l = Optics.SphericalLens(R1, R2, d, t, nl)
end
