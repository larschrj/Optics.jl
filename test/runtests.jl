using Base.Test
import Optics

# write your own tests here
@testset "SphericalLens Tests" begin
    println("Biconcave spherical lens")
    @show R1 = 2
    @show R2 = -2
    @show d = 1
    @show t = 1
    @show nl = 1.5

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    @test xR1 ≈ 2
    @test xR2 ≈ -1.872983346207417
    @test xV1 ≈ 0
    @test xV2 ≈ 0.12701665379258298

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d)
    @test xe1 ≈ 0.06350832689629149
    @test xe2 ≈ 0.06350832689629149

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

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexThick(R1, R2, t)
    @test xR1 ≈ 2
    @test xR2 ≈ -1
    @test xV1 ≈ 0
    @test xV2 ≈ 1

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d, t)
    @test xe1 ≈ 0.06350832689629149
    @test xe2 ≈ 0.9364916731037085

    l = Optics.SphericalLens(R1, R2, t, d, nl)
    @test l.R1 ≈ R1
    @test l.R2 ≈ R2
    @test l.t ≈ t
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



    println("R1 > 0, R2 > 0, spherical lens")
    @show R1 = 2
    @show R2 = 2
    @show d = 1
    @show t = 1
    @show nl = 1.5

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    @test xR1 ≈ 2.0
    @test xR2 ≈ 2.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 0.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d)
    @test xe1 ≈ 0.06350832689629149
    @test xe2 ≈ 0.06350832689629149

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
    @test l.xmax ≈ xe2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexThick(R1, R2, t)
    @test xR1 ≈ 2.0
    @test xR2 ≈ 3.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 1.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d, t)
    @test xe1 ≈ 0.06350832689629149
    @test xe2 ≈ 1.06350832689629149

    l = Optics.SphericalLens(R1, R2, t, d, nl)
    @test l.R1 ≈ R1
    @test l.R2 ≈ R2
    @test l.t ≈ t
    @test l.d ≈ d
    @test l.nl ≈ nl
    @test l.xV1 ≈ xV1
    @test l.xV2 ≈ xV2
    @test l.xR1 ≈ xR1
    @test l.xR2 ≈ xR2
    @test l.xe1 ≈ xe1
    @test l.xe2 ≈ xe2
    @test l.xmin ≈ xV1
    @test l.xmax ≈ xe2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2



    println("R1 < 0, R2 > 0, spherical lens")
    @show R1 = -2
    @show R2 = 2
    @show d = 1
    @show t = 1
    @show nl = 1.5

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    @test xR1 ≈ -2.0
    @test xR2 ≈ 2.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 0.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d)
    @test xe1 ≈ -0.06350832689629149
    @test xe2 ≈ 0.06350832689629149

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
    @test l.xmin ≈ xe1
    @test l.xmax ≈ xe2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexThick(R1, R2, t)
    @test xR1 ≈ -2.0
    @test xR2 ≈ 3.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 1.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d, t)
    @test xe1 ≈ -0.06350832689629149
    @test xe2 ≈ 1.06350832689629149

    l = Optics.SphericalLens(R1, R2, t, d, nl)
    @test l.R1 ≈ R1
    @test l.R2 ≈ R2
    @test l.t ≈ t
    @test l.d ≈ d
    @test l.nl ≈ nl
    @test l.xV1 ≈ xV1
    @test l.xV2 ≈ xV2
    @test l.xR1 ≈ xR1
    @test l.xR2 ≈ xR2
    @test l.xe1 ≈ xe1
    @test l.xe2 ≈ xe2
    @test l.xmin ≈ xe1
    @test l.xmax ≈ xe2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2



    println("R1 < 0, R2 < 0, spherical lens")
    @show R1 = -2
    @show R2 = -2
    @show d = 1
    @show t = 1
    @show nl = 1.5

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexDiam(R1, R2, d)
    @test xR1 ≈ -2.0
    @test xR2 ≈ -2.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 0.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d)
    @test xe1 ≈ -0.06350832689629149
    @test xe2 ≈ -0.06350832689629149

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
    @test l.xmin ≈ xe1
    @test l.xmax ≈ xV2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2

    (xR1, xR2, xV1, xV2) = Optics.sphericalLensVertexThick(R1, R2, t)
    @test xR1 ≈ -2.0
    @test xR2 ≈ -1.0
    @test xV1 ≈ 0.0
    @test xV2 ≈ 1.0

    (xe1, xe2) = Optics.sphericalLensxEdge(R1, R2, d, t)
    @test xe1 ≈ -0.06350832689629149
    @test xe2 ≈ 0.9364916731037085

    l = Optics.SphericalLens(R1, R2, t, d, nl)
    @test l.R1 ≈ R1
    @test l.R2 ≈ R2
    @test l.t ≈ t
    @test l.d ≈ d
    @test l.nl ≈ nl
    @test l.xV1 ≈ xV1
    @test l.xV2 ≈ xV2
    @test l.xR1 ≈ xR1
    @test l.xR2 ≈ xR2
    @test l.xe1 ≈ xe1
    @test l.xe2 ≈ xe2
    @test l.xmin ≈ xe1
    @test l.xmax ≈ xV2
    @test l.xmid ≈ (xV1 + xV2)/2
    @test l.yR ≈ 0
    @test l.ymin ≈ -d/2
    @test l.ymax ≈ d/2
end
