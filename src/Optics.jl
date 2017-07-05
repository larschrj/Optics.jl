"""
Geometric (Gaussian) optics module
* fo: primary (object) focal length from the primary principal plane
* fi: secondary (image) focal length from the secondary principal plane
* nm1: Left medium (surface 1 side) index of refraction
* nm2: Right medium (surface 2 side) index of refraction
"""
module Optics
 
using Core
using Base
using PyCall
using PyPlot
using Colors

import Base.sort
import Base.sort!
import Base.copy
#import Base.reverse

export Aperture,
CircleAperture,
Optic,
Lens,
SphericalLens,
OpticsSystem,
OpticsSystemMatrix,
RGB2XYZ,
RGB2sRGB,
XYZ2xyz,
circleNormAngle,
draw,
lineCircleIntersect,
lineCircleNormAngle,
lineLineIntersect,
matrixRayTrace,
move,
move!,
overlap,
pyplot,
rayLensAngle,
rayLensIntersect,
rayLensRefract,
rayMatrix,
rayTrace,
refractMatrix,
snell,
spectrum2RGB,
spectrum2RGBStilesBurch,
spectrum2XYZ,
spectrum2sRGB,
sphericalLensMaxDiam,
sphericalLensMinThick,
sphericalLensR1,
sphericalLensR2,
sphericalLensRfi,
sphericalLensRfo,
sphericalLensVertexDiam,
sphericalLensVertexThick,
sphericalLensfi,
sphericalLensfi2fo,
sphericalLensfo,
sphericalLensfo2fi,
sphericalLenssi,
sphericalLensso,
sphericalLensx,
sphericalLensxH1,
sphericalLensxH2,
sphericalLensxV1xH1,
sphericalLensxV2xH2,
transferMatrix



abstract Optic <: Any
abstract Lens <: Optic
abstract Aperture <: Optic



"""
Spherical lens
* R1::Real Radius of curvature of the first surface of the lens, positive with the center to the right of the first surface
* R2::Real Radius of curvature of the second surface of the lens, positive with the center to the right of the second surface
* t::Real Thickness on optical axis
* d::Real Diameter
* nl::Number Lens index of refraction
* xV1::Real Vertex 1 x coordinate
* xV2::Real Vertex 2 x coordinate
* xR1::Real Surface 1 center of curvature x coordinate
* xR2::Real Surface 2 center of curvature x coordinate
* xe1::Real Surface 1 top and bottom edge x coordinate
* xe2::Real Surface 2 top and bottom edge x coordinate
* xmin::Real Surface 1 minimum x coordinate
* xmax::Real Surface 2 maximum x coordinate
* yR::Real y coordinate of the radius curvature centers (yR1 = yR2 = yR)
* ymin::Real Minimum y coordinate
* ymax::Real Maximum y coordinate
"""
type SphericalLens <: Lens
    R1::Real # surface 1 radius of curvature
    R2::Real # surface 2 radius of curvature
    t::Real # thickness on optical axis
    d::Real # diameter
    nl::Number # index of refraction
    xV1::Real # surface 1 vertex x coord 
    xV2::Real # surface 2 vertex x coord
    xR1::Real # surface 1 center of curvature x coord
    xR2::Real # surface 2 center of curvature x coord 
    xe1::Real # surface 1 edge x coord
    xe2::Real # surface 2 edge x coord
    xmin::Real # surface 1 min x coord
    xmax::Real # surface 2 max x coord
    xmid::Real # lens middle x coord
    yR::Real # y coord of the radius curvature centers (yR1 = yR2 = yR)
    ymin::Real # min y coord
    ymax::Real # max y coord 
end


"""
    SphericalLens(R1::Real, R2::Real, d::Real, nl::Real)::SphericalLens
"""
function SphericalLens(R1::Real, R2::Real, d::Real, nl::Real)::SphericalLens
    (xR1, xR2, xV1, xV2) = sphericalLensVertexDiam(R1, R2, d);
    t = abs(xV1 - xV2);
    (xe1, xe2) = sphericalLensxEdge(R1, R2, d)
    xmin = minimum((xV1, xe1));
    xmax = maximum((xV2, xe2));
    xmid = (xV1 + xV2)/2
    yR = 0;
    ymin = yR - d/2;
    ymax = yR + d/2;
    return SphericalLens(R1, R2, t, d, nl, xV1, xV2, xR1, xR2, xe1, xe2, xmin, xmax, xmid, yR, ymin, ymax)
end


"""
    SphericalLens(R1::Real, R2::Real, t::Real, d::Real, nl::Real)::SphericalLens
"""
function SphericalLens(R1::Real, R2::Real, t::Real, d::Real, nl::Real)::SphericalLens
    dmax = sphericalLensMaxDiam(R1, R2, t);
    @assert d < dmax "d > dmax = $dmax";
    (xR1, xR2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t);
    (xe1, xe2) = sphericalLensxEdge(R1, R2, d, t)
    xmin = minimum((xV1, xe1));
    xmax = maximum((xV2, xe2));
    xmid = (xV1 + xV2)/2
    yR = 0;
    ymin = yR - d/2;
    ymax = yR + d/2;
    return SphericalLens(R1, R2, t, d, nl, xV1, xV2, xR1, xR2, xe1, xe2, xmin, xmax, xmid, yR, ymin, ymax)
end


"""
    SphericalLens(l::SphericalLens)::SphericalLens

Create a copy of a spherical lens l
"""
function SphericalLens(l::SphericalLens)::SphericalLens
    return SphericalLens(l.R1, l.R2, l.t, l.d, l.nl)
end 


"""
    copy(l::SphericalLens)::SphericalLens
"""
function copy(l::SphericalLens)::SphericalLens
    lensCopy = SphericalLens(l)
    return lensCopy
end


"""
Circular aperture
* d1::Real Diameter of opening 
* d2::Real Diameter of aperture (d2 > d1)
* x::Real x coordinate of aperture center
* y::Real y coordinate of aperture center (center of the opening)
"""
type CircleAperture <:Aperture
    d1::Real
    d2::Real
    x::Real
    y::Real
    xmin::Real
    xmax::Real
    xmid::Real
    ymin::Real
    ymax::Real

    """
        CircleAperture(d1::Real, d2::Real, x::Real, y::Real)
    """
    function CircleAperture(d1::Real, d2::Real, x::Real, y::Real)
        @assert d1 < d2 "d1 < d2"
        new(d1, d2, x, y, x, x, x, y - d2/2, y + d2/2)
    end

    """
        CircleAperture(a::CircleAperture)
    """
    function CircleAperture(a::CircleAperture)
        new(a.d1, a.d2, a.x, a.y, a.x, a.x, a.x, a.ymin, a.ymax)
    end
end 


"""
    overlap(l1::SphericalLens, l2::SphericalLens)
"""
function overlap(l1::SphericalLens, l2::SphericalLens) 
    if l1.xV1 > l2.xV1
        temp = l1
        l1 = l2
        l2 = temp
    end
    
    if !isinf(l1.R2) & !isinf(l2.R1) 
        if (l1.xR2 - l2.xR1) == 0
            if l1.R2 < l2.R1
                return true
            else
                return false
            end
        end

        if l1.xR2 <= l2.xR1
            x1 = l1.xR2
            x2 = l2.xR1
            R1 = abs(l1.R2)
            R2 = abs(l2.R1)
            y1 = l1.yR
            y2 = l2.yR
            d1 = l1.d
            d2 = l2.d
            x1min = l1.xmin
            x2min = l2.xmin
            x1max = l1.xmax
            x2max = l2.xmax
            y1min = l1.ymin
            y2min = l2.ymin
            y1max = l1.ymax
            y2max = l2.ymax
        else 
            x1 = l2.xR1
            x2 = l1.xR2
            R1 = abs(l2.R1)
            R2 = abs(l1.R2)
            y1 = l2.yR
            y2 = l1.yR
            d1 = l2.d
            d2 = l1.d
            x1min = l2.xmin
            x2min = l1.xmin
            x1max = l2.xmax
            x2max = l1.xmax
            y1min = l2.ymin
            y2min = l1.ymin
            y1max = l2.ymax
            y2max = l1.ymax
        end
        a = R1
        b = R2
        c = sqrt((x2 - x1)^2 + (y2 - y1)^2)
        d = (a^2 + c^2 - b^2)/(2c)
        f = sqrt(complex(a^2 - d^2))
        if !isreal(f)
            return false
        else
            f = Real(f)
        end

        f = Real(f)
        theta = abs(atan(f/d))
        phi = atan((y2 - y1)/(x2 - x1))
        gamma1 = phi + theta
        gamma2 = phi - theta
        xi1 = x1 + R1*cos(gamma1)
        xi2 = x1 + R1*cos(gamma2)
        yi1 = y1 + R1*sin(gamma1)
        yi2 = y1 + R1*sin(gamma2)

        yInRange = ((Real(yi1) < y1max) & (Real(yi1) > y1min) & (Real(yi1) < y2max) & (Real(yi1) > y2min)) | ((Real(yi2) < y1max) & (Real(yi2) > y1min) & (Real(yi2) < y2max) & (Real(yi2) > y2min)) 
        if yInRange
            return true
        else
            return false
        end
    elseif isinf(l1.R2)
        y1 = l1.yR
        d1 = l1.d
        x2 = l2.xR1
        y2 = l2.yR
        R2 = l2.R1
        d2 = l2.d
        xi = l1.xV2
        yi1 = sqrt(complex(R2^2 - (xi - x2)^2)) + y2
        yi2 = -sqrt(complex(R2^2 - (xi - x2)^2)) + y2
        if isreal(yi1)
            yInRange = ((Real(yi1) < y1max) & (Real(yi1) > y1min) & (Real(yi1) < y2max) & (Real(yi1) > y2min)) | ((Real(yi2) < y1max) & (Real(yi2) > y1min) & (Real(yi2) < y2max) & (Real(yi2) > y2min)) 
            if yInRange
                return true
            else
                return false
            end
        else
            return false
        end
    elseif isinf(l2.R1)
        x1 = l1.xR2
        y1 = l1.yR
        d1 = l1.d
        R1 = l1.R2
        y2 = l2.yR
        d2 = l2.d
        xi = l2.xV1
        yi1 = sqrt(complex(R1^2 - (xi - x1)^2)) + y1
        yi2 = -sqrt(complex(R1^2 - (xi - x1)^2)) + y1
        if isreal(yi)
            if ((Real(yi1) < y1 + d1/2) & (Real(yi1) > y1 - d1/2) & (Real(yi1) <= y2 + d2/2) & (Real(yi2) > y2 - d2/2)) | ((Real(yi2) < y1 + d1/2) & (Real(yi2) > y1 - d1/2) & (Real(yi2) <= y2 + d2/2) & (Real(yi2) > y2 - d2/2) )
                return true
            else
                return false
            end
        else
            return false
        end
    else
        return false
    end
end


"""
    overlap(l::SphericalLens...)
"""
function overlap(lenses::SphericalLens...)
    a = Array{Tuple{Real,Real},1}()
    for (i, lensi) in enumerate(lenses)
        for (j, lensj) in enumerate(lenses)
            if (i != j) & !((j, i) in a)
                if overlap(lensi, lensj)
                    push!(a, (i, j)) 
                end
            end
        end
    end
    return a
end


"""
    move(lens::SphericalLens, x::Real, y::Real)::SphericalLens
"""
function move(lens::SphericalLens, x::Real, y::Real)::SphericalLens
    lensCopy = SphericalLens(lens);
    dx = x - lensCopy.xV1
    dy = y - lensCopy.yR
    lensCopy.xV1 += dx;
    lensCopy.xV2 += dx;
    lensCopy.xR1 += dx;
    lensCopy.xR2 += dx;
    lensCopy.xe1 += dx;
    lensCopy.xe2 += dx;
    lensCopy.xmax += dx;
    lensCopy.xmin += dx;
    lensCopy.xmid += dx;
    lensCopy.yR += dy;
    lensCopy.ymin += dy;
    lensCopy.ymax += dy;
    return lensCopy
end


"""
    move(aperture::CircleAperture, x::Real, y::Real)::CircleAperture
"""
function move(aperture::CircleAperture, x::Real, y::Real)::CircleAperture
    aperutreCopy = CircleAperture(aperture);
    dx = x - apertureCopy.x
    dy = y - apertureCopy.y
    aperutreCopy.x += dx
    aperutreCopy.xmin += dx
    aperutreCopy.xmax += dx
    aperutreCopy.xmid += dx
    aperutreCopy.y += dy
    aperutreCopy.ymin += dy
    aperutreCopy.ymax += dy
    return aperutreCopy
end


"""
    move!(lens::SphericalLens, x::Real, y::Real)::SphericalLens
"""
function move!(lens::SphericalLens, x::Real, y::Real)::SphericalLens
    dx = x - lens.xV1
    dy = y - lens.yR
    lens.xV1 += dx;
    lens.xV2 += dx;
    lens.xR1 += dx;
    lens.xR2 += dx;
    lens.xe1 += dx;
    lens.xe2 += dx;
    lens.xmax += dx;
    lens.xmin += dx;
    lens.xmid += dx;
    lens.yR += dy;
    lens.ymin += dy;
    lens.ymax += dy;
    return lens
end


"""
    move!(aperture::CircleAperture, x::Real, y::Real)::CircleAperture
"""
function move!(aperture::CircleAperture, x::Real, y::Real)::CircleAperture
    dx = x - aperture.x
    dy = y - aperture.y
    aperutre.x += dx
    aperutre.xmin += dx
    aperutre.xmax += dx
    aperutre.xmid += dx
    aperutre.y += dy
    aperutre.ymin += dy
    aperutre.ymax += dy
    return aperutre
end


"""
    sphericalLensx{T<:Real}(l::SphericalLens, y::AbstractArray{T,1})
"""
function sphericalLensx{T<:Real}(l::SphericalLens, y::AbstractArray{T,1})
    @assert all(abs(y) .<= l.d/2) "any(abs(y) .> d/2) = $(l.d/2)"

    # x coords for surface 1 and 2
    if !isinf(l.R1)
        x1 = -sign(l.R1)*sqrt(l.R1^2 - y.^2) + l.xR1
    else
        x1 = zeros(y) + l.xV1
    end
    if !isinf(l.R2)
        x2 = -sign(l.R2)*sqrt(l.R2^2 - y.^2) + l.xR2 
    else
        x2 = zeros(y) + l.xV2
    end

    return (x1, x2)
end


"""
    sphericalLensVertexDiam(R1::Real, R2::Real, d::Real)::Real
Curvature center x coordinate for surface 1 with radii R1 and surface 2 with radii R2 and respective vertices for a lens with diameter d
"""
function sphericalLensVertexDiam(R1::Real, R2::Real, d::Real)::Tuple{Real,Real,Real,Real}
    if (R1 > 0) & (R2 < 0) # bi convex
        @assert R1^2 - d^2/4 > 0 "R1^2 - d^2/4 <= 0"
        @assert R2^2 - d^2/4 > 0 "R2^2 - d^2/4 <= 0"
        xR1 = sqrt(R1^2 - d^2/4)
        xR2 = -sqrt(R2^2 - d^2/4) 
    elseif (R1 < 0) & (R2 > 0) # bi concave
        xR1 = R1
        xR2 = R2
    elseif (R1 > 0) & (R2 > 0)
        @assert R1^2 - d^2/4 > 0 "R1^2 - d^2/4 <= 0"
        xR1 = sqrt(R1^2 - d^2/4)
        if abs(R2) < abs(R1)
            xR2 = xR1 - R1 + R2
        else
            xR2 = sqrt(R2^2 - d^2/4)
        end
    elseif (R1 < 0) & (R2 < 0)
        @assert R2^2 - d^2/4 > 0 "R2^2 - d^2/4 <= 0"
        xR1 = -sqrt(R1^2 - d^2/4)
        if abs(R2) < abs(R1)
            xR2 = -sqrt(R2^2 - d^2/4)
        else
            xR2 = xR1 - R1 + R2
        end
    end 
    
    if !isinf(xR1)
        xV1 = xR1 - R1
    else
        xV1 = 0
    end
    if !isinf(xR2) 
        xV2 = xR2 - R2
    else
        xV2 = 0
    end
    xR1 -= xV1
    xR2 -= xV1
    xV2 -= xV1
    return (xR1, xR2, 0, xV2)
end # function radiusCenter


"""
    sphericalLensVertexThick(R1::Real, R2::Real, t::Real)::Real
Curvature center x coordinate for surface 1 with radii R1 and surface 2 with radii R2 and respective vertices for a lens with thickness t
"""
function sphericalLensVertexThick(R1::Real, R2::Real, t::Real)::Tuple{Real,Real,Real,Real}
    xR1 = R1
    xR2 = t + R2
    xV1 = 0

    if !isinf(xR2) 
        xV2 = xR2 - R2
    else
        xV2 = t
    end

    return (xR1, xR2, xV1, xV2)
end


"""
    sphericalLensMinThick(R1::Real, R2::Real, d::Real)::Real
Minimum thickness of a lens with surface radii R1 and R2 and diameter d
"""
function sphericalLensMinThick(R1::Real, R2::Real, d::Real)
    (xR1, xR2, xV1, xV2) = sphericalLensVertexDiam(R1, R2, d)
    return abs(xV1 - xV2);
end # function sphericalLensMinThick


"""
    sphericalLensMaxDiam(R1::Real, R2::Real, t::Real)::Real
Maximum diameter of a lens with surface radii R1 and R2 and thickness t
"""
function sphericalLensMaxDiam(R1::Real, R2::Real, t::Real)
    (xR1, xR2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t)
    if (R1 < 0) & (R2 > 0) # bi-concave
        return maximum((abs(R1), R2))
    elseif isinf(R1) & isinf(R2) # planar lens
        return Inf
    elseif isinf(R1) & (R2 > 0) # plane-concave
        return R2
    elseif isinf(R2) & (R1 < 0) # concave-plane
        return R1
    elseif isinf(R1) & (R2 < 0) # plane-convex
        d = 2*sqrt(R2^2 - (xV1 - xR2)^2)
        return d
    elseif isinf(R2) & (R1 > 0) # convex-plane
        d = 2*sqrt(R1^2 - (xV2 - xR1)^2)
        return d
    elseif (R1 > 0) & (R2 < 0) # bi-convex
        xi = (R2^2 - R1^2 + xR1^2 - xR2^2)/(2*(xR1 - xR2))
        d = 2*sqrt(R1^2 - (xi - xR1)^2)
        return d
    elseif sign(R1)*sign(R2) > 0
        if R1 > 0
            if R2 > R1 
                xi = (R2^2 - R1^2 + xR1^2 - xR2^2)/(2*(xR1 - xR2))
                d = 2*sqrt(R1^2 - (xi - xR1)^2) 
                return d
            else
                return abs(R2)
            end
        elseif R1 < 0
            if R1 < R2 
                xi = (R2^2 - R1^2 + xR1^2 - xR2^2)/(2*(xR1 - xR2))
                d = 2*sqrt(R1^2 - (xi - xR1)^2) 
                return d
            else
                return abs(R1)
            end 
        end
    end
end # function sphericalLensMaxDiam


"""
    sphericalLensxEdge(R1::Real, R2::Real, d::Real)::Tuple{Real, Real} 
"""
function sphericalLensxEdge(R1::Real, R2::Real, d::Real)::Tuple{Real, Real} 
    (xR1, xR2, xV1, xV2) = sphericalLensVertexDiam(R1, R2, d);
    t = abs(xV1 - xV2);
    if !isinf(R1)
        xe1 = -sign(R1)*sqrt(R1^2 - d^2/4) + xR1;
    else
        xe1 = xV1;
    end
    if !isinf(R2)
        xe2 = -sign(R2)*sqrt(R2^2 - d^2/4) + xR2;
    else
        xe2 = xV2;
    end
    return (xe1, xe2)
end


"""
    sphericalLensxEdge(R1::Real, R2::Real, d::Real, t::Real)::Tuple{Real, Real} 
"""
function sphericalLensxEdge(R1::Real, R2::Real, d::Real, t::Real)::Tuple{Real, Real} 
    (xR1, xR2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t);
    t = abs(xV1 - xV2);
    if !isinf(R1)
        xe1 = -sign(R1)*sqrt(R1^2 - d^2/4) + xR1;
    else
        xe1 = xV1;
    end
    if !isinf(R2)
        xe2 = -sign(R2)*sqrt(R2^2 - d^2/4) + xR2;
    else
        xe2 = xV2;
    end
    return (xe1, xe2)
end


"""
    sphericalLensfo2fi(fo::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensfo2fi(fo::Real, nm1::Real, nm2::Real)::Real
    fi = nm2*fo/nm1
    return fi
end


"""
    sphericalLensfo2fi(fi::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensfi2fo(fi::Real, nm1::Real, nm2::Real)::Real
    fo = nm1*fi/nm2
    return fo
end


"""
    sphericalLensfo(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real

Object focal length fo from primary principal plane H1 from surface 1 radius of curvature R1, surface 2 radius of curvature R2, lens thickness t, lens index of refraction nl, medium 1 index of refraction nm1 and medium 2 index of refraction.
"""
function sphericalLensfo(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    n12 = nl - nm1;
    n23 = nm2 - nl;
    x = n12/R1 + n23/R2 - t*n12*n23/(nl*R1*R2)
    fo = nm1/x
    return fo
end


"""
    sphericalLensfo(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensfo(l::SphericalLens, nm1::Real, nm2::Real)::Real
    fo = sphericalLensfo(l.R1, l.R2, l.t, l.nl, nm1, nm2)
    return fo
end


"""
    sphericalLensfi(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real

Image focal length fo from primary principal plane H1 from surface 1 radius of curvature R1, surface 2 radius of curvature R2, lens thickness t, lens index of refraction nl, medium 1 index of refraction nm1 and medium 2 index of refraction.
"""
function sphericalLensfi(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    fo = sphericalLensfo(R1, R2, t, nl, nm1, nm2)
    fi = nm1*fo/nm2
    return fi
end 


"""
    sphericalLensfi(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensfi(l::SphericalLens, nm1::Real, nm2::Real)::Real
    fi = sphericalLensfi(l.R1, l.R2, l.t, l.nl, nm1, nm2)
    return fi
end


"""
    sphericalLensRfo(fo::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real

Radius of curvature for surface 1 and 2 where R1 = R and R2 = -R from primary (object) focal length fo, thickness t, lens index of refraction nl, medium 1 index of refraction nm1 and medium 2 index of refraction.
"""
function sphericalLensRfo(fo::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    n12 = nl - nm1;
    n23 = nm2 - nl;
    a = t/nl*n12*n23;
    b = n12 - n23;
    c = -nm1/fo;
    if a == 0
        x = -c/b
    else
        @assert b^2 > 4*a*c "(2nl - nm1 - nm2)^2 <= 4*t/nl*(nl - nm1)*(nm2 - nl)*nm1/f"
        x = (sqrt(b^2/4 - c*a) - b/2)/a
    end
    R = 1/x
    return R
end


function sphericalLensRfi(fi::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    n12 = nl - nm1;
    n23 = nm2 - nl;
    a = t/nl*n12*n23;
    b = n12 - n23;
    c = -nm2/fi;
    if a == 0
        x = -c/b
    else
        @assert b^2 > 4*a*c "(2nl - nm1 - nm2)^2 <= 4*t/nl*(nl - nm1)*(nm2 - nl)*nm1/f"
        x = (sqrt(b^2/4 - c*a) - b/2)/a
    end
    R = 1/x
    return R
end


function sphericalLensR1(fo::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)
    n12 = nl - nm1;
    n23 = nm2 - nl;
    a = (nm1/fo - n23/R2)
    b = n12*(1 - t*n23/(R2*nl))
    R1 = b/a
    return R1
end


function sphericalLensR2(fo::Real, R1::Real, t::Real, nl::Real, nm1::Real, nm2::Real)
    n12 = nl - nm1;
    n23 = nm2 - nl;
    a = (nm1/fo - n12/R1)
    b = n23*(1 - t*n12/(R1*nl))
    R1 = b/a
    return R1
end


"""
    sphericalLensxV1xH1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real

Primary principal plane distace from the primary (surface 1) vertex
"""
function sphericalLensxV1xH1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    fo = sphericalLensfo(R1, R2, t, nl, nm1, nm2)
    n23 = nm2 - nl
    f2p = R2*nl/n23
    xV1xH1 = fo*t/f2p
    return xV1xH1
end


"""
    sphericalLensxV1xH1(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxV1xH1(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xV1xH1 = sphericalLensxV1xH1(l.R1, l.R2, l.t, l.nl, nm1, nm2)
    return xV1xH1
end


"""
    sphericalLensxV2xH2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real

Secondary principal plane distace from the secondary (surface 2) vertex
"""
function sphericalLensxV2xH2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    fi = sphericalLensfi(R1, R2, t, nl, nm1, nm2)
    n12 = nl - nm1
    f1s = R1*nl/n12
    xV2xH2 = -fi*t/f1s
    return xV2xH2
end


"""
    sphericalLensxV2xH2(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxV2xH2(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xV2xH2 = sphericalLensxV2xH2(l.R1, l.R2, l.t, l.nl, nm1, nm2)
    return xV2xH2
end


"""
    sphericalLensxH1xN1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH1xN1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    fi = sphericalLensfi(R1, R2, t, nl, nm1, nm2)
    xH1xN1 = fi*(nm2 - nm1)/nm2
    return xH1xN1
end 


"""
    sphericalLensxH1xN1(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH1xN1(l::SphericalLens, nm1::Real, nm2::Real)::Real
    fi = sphericalLensfi(l.R1, l.R2, l.t, l.nl, nm1, nm2)
    xH1xN1 = fi*(nm2 - nm1)/nm2
    return xH1xN1
end 


"""
    sphericalLensxH2xN2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH2xN2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real 
    xH2xN2 = sphericalLensxH1xN1(R1, R2, t, nl, nm1, nm2)
    return xH2xN2
end 


"""
    sphericalLensxH2xN2(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH2xN2(l::SphericalLens, nm1::Real, nm2::Real)::Real 
    xH2xN2 = sphericalLensxH1xN1(l, nm1, nm2)
    return xH2xN2
end 


"""
    sphericalLensxH1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    (R1, R2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t)
    xV1xH1 = sphericalLensxV1xH1(R1, R2, t, nl, nm1, nm2)
    xH1 = xV1xH1 + xV1
    return xH1
end


"""
    sphericalLensxH1(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH1(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xV1xH1 = sphericalLensxV1xH1(l, nm1, nm2)
    xH1 = xV1xH1 + l.xV1
    return xH1
end


"""
    sphericalLensxH2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    (R1, R2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t)
    xV2xH2 = sphericalLensxV2xH2(R1, R2, t, nl, nm1, nm2)
    xH2 = xV2xH2 + xV2
    return xH2
end


"""
    sphericalLensxH2(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxH2(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xV2xH2 = sphericalLensxV2xH2(l, nm1, nm2)
    xH2 = xV2xH2 + l.xV2
    return xH2
end


"""
    sphericalLensxN1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxN1(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    xH1xN1 = sphericalLensxH1xN1(R1, R2, t, nl, nm1, nm2)
    xH1 = sphericalLensxH1(R1, R2, t, nl, nm1, nm2)
    xN1 = xH1 + xH1xN1
    return xN1
end 


"""
    sphericalLensxN1(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxN1(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xH1xN1 = sphericalLensxH1xN1(l, nm1, nm2)
    xH1 = sphericalLensxH1(l, nm1, nm2)
    xN1 = xH1 + xH1xN1
    return xN1
end 


"""
    sphericalLensxN2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxN2(R1::Real, R2::Real, t::Real, nl::Real, nm1::Real, nm2::Real)::Real
    xH2xN2 = sphericalLensxH1xN1(R1, R2, t, nl, nm1, nm2)
    xH2 = sphericalLensxH2(R1, R2, t, nl, nm1, nm2)
    xN2 = xH2 + xH2xN2
    return xN2
end 


"""
    sphericalLensxN2(l::SphericalLens, nm1::Real, nm2::Real)::Real
"""
function sphericalLensxN2(l::SphericalLens, nm1::Real, nm2::Real)::Real
    xH2xN2 = sphericalLensxH2xN2(l, nm1, nm2)
    xH2 = sphericalLensxH2(l, nm1, nm2)
    xN2 = xH2 + xH2xN2
    return xN2
end 


"""
    sphericalLenssi(so::Real, fo::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLenssi(so::Real, fo::Real, nl::Real, nm1::Real, nm2::Real)::Real
    x = nm1*(1/fo - 1/so);
    si = nm2/x;
    return si
end


"""
    sphericalLensso(si::Real, fi::Real, nl::Real, nm1::Real, nm2::Real)::Real
"""
function sphericalLensso(si::Real, fi::Real, nl::Real, nm1::Real, nm2::Real)::Real
    x = nm2*(1/fi - 1/si);
    so = nm1/x;
    return so
end


"""
    sphericalLensx(l::SphericalLens, y::Real)

x coordinates of the first and second lens surface corresponding to the y coordinate given
"""
function sphericalLensx(l::SphericalLens, y::Real)
    @assert abs(y) <= l.d/2 "abs(y) > d/2 = $(l.d/2)"

    # x coords for surface 1 and 2
    if !isinf(l.R1)
        x1 = -sign(l.R1)*sqrt(l.R1^2 - y.^2) + l.xR1
    else
        x1 = zero(y) + l.xV1
    end
    if !isinf(l.R2)
        x2 = -sign(l.R2)*sqrt(l.R2^2 - y.^2) + l.xR2 
    else
        x2 = zero(y) + l.xV2
    end

    return (x1, x2)
end


"""
    lineLineIntersect(x0, y0, angle0, x1, y1, angle1)

Intersection (x, y) coordinates of two lines passing through (x0, y0) and (x1, y1) with angle a0 and a1 respectively.
"""
function lineLineIntersect(x0, y0, angle0, x1, y1, angle1)
    @assert !isapprox(angle0, angle1) "angle0 ≈ angle1"
    s0 = tan(angle0)
    s1 = tan(angle1)

    if isapprox(angle0, pi/2) | isapprox(angle0, -pi/2)
        return lineLineIntersect(x1, y1, angle1, x0)
    elseif isapprox(angle1, pi/2) | isapprox(angle1, -pi/2)
        return lineLineIntersect(x0, y0, angle0, x1) 
    else
        a0 = -s0
        a1 = -s1
        b0 = 1
        b1 = 1
        c0 = y0 - a0*x0
        c1 = y1 - a1*x1
        A = [a0 b0; a1 b1]
        b = [c0; c1]
        (xi, yi) = A\b
        return (xi, yi)
    end
end


"""
    lineLineIntersect(x0, y0, angle0, x1)

Intersection (x, y) coordinates of a line passing through (x0, y0) with angle a0 and a vertical line passing through x1 respectively.
"""
function lineLineIntersect(x0, y0, angle0, x1)
    s0 = tan(angle0)
    xi = x1
    yi = s0*(xi - x0) + y0
    return (xi, yi)
end


#function lineCircleIntersect(x0, y0, angle0, xR, yR, R)
#    s = tan(angle0)
#    if isapprox(angle0, pi/2)
#        x = x0
#        y1 = sqrt((R^2 - (x - xR)^2)) + yR
#        y2 = -sqrt((R^2 - (x - xR)^2)) + yR
#        return ((x, y1), (x, y2))
#    else
#        a = 1 + s^2
#        b = -2*(s^2*x0 + s*y0 + s*yR + xR)
#        c = -R^2 + s^2*x0^2 + 2*s*x0*y0 + 2*s*x0*yR + xR^2 + y0^2 + 2*y0*yR + yR^2 
#        f = sqrt((-4*a*c + b^2))
#        x1 = -(b + f)/(2*a)
#        x2 = (-b + f)/(2*a)
#        y1 = s*(x1 - x0) + y0
#        y2 = s*(x2 - x0) + y0
#        return ((x1, y1), (x2, y2))
#    end
#end

"""
    lineCircleIntersect(x0, y0, angle0, xR, yR, R)

Intersection coordinates of a line and circle. ((x1, y1), (x2, y2)) is returned. The line has angle angle0 (in radians) and passes through point at (x0, y0). Circle with radius R has center at (xR, yR). Two, one or no intersection points can exist between the line and circle. A missing point of intersection is represented as (nothing, nothing).

Source: https://en.wikipedia.org/wiki/Intersection_(Euclidean_geometry)#A_line_and_a_circle

# Examples
```jldoctest
julia> lineCircleIntersect(0, 0, 0, 0, 0, 1)
((-1.0,0.0),(1.0,0.0))

```
"""
function lineCircleIntersect(x0, y0, angle0, xR, yR, R)
    if abs(angle0) > 2pi
        # remainder of angle0/(2pi) 
        angle0 = angle0%(2pi)
    end

    xa = x0 - xR
    ya = y0 - yR
    xb = xa + 1
    yb = tan(angle0)*1 + ya
    dx = 1
    dy = yb - ya
    dr = sqrt(dx^2 + dy^2)
    D = xa*yb - xb*ya
    
    a = sqrt(complex(R^2*dr^2 - D^2))
    if !isreal(a)
        x1 = nothing
        y1 = nothing
        x2 = nothing
        y2 = nothing
        return ((x1, y1), (x2, y2))
    else
        a = Real(a)
    end

    sgn = sign(dy)
    if sgn == zero(sgn)
        sgn += 1
    end
    b = sign(dy)*dx*a
    c = abs(dy)*a

    if (angle0 >= 0)
        x1 = (D*dy - a)/dr^2 + xR
        x2 = (D*dy + a)/dr^2 + xR
        y1 = (-D*dx - c)/dr^2 + yR
        y2 = (-D*dx + c)/dr^2 + yR
    else
        x1 = (D*dy - a)/dr^2 + xR
        x2 = (D*dy + a)/dr^2 + xR
        y1 = (-D*dx + c)/dr^2 + yR
        y2 = (-D*dx - c)/dr^2 + yR
    end

    # Ray in positive x direction
    if (angle0 <= pi/2) & (angle0 >= -pi/2)
        return ((x1, y1), (x2, y2))
    # Ray in negative x direction
    else  
        return ((x2, y2), (x1, y1))
    end
end


""" 
    circleNormAngle(y, xR, yR, R) 

Normal angles (in radians) clockwise from x axis at y for circle with centre (xR, yR) and radius R. (normalAngle1, normalAngle2) is returned. normalAngle1 is the left side normal angle and normalAngle2 is the right side normal angle
"""
function circleNormAngle(y, xR, yR, R)
    normAngle1 = -asin((y - yR)/abs(R))
    normAngle2 = -normAngle1
    return (normAngle1, normAngle2)
end


"""
    lineCircleNormAngle(x0, y0, angle0, xR, yR, R)

Angle between line and circle normal at the two points of intersection between line and circle. (angle1, angle2, normalAngle1, normalAngle2) is returned. The line has angle angle0 (in radians) and passes through point at (x0, y0). Circle with radius R has center at (xR, yR). Two, one or no intersection points can exist between the line and circle. A missing point of intersection is represented as (nothing, nothing).

# Examples
```jldoctest
julia> lineCircleNormAngle()

```
"""
function lineCircleNormAngle(x0, y0, angle0, xR, yR, R)
    ((x1, y1), (x2, y2)) = lineCircleIntersect(x0, y0, angle0, xR, yR, R)
    if !isa(y1, Void)
        (normAngle11, normAngle12) = circleNormAngle(y1, xR, yR, R)
        normAngle1 = normAngle11
        angle1 = angle0 - normAngle1
    else
        normAngle1 = nothing
        angle1 = nothing
    end

    if !isa(y2, Void)
        (normAngle21, normAngle22) = circleNormAngle(y2, xR, yR, R)
        normAngle2 = normAngle22
        angle2 = angle0 - normAngle2
    else
        normAngle2 = nothing
        angle2 = nothing
    end

    return (angle1, angle2, normAngle1, normAngle2)
end


"""
    rayLensIntersect(x0::Real, y0::Real, angle0::Real, l::SphericalLens)

Intersection coordinates for a ray and lens. The ray originates at (x0, y0) with angle angle0 (in radians). The intersection coordinates are returned as ((x1, y1), (x2, y2)). Two, one or no intersection points can exist between the ray and lens. A missing intersection is represented and returned as (nothing, nothing).
"""
function rayLensIntersect(x0::Real, y0::Real, angle0::Real, l::SphericalLens)
    yR = l.yR
    xR1 = l.xR1
    R1 = l.R1
    xR2 = l.xR2
    R2 = l.R2
    d = l.d
    
    # left surface is curved
    if !isinf(R1)
        points = lineCircleIntersect(x0, y0, angle0, xR1, yR, R1)
        if R1 > 0
            # if radius positive use left intersection
            (x1, y1) = points[1]
        else
            # else use right intersection
            (x1, y1) = points[2]
        end
    # left surface is plane
    else
        x1 = l.xV1
        (x1, y1) = lineLineIntersect(x0, y0, angle0, x1)
    end

    # right surface is curved
    if !isinf(R2)
        points = lineCircleIntersect(x0, y0, angle0, xR2, yR, R2)
        if R2 > 0
            # if radius positive use left intersection 
            (x2, y2) = points[1]
        else
            # else use right intersection
            (x2, y2) = points[2]
        end
    # right surface is plane
    else
        x2 = l.xV2
        (x2, y2) = lineLineIntersect(x0, y0, angle0, x2)
    end

    # check intersection is within lens 
    if (x1, y1) != (nothing, nothing) 
        if (y1 > l.ymax) | (y1 < l.ymin) | (x1 < l.xmin) | (x1 > l.xmax)
            x1 = nothing
            y1 = nothing
        end
    end
    if (x2, y2) != (nothing, nothing) 
        if (y2 > l.ymax) | (y2 < l.ymin) | (x2 < l.xmin) | (x2 > l.xmax)
            x2 = nothing
            y2 = nothing
        end
    end

    return ((x1, y1), (x2, y2)) 
end


"""
    rayLensAngle(x0::Real, y0::Real, angle0::Real, l::SphericalLens)

Angle between ray (line) and lens surface normal at each surface intersection. Order is left surface then right surface. Returns (angle 1, angle 2, normal angle 1, normal angle 2). If an intersection with one of the surface does not exist then nothing is returned for the corresponding angle and normal angle.
"""
function rayLensAngle(x0::Real, y0::Real, angle0::Real, l::SphericalLens)
    ((x1, y1), (x2, y2)) = rayLensIntersect(x0, y0, angle0, l)
    if !isa(x1, Void) 
        if !isinf(l.R1)
            (angle11, angle12, normAngle11, normAngle12) = lineCircleNormAngle(x0, y0, angle0, l.xR1, l.yR, l.R1) 
            if l.R1 > 0
                angle1 = angle11
                normAngle1 = normAngle11
            else
                angle1 = angle12
                normAngle1 = normAngle12
            end
        else
            angle1 = angle0
            normAngle1 = 0
        end
    else
        angle1 = nothing
        normAngle1 = nothing
    end

    if !isa(x2, Void)
        if !isinf(l.R2)
            (angle21, angle22, normAngle21, normAngle22) = lineCircleNormAngle(x0, y0, angle0, l.xR2, l.yR, l.R2) 
            if l.R2 > 0
                angle2 = angle21
                normAngle2 = normAngle21
            else
                angle2 = angle22
                normAngle2 = normAngle22
            end
        else
            angle2 = angle0
            normAngle2 = 0
        end
    else
        angle2 = nothing
        normAngle2 = nothing 
    end

    return (angle1, angle2, normAngle1, normAngle2)
end


"""
    rayApertureIntersect(x0::Real, y0::Real, angle0::Real, a::CircleAperture)
"""
function rayApertureIntersect(x0::Real, y0::Real, angle0::Real, a::CircleAperture)
    (x1, y1) = lineLineIntersect(x0, y0, angle0, a.x)
    if (((y1 >= (a.y + a.d1/2)) & (y1 <= (a.y + a.d2/2))) | ((y1 <= (a.y - a.d1/2)) & (y1 >= (a.y - a.d2/2))))
        return (x1, y1)
    else
        (x1, y1) = (nothing, nothing)
        return (x1, y1)
    end
end


"""
    OpticsSystem
"""
typealias OpticsSystem Array{Optic,1}


"""
    OpticsSystem(os::Array{Optic,1})::OpticsSystem
"""
function OpticsSystem(os::Array{Optic,1})::OpticsSystem
    return sort(os)
end


"""
    sort(os::OpticsSystem)::OpticsSystem

Optics system sorted by optic (min x, min y) coordinate tuple
"""
function sort(os::OpticsSystem)::OpticsSystem
    sort(os, by = l -> (l.xmin, l.ymin))
end


"""
    sort!(os::OpticsSystem)

Optics system sorted by optic (min x, min y) coordinate tuple
"""
function sort!(os::OpticsSystem)
    sort!(os, by = l -> (l.xmin, l.ymin))
end


"""
    flipVert(l::SphericalLens)
"""
function flipVert(l::SphericalLens)
    lr = SphericalLens(-l.R2, -l.R1, l.t, l.d, l.nl)
    lr = move(lr, l.xV1, 0)
    return lr
end


"""
    flipVert(a::CircleAperture)
"""
function flipVert(a::CircleAperture)
    return CircleAperture(a)
end


"""
    flipx(os::OpticsSystem, x::Real=0)::OpticsSystem
"""
function flipx(os::OpticsSystem, x::Real=0)::OpticsSystem
    os = sort(os)
    osn = OpticsSystem(length(os))
    xMidSysMax = os[1].xmid
    xMidSysMin = os[end].xmid
    xMidPosList = Array{Real,1}(length(os))
    for (i, o) in enumerate(os)
        xMidPosList[i] = o.xmid
        if o.xmid > xMidSysMax
            xMidSysMax = o.xmid
        end
        if o.xmid < xMidSysMin
            xMidSysMin = o.xmid
        end
    end

    xMidSys = (xMidSysMin + xMidSysMax)/2

    for (i, (o, xMidPos)) in enumerate(zip(os, xMidPosList))
        on = flipVert(o)
        if isa(on, SphericalLens)
            osn[i] = move(on, xMidSys - (xMidPos - xMidSys) - on.t/2, 0)
        elseif isa(on, CircleAperture)
            osn[i] = move(on, xMidSys - (xMidPos - xMidSys), 0)
        end
    end

    osn = sort(osn)

    return osn
end


"""
    OpticsSystemMatrix
"""
typealias OpticsSystemMatrix Array{Array{Number,2},1}


"""
    OpticsSystemMatrix(l::SphericalLens, nm1, nm2)
"""
function OpticsSystemMatrix(l::SphericalLens, nm1, nm2)
    A1 = refractMatrix(l.R1, nm1, l.nl)
    A2 = transferMatrix(l.t)
    A3 = refractMatrix(l.R2, l.nl, nm2)
    A = OpticsSystemMatrix([A3, A2, A1]) 
end


"""
    OpticsSystemMatrix{N<:Number}(os::OpticsSystem, n::Array{N,1})
"""
function OpticsSystemMatrix{N<:Number}(os::OpticsSystem, n::Array{N,1})
    SM = OpticsSystemMatrix()
    transferOpticsInd = find([isa(o, Lens) for o=os])
    lastInd = transferOpticsInd[end]
    for j=1:length(transferOpticsInd)
        i = transferOpticsInd[j]
        if i != lastInd
            iNext = transferOpticsInd[j+1]
        end
        if isa(os[i], Lens)
            if (i != lastInd) # not the last optic in system
                M = OpticsSystemMatrix(os[i], n[i], n[i + 1]) # optic matrix
                T = transferMatrix(os[iNext].xV1 - os[i].xV2) # transfer matrix to next optic
                prepend!(SM, M)
                prepend!(SM, [T])
            elseif (i == lastInd) # last optic in system
                M = OpticsSystemMatrix(os[i], n[i], n[i + 1]) # optic matrix
                prepend!(SM, M)
            end
        end
    end
    return SM
end


"""
    rayMatrix(y::Real, angle::Real)::Array{Number} 

Creates a ray vector r=[angle, y] for a ray with an angle measured clockwise from the optical axis (in radians) emanating a height y above the optical axis

# Example
```jldoctest
julia> rayMatrix(1, pi/8)
2-element Array{Real,1}:
 0.392699
 1.0

```
"""
function rayMatrix(y::Real, angle::Real)::Array{Real}
    #A = [n*angle; y]
    A = [angle; y]
    return A
end


"""
    refractMatrix(R::Real, n1::Real, n2::Real)::Array{Number} 

Refraction matrix A for a spherical surface with radius R and indices of refraction n1 and n2 on each side respectively 

A = [n1/n2 -(n2 - n1)/(R*n2); 0 1]

# Example
```jldoctest
julia> Optics.refractMatrix(1, 1, 1.5)
2×2 Array{Number,2}:
 0.666667  -0.333333
 0.0        1.0

```
"""
function refractMatrix(R::Real, n1::Real, n2::Real)::Array{Number}
    D = (n2 - n1)/R
    #A = [1 -D; 0 1] # Hecht definition
    A = [n1/n2 -D/n2; 0 1]
    return A
end


"""
    transferMatrix(t::Real)::Array{Number}

Ray transfer matrix for a distance t parallel to the optical axis

A = [1 0; t 1]

# Example
```jldoctest
julia> Optics.transferMatrix(0.1)
2×2 Array{Number,2}:
 1.0  0.0
 0.1  1.0

```
"""
function transferMatrix(t::Real)::Array{Number}
    #A = [1 0; t/n 1] # Hecht definition
    A = [1 0; t 1]
    return A
end


"""
    opticsSystemfo{T<:Number}(os::OpticsSystem, nm::Array{T,1})

Object (primary) focal length for optics system os with indices of refraction between optics given by nm. Focal length is measured from the primary principal point H1 of the optics system. length(nm) = length(os) + 1.

# Example
```jldoctest
julia> l = SphericalLens(1, -1, 1, 1.5)
Optics.SphericalLens(1,-1,0.2679491924311228,1,1.5,0,0.2679491924311228,1.0,-0.7320508075688772,0.1339745962155614,0.1339745962155614,0.0,0.2679491924311228,0.1339745962155614,0,-0.5,0.5)

julia> os = OpticsSystem([l])
1-element Array{Optics.Optic,1}:
 Optics.SphericalLens(1,-1,0.267949,1,1.5,0,0.267949,1.0,-0.732051,0.133975,0.133975,0.0,0.267949,0.133975,0,-0.5,0.5)

julia> Optics.opticsSystemfo(os, [1, 1])
-1.0467457811220566

```
"""
function opticsSystemfo{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    fo = nm[1]/A[1,2]
    return fo
end


"""
    opticsSystemfi{T<:Number}(os::OpticsSystem, nm::Array{T,1})

Image (secondary) focal length for optics system os with indices of refraction between optics given by nm. Focal length is measured from the secondary principal point H2 of the optics system. length(nm) = length(os) + 1.

# Example
```jldoctest
julia> l = SphericalLens(1, -1, 1, 1.5)
Optics.SphericalLens(1,-1,0.2679491924311228,1,1.5,0,0.2679491924311228,1.0,-0.7320508075688772,0.1339745962155614,0.1339745962155614,0.0,0.2679491924311228,0.1339745962155614,0,-0.5,0.5)

julia> os = OpticsSystem([l])
1-element Array{Optics.Optic,1}:
 Optics.SphericalLens(1,-1,0.267949,1,1.5,0,0.267949,1.0,-0.732051,0.133975,0.133975,0.0,0.267949,0.133975,0,-0.5,0.5)

julia> Optics.opticsSystemfi(os, [1, 1])
1.0467457811220566

```
"""
function opticsSystemfi{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    fi = -nm[end]/A[1,2]
    return fi
end


"""
    opticsSystemxV1xH1{T<:Number}(os::OpticsSystem, nm::Array{T,1})

Distance from the first (left most) vertex V1 to the primary principal point H1 of the optics system.

# Example
```jldoctest
julia> l = SphericalLens(1, -1, 1, 1.5)
Optics.SphericalLens(1,-1,0.2679491924311228,1,1.5,0,0.2679491924311228,1.0,-0.7320508075688772,0.1339745962155614,0.1339745962155614,0.0,0.2679491924311228,0.1339745962155614,0,-0.5,0.5)

julia> os = OpticsSystem([l])
1-element Array{Optics.Optic,1}:
 Optics.SphericalLens(1,-1,0.267949,1,1.5,0,0.267949,1.0,-0.732051,0.133975,0.133975,0.0,0.267949,0.133975,0,-0.5,0.5)

julia> Optics.opticsSystemxV1xH1(os, [1, 1])
0.09349156224411338

```
"""
function opticsSystemxV1xH1{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    xV1xH1 = -nm[1]*(1 - A[1,1])/A[1,2]
    return xV1xH1
end


"""
    opticsSystemxV2xH2{T<:Number}(os::OpticsSystem, nm::Array{T,1})
"""
function opticsSystemxV2xH2{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    xV2xH2 = -nm[end]*(A[2,2] - 1)/A[1,2]
    return xV2xH2
end


"""
    opticsSystemffl{T<:Number}(os::OpticsSystem, nm::Array{T,1})
"""
function opticsSystemffl{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    fo = opticsSystemfo(os, nm) 
    ffl = A[1,1]*fo
    return ffl
end


"""
    opticsSystembfl{T<:Number}(os::OpticsSystem, nm::Array{T,1})
"""
function opticsSystembfl{T<:Number}(os::OpticsSystem, nm::Array{T,1})
    A = *(OpticsSystemMatrix(os, nm)...)
    fi = opticsSystemfi(os, nm) 
    bfl = A[2,2]*fi
    return bfl
end


"""
    matrixRayTrace{T<:Number}(os::OpticsSystem, nm::Array{T,1}, x0, y0, angle0, xEnd)
"""
function matrixRayTrace{T<:Number}(os::OpticsSystem, nm::Array{T,1}, x0, y0, angle0, xEnd)
    xs = Array{Number,1}()
    ys = Array{Number,1}()
    as = Array{Number,1}()
    rays = Array{Array{Number},1}()
    push!(xs, x0)
    push!(ys, y0)
    push!(as, angle0)

    i = 1
    ray = rayMatrix(y0, angle0)
    push!(rays, ray)
    for o in os
        if isa(o, SphericalLens)
            xV1 = o.xV1
            xV2 = o.xV2
            #xH1 = sphericalLensxH1(o, nm[i], nm[i + 1])
            #xH2 = sphericalLensxH2(o, nm[i], nm[i + 1])

            (x1, y1) = lineLineIntersect(xs[end], ys[end], as[end], xV1)
            if (y1 > (o.yR + o.d/2)) | (y1 < (o.yR - o.d/2))
                (x1, y1) = (nothing, nothing)
            end

            if (x1, y1) != (nothing, nothing) 
                # transfer to primary vertex plane
                ray = transferMatrix(xV1 - xs[end])*ray
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 

                # refraction at the primary vertex plane (first lens surface)
                ray = refractMatrix(o.R1, nm[i], o.nl)*ray
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 

                # transfer from first vertex plane to second vertex plane
                ray = transferMatrix(xV2 - xV1)*ray
                (x1, y1) = lineLineIntersect(xs[end], ys[end], as[end], xV2)
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 

                # refraction at the vertex principal plane (second lens surface)
                ray = refractMatrix(o.R2, o.nl, nm[i+1])*ray
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 
            else
                # transfer to primary vertex plane
                (x1, y1) = lineLineIntersect(xs[end], ys[end], as[end], xV1)
                ray = transferMatrix(xV1 - xs[end])*ray
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 

                # transfer to second vertex plane
                ray = transferMatrix(xV2 - xV1)*ray
                (x1, y1) = lineLineIntersect(xs[end], ys[end], as[end], xV2)
                push!(xs, x1)
                push!(ys, ray[2])
                push!(as, ray[1])
                push!(rays, ray) 
            end
        elseif isa(o, Aperture)
            (x1, y1) = rayApertureIntersect(xs[end], ys[end], as[end], o)
            if (x1, y1) != (nothing, nothing)
                push!(xs, x1)
                push!(ys, y1)
                push!(as, as[end])
                return (xs, ys, as)
            end
        end 
    end

    # transfer ray to the end
    ray = transferMatrix(xEnd - xs[end])*ray
    push!(xs, xEnd)
    push!(ys, ray[2])
    push!(as, ray[1])
    push!(rays, ray) 

    return (xs, ys, as)
end 


"""
    snell(angle1, n1, n2)
"""
function snell(angle1, n1, n2)
    return asin(n1*sin(angle1)/n2)
end


"""
    rayLensRefract(x0, y0, angle0, l::SphericalLens, n1, n2)
"""
function rayLensRefract(x0, y0, angle0, l::SphericalLens, n1, n2)
    xs = Array{Number,1}([x0])
    ys = Array{Number,1}([y0])
    as = Array{Number,1}([angle0])

    ((x1a, y1a), (x2a, y2a)) = rayLensIntersect(xs[end], ys[end], as[end], l)
    if !isa(x1a, Void) 
        (angle1a, angle2a, normAngle1a, normAngle2a) = rayLensAngle(x1a, y1a, as[end], l)
        a1 = snell(angle1a, n1, l.nl) + normAngle1a


        ((x1b, y1b), (x2b, y2b)) = rayLensIntersect(x1a, y1a, a1, l)
        if !isa(x2b, Void)
            push!(xs, x1a)
            push!(ys, y1a)
            push!(as, a1)
            (angle1b, angle2b, normAngle1b, normAngle2b) = rayLensAngle(xs[end], ys[end], as[end], l)
            a2 = snell(angle2b, l.nl, n2) + normAngle2b
            push!(xs, x2b)
            push!(ys, y2b)
            push!(as, a2)
        else
            push!(xs, l.xe2)
            push!(ys, tan(as[end])*(xs[end] - xs[end-1]) + ys[end])
            push!(as, as[end])
            return (xs, ys, as) 
        end

        return (xs, ys, as) 
    else
        push!(xs, l.xe2)
        push!(ys, tan(as[end])*(xs[end] - xs[end-1]) + ys[end])
        push!(as, as[end])
        return (xs, ys, as) 
    end 
end


"""
    rayTrace{T<:Number}(os::OpticsSystem, nm::Array{T,1}, x0::Real, y0::Real, angle0::Real, xEnd::Real)
"""
function rayTrace{T<:Number}(os::OpticsSystem, nm::Array{T,1}, x0::Real, y0::Real, angle0::Real, xEnd::Real)
    xs = Array{Number,1}()
    ys = Array{Number,1}()
    as = Array{Number,1}()
    push!(xs, x0)
    push!(ys, y0)
    push!(as, angle0)
    for (i, o) in enumerate(os) 
        if isa(o, Lens)
            ((x1, y1), (x2, y2)) = rayLensIntersect(xs[end], ys[end], as[end], o)
            if (x1, y1) != (nothing, nothing)
                if (ys[end] < o.ymax) & (ys[end] > o.ymin)
                    (x, y, a) = rayLensRefract(xs[end], ys[end], as[end], o, nm[i], nm[i + 1])

                    # store results
                    push!(xs, x...)
                    push!(ys, y...)
                    push!(as, a...)
                end
            end
        elseif isa(o, Aperture)
            (x1, y1) = rayApertureIntersect(xs[end], ys[end], as[end], o)
            if (x1, y1) != (nothing, nothing)
                push!(xs, x1)
                push!(ys, y1)
                push!(as, as[end])
                return (xs, ys, as)
            end
        end
    end

    if xs[end] < xEnd
        x = xEnd
        y = tan(as[end])*(x - xs[end]) + ys[end]
        a = as[end]
        push!(xs, x)
        push!(ys, y)
        push!(as, a)
    end

    return (xs, ys, as)
end


"""
    draw(R1::Real, R2::Real, d::Real, t::Real=0, x::Real=0, y::Real=0, ny::Integer=100)

Spherical lens surface coordinates for plotting
"""
function draw(R1::Real, R2::Real, d::Real, t::Real=0, x::Real=0, y::Real=0, ny::Integer=100)
    @assert d <= sphericalLensMaxDiam(R1, R2, t) "d > dmax"
    (xR1, xR2, xV1, xV2) = sphericalLensVertexThick(R1, R2, t)

    # x and y coords for surface 1 and 2
    y1 = linspace(-d/2, d/2, ny)
    if !isinf(R1)
        x1 = -sign(R1)*sqrt(R1^2 - y1.^2) + xR1
    else
        x1 = zeros(y1) + xV1
    end
    if !isinf(R2)
        x2 = -sign(R2)*sqrt(R2^2 - y1.^2) + xR2 
    else
        x2 = zeros(y1) + xV2
    end

    # x coordinates of lens top and bottom edges
    xe = zeros(2);
    if !isinf(R1) 
        xe[1] = -sign(R1)*sqrt(R1^2 - d^2/4) + xR1
    else
        xe[1] = xV1
    end
    if !isinf(R2) 
        xe[2] = -sign(R2)*sqrt(R2^2 - d^2/4) + xR2
    else
        xe[2] = xV2
    end

    xs = [x1; xe; reverse(x2); reverse(xe)] + x
    ys = [y1; [d/2, d/2]; reverse(y1); [-d/2, -d/2]] + y
    return (xs, ys) 
end # function draw


"""
    draw(l::SphericalLens, ny=100)

Spherical lens surface coordinates for plotting
"""
function draw(l::SphericalLens, ny::Integer=100)
    return draw(l.R1, l.R2, l.d, l.t, l.xV1, l.yR, ny)
end # function draw


"""
    draw(a::CircleAperture, ny=100)

Circular aperture surface coordinates for plotting
"""
function draw(a::CircleAperture, ny::Integer=100)
    xs = a.x*ones(div(ny, 2))
    ys1 = [linspace(-a.d2/2, -a.d1/2, div(ny, 2))...]
    ys2 = [linspace(a.d2/2, a.d1/2, div(ny, 2))...]
    ys1 += a.y
    ys2 += a.y
    return ((xs, ys1), (xs, ys2))
end # function draw 


"""
    draw(os::OpticalSystem, ny=100)

Array of (x array, y array) tuples with surface coordinates for each element in the optics system os for plotting
"""
function draw(os::OpticsSystem, ny::Real=100)
    xys = Array{Tuple{Array{Number,1}, Array{Number,1}},1}()
    for o in os
        xy = draw(o, ny)
        push!(xys, xy)
    end
    return xys
end # function draw


"""
    pyplot(ax::PyCall.PyObject, os::OpticsSystem, ny=100)

Patches for each optic in the optics system os are added to the PyPlot axis ax
"""
function pyplot(ax::PyCall.PyObject, os::OpticsSystem, ny=100)
    
    ps = Array{PyCall.PyObject,1}(length(os))
    for o in os
        if isa(o, Lens)
            xy = draw(o)
            p = PyPlot.matplotlib[:patches][:Polygon]([xy[1] xy[2]], closed=true, alpha=0.25, linestyle="None");
            push!(ps, p)
            ax[:add_patch](p)
        elseif isa(o, Aperture)
            xy = draw(o)
            ax[:plot](xy[1][1], xy[1][2], "-k", lw=1)
            ax[:plot](xy[2][1], xy[2][2], "-k", lw=1)
        end
    end 
end


"""
CIE 1964 10-deg colour matching function xbar, ybar, zbar values at 5 nm intervals from 360 nm to 830 nm.

Source: http://cvrl.ioo.ucl.ac.uk/cie.htm
"""
const cie_colour_match = Float32[
0.000129900000 0.000003917000 0.000606100000;
0.000232100000 0.000006965000 0.001086000000;
0.000414900000 0.000012390000 0.001946000000;
0.000741600000 0.000022020000 0.003486000000;
0.001368000000 0.000039000000 0.006450001000;
0.002236000000 0.000064000000 0.010549990000;
0.004243000000 0.000120000000 0.020050010000;
0.007650000000 0.000217000000 0.036210000000;
0.014310000000 0.000396000000 0.067850010000;
0.023190000000 0.000640000000 0.110200000000;
0.043510000000 0.001210000000 0.207400000000;
0.077630000000 0.002180000000 0.371300000000;
0.134380000000 0.004000000000 0.645600000000;
0.214770000000 0.007300000000 1.039050100000;
0.283900000000 0.011600000000 1.385600000000;
0.328500000000 0.016840000000 1.622960000000;
0.348280000000 0.023000000000 1.747060000000;
0.348060000000 0.029800000000 1.782600000000;
0.336200000000 0.038000000000 1.772110000000;
0.318700000000 0.048000000000 1.744100000000;
0.290800000000 0.060000000000 1.669200000000;
0.251100000000 0.073900000000 1.528100000000;
0.195360000000 0.090980000000 1.287640000000;
0.142100000000 0.112600000000 1.041900000000;
0.095640000000 0.139020000000 0.812950100000;
0.057950010000 0.169300000000 0.616200000000;
0.032010000000 0.208020000000 0.465180000000;
0.014700000000 0.258600000000 0.353300000000;
0.004900000000 0.323000000000 0.272000000000;
0.002400000000 0.407300000000 0.212300000000;
0.009300000000 0.503000000000 0.158200000000;
0.029100000000 0.608200000000 0.111700000000;
0.063270000000 0.710000000000 0.078249990000;
0.109600000000 0.793200000000 0.057250010000;
0.165500000000 0.862000000000 0.042160000000;
0.225749900000 0.914850100000 0.029840000000;
0.290400000000 0.954000000000 0.020300000000;
0.359700000000 0.980300000000 0.013400000000;
0.433449900000 0.994950100000 0.008749999000;
0.512050100000 1.000000000000 0.005749999000;
0.594500000000 0.995000000000 0.003900000000;
0.678400000000 0.978600000000 0.002749999000;
0.762100000000 0.952000000000 0.002100000000;
0.842500000000 0.915400000000 0.001800000000;
0.916300000000 0.870000000000 0.001650001000;
0.978600000000 0.816300000000 0.001400000000;
1.026300000000 0.757000000000 0.001100000000;
1.056700000000 0.694900000000 0.001000000000;
1.062200000000 0.631000000000 0.000800000000;
1.045600000000 0.566800000000 0.000600000000;
1.002600000000 0.503000000000 0.000340000000;
0.938400000000 0.441200000000 0.000240000000;
0.854449900000 0.381000000000 0.000190000000;
0.751400000000 0.321000000000 0.000100000000;
0.642400000000 0.265000000000 0.000049999990;
0.541900000000 0.217000000000 0.000030000000;
0.447900000000 0.175000000000 0.000020000000;
0.360800000000 0.138200000000 0.000010000000;
0.283500000000 0.107000000000 0.000000000000;
0.218700000000 0.081600000000 0.000000000000;
0.164900000000 0.061000000000 0.000000000000;
0.121200000000 0.044580000000 0.000000000000;
0.087400000000 0.032000000000 0.000000000000;
0.063600000000 0.023200000000 0.000000000000;
0.046770000000 0.017000000000 0.000000000000;
0.032900000000 0.011920000000 0.000000000000;
0.022700000000 0.008210000000 0.000000000000;
0.015840000000 0.005723000000 0.000000000000;
0.011359160000 0.004102000000 0.000000000000;
0.008110916000 0.002929000000 0.000000000000;
0.005790346000 0.002091000000 0.000000000000;
0.004109457000 0.001484000000 0.000000000000;
0.002899327000 0.001047000000 0.000000000000;
0.002049190000 0.000740000000 0.000000000000;
0.001439971000 0.000520000000 0.000000000000;
0.000999949300 0.000361100000 0.000000000000;
0.000690078600 0.000249200000 0.000000000000;
0.000476021300 0.000171900000 0.000000000000;
0.000332301100 0.000120000000 0.000000000000;
0.000234826100 0.000084800000 0.000000000000;
0.000166150500 0.000060000000 0.000000000000;
0.000117413000 0.000042400000 0.000000000000;
0.000083075270 0.000030000000 0.000000000000;
0.000058706520 0.000021200000 0.000000000000;
0.000041509940 0.000014990000 0.000000000000;
0.000029353260 0.000010600000 0.000000000000;
0.000020673830 0.000007465700 0.000000000000;
0.000014559770 0.000005257800 0.000000000000;
0.000010253980 0.000003702900 0.000000000000;
0.000007221456 0.000002607800 0.000000000000;
0.000005085868 0.000001836600 0.000000000000;
0.000003581652 0.000001293400 0.000000000000;
0.000002522525 0.000000910930 0.000000000000;
0.000001776509 0.000000641530 0.000000000000;
0.000001251141 0.000000451810 0.000000000000;
];


"""
Stiles and Burch (1959) 10-deg colour matching function values from 390 nm to 830 nm at 5 nm intervals

Source: http://cvrl.ioo.ucl.ac.uk/cmfs.htm
"""
const stiles_burch_10deg = Float32[
1.5000E-03 -4.0000E-04 6.2000E-03;
3.8000E-03 -1.0000E-03 1.6100E-02;
8.9000E-03 -2.5000E-03 4.0000E-02;
1.8800E-02 -5.9000E-03 9.0600E-02;
3.5000E-02 -1.1900E-02 1.8020E-01;
5.3100E-02 -2.0100E-02 3.0880E-01;
7.0200E-02 -2.8900E-02 4.6700E-01;
7.6300E-02 -3.3800E-02 6.1520E-01;
7.4500E-02 -3.4900E-02 7.6380E-01;
5.6100E-02 -2.7600E-02 8.7780E-01;
3.2300E-02 -1.6900E-02 9.7550E-01;
-4.4000E-03 2.4000E-03 1.0019E+00;
-4.7800E-02 2.8300E-02 9.9960E-01;
-9.7000E-02 6.3600E-02 9.1390E-01;
-1.5860E-01 1.0820E-01 8.2970E-01;
-2.2350E-01 1.6170E-01 7.4170E-01;
-2.8480E-01 2.2010E-01 6.1340E-01;
-3.3460E-01 2.7960E-01 4.7200E-01;
-3.7760E-01 3.4280E-01 3.4950E-01;
-4.1360E-01 4.0860E-01 2.5640E-01;
-4.3170E-01 4.7160E-01 1.8190E-01;
-4.4520E-01 5.4910E-01 1.3070E-01;
-4.3500E-01 6.2600E-01 9.1000E-02;
-4.1400E-01 7.0970E-01 5.8000E-02;
-3.6730E-01 7.9350E-01 3.5700E-02;
-2.8450E-01 8.7150E-01 2.0000E-02;
-1.8550E-01 9.4770E-01 9.5000E-03;
-4.3500E-02 9.9450E-01 7.0000E-04;
1.2700E-01 1.0203E+00 -4.3000E-03;
3.1290E-01 1.0375E+00 -6.4000E-03;
5.3620E-01 1.0517E+00 -8.2000E-03;
7.7220E-01 1.0390E+00 -9.4000E-03;
1.0059E+00 1.0029E+00 -9.7000E-03;
1.2710E+00 9.6980E-01 -9.7000E-03;
1.5574E+00 9.1620E-01 -9.3000E-03;
1.8465E+00 8.5710E-01 -8.7000E-03;
2.1511E+00 7.8230E-01 -8.0000E-03;
2.4250E+00 6.9530E-01 -7.3000E-03;
2.6574E+00 5.9660E-01 -6.3000E-03;
2.9151E+00 5.0630E-01 -5.3700E-03;
3.0779E+00 4.2030E-01 -4.4500E-03;
3.1613E+00 3.3600E-01 -3.5700E-03;
3.1673E+00 2.5910E-01 -2.7700E-03;
3.1048E+00 1.9170E-01 -2.0800E-03;
2.9462E+00 1.3670E-01 -1.5000E-03;
2.7194E+00 9.3800E-02 -1.0300E-03;
2.4526E+00 6.1100E-02 -6.8000E-04;
2.1700E+00 3.7100E-02 -4.4200E-04;
1.8358E+00 2.1500E-02 -2.7200E-04;
1.5179E+00 1.1200E-02 -1.4100E-04;
1.2428E+00 4.4000E-03 -5.4900E-05;
1.0070E+00 7.8000E-05 -2.2000E-06;
7.8270E-01 -1.3680E-03 2.3700E-05;
5.9340E-01 -1.9880E-03 2.8600E-05;
4.4420E-01 -2.1680E-03 2.6100E-05;
3.2830E-01 -2.0060E-03 2.2500E-05;
2.3940E-01 -1.6420E-03 1.8200E-05;
1.7220E-01 -1.2720E-03 1.3900E-05;
1.2210E-01 -9.4700E-04 1.0300E-05;
8.5300E-02 -6.8300E-04 7.3800E-06;
5.8600E-02 -4.7800E-04 5.2200E-06;
4.0800E-02 -3.3700E-04 3.6700E-06;
2.8400E-02 -2.3500E-04 2.5600E-06;
1.9700E-02 -1.6300E-04 1.7600E-06;
1.3500E-02 -1.1100E-04 1.2000E-06;
9.2400E-03 -7.4800E-05 8.1700E-07;
6.3800E-03 -5.0800E-05 5.5500E-07;
4.4100E-03 -3.4400E-05 3.7500E-07;
3.0700E-03 -2.3400E-05 2.5400E-07;
2.1400E-03 -1.5900E-05 1.7100E-07;
1.4900E-03 -1.0700E-05 1.1600E-07;
1.0500E-03 -7.2300E-06 7.8500E-08;
7.3900E-04 -4.8700E-06 5.3100E-08;
5.2300E-04 -3.2900E-06 3.6000E-08;
3.7200E-04 -2.2200E-06 2.4400E-08;
2.6500E-04 -1.5000E-06 1.6500E-08;
1.9000E-04 -1.0200E-06 1.1200E-08;
1.3600E-04 -6.8800E-07 7.5300E-09;
9.8400E-05 -4.6500E-07 5.0700E-09;
7.1300E-05 -3.1200E-07 3.4000E-09;
5.1800E-05 -2.0800E-07 2.2700E-09;
3.7700E-05 -1.3700E-07 1.5000E-09;
2.7600E-05 -8.8000E-08 9.8600E-10;
2.0300E-05 -5.5300E-08 6.3900E-10;
1.4900E-05 -3.3600E-08 4.0700E-10;
1.1000E-05 -1.9600E-08 2.5300E-10;
8.1800E-06 -1.0900E-08 1.5200E-10;
6.0900E-06 -5.7000E-09 8.6400E-11;
4.5500E-06 -2.7700E-09 4.4200E-11;
]


function XYZ2xyz(XYZ)
    if XYZ == zeros(XYZ)
        xyz = ones(3,1)/3
        return xyz
    else
        xyz = XYZ./+(XYZ...)
        return xyz
    end
end


function XYZ2RGB(XYZ)
    MXYZ2RGB = [ 3.240479 -1.537150 -0.498535;
                -0.969256  1.875992  0.041556;
                 0.055648 -0.204043  1.057311;]
    RGB = MXYZ2RGB*XYZ
    if minimum(RGB) < 0
        RGB -= minimum(RGB)
    end
    if maximum(RGB) > 1
        RGB /= maximum(RGB)
    end
    return RGB
end


function RGB2XYZ(RGB)
    MRGB2XYZ = [ 0.412453  0.357580  0.180423;
                 0.212671  0.715160  0.072169;
                 0.019334  0.119193  0.950227;]
    XYZ = MRGB2XYZ*RGB
    return XYZ
end


function RGB2sRGB(RGB)
    ind0 = RGB .<= 0.0031308 
    ind1 = RGB .>  0.0031308 
    sRGB = zeros(RGB)
    sRGB[ind0] = 12.92*RGB[ind0]
    sRGB[ind1] = 1.055*RGB[ind1].^(1/2.4) - 0.055
    return sRGB
end


function XYZ2sRGB(XYZ)
    RGB = XYZ2RGB(XYZ)
    sRGB = RGB2sRGB(RGB)
    return sRGB
end


"""
    spectrum2XYZ(wvl, val)

Convert a spectal power distribution to CIE 1931 XYZ coordinates
* wvl::Real: Wavelength in meters
* val::Real: Absolute or relative spectral power

# Examples
```jldoctest
julia> spectrum2XYZ([400e-9, 500e-9, 600e-9], [1, 3, 2])
3-element Array{Float32,1}:
 10.767
 11.157
  4.42725

```
"""
function spectrum2XYZ(wvl, val)
    if !isa(wvl, AbstractArray) & isa(wvl, Number)
        wvl = [wvl]
    end

    if !isa(val, AbstractArray) & isa(val, Number)
        val = [val]
    end

    @assert size(wvl) == size(val) "size(wvl) != size(val)"

    ind = Array{Int}(div(wvl.*1e9 - 360, 5) + 1)
    indOutMinRange = ind .< 1 
    indOutMaxRange = ind .> size(cie_colour_match)[1]
    indInRange = !(indOutMinRange | indOutMaxRange)
    ind = ind[indInRange]
    XYZ = transpose(5*transpose(val[indInRange])*cie_colour_match[ind,:]) + zeros(3,1)
    return XYZ
end # function spectrum2xyz


"""
    spectrum2RGB(wvl, val)
"""
function spectrum2RGB(wvl, val)
    XYZ = spectrum2XYZ(wvl, val)
    RGB = XYZ2RGB(XYZ)
    return RGB
end


"""
    spectrum2sRGB(wvl, val)
"""
function spectrum2sRGB(wvl, val)
    RGB = spectrum2RGB(wvl, val)
    sRGB = RGB2sRGB(RGB)
    return sRGB
end


"""
    wvl2RGBStilesBurch(wvl, val)

Convert a spectral power distribution to RGB coordinates
* wvl::Real: Wavelength in meters
* val::Real: Absolute or relative spectral power

# Examples
```jldoctest
julia> spectrum2RGB([400e-9, 500e-9, 600e-9], [1, 3, 2])
3-element Array{Float32,1}:
 1.0
 0.440969
 0.0

```
"""
function spectrum2RGBStilesBurch(wvl, val) 
    if !isa(wvl, AbstractArray) & isa(wvl, Number)
        wvl = [wvl]
    end

    if !isa(val, AbstractArray) & isa(val, Number)
        val = [val]
    end

    @assert size(wvl) == size(val) "size(wvl) != size(val)"

    ind = Array{Int}(div(wvl.*1e9 - 390, 5) + 1)
    RGB = 5*transpose(val)*stiles_burch_10deg[ind,:]
    RGB -= minimum(RGB)
    RGB /= maximum(RGB)
    R, G, B = RGB
    return [R, G, B]
end


end # module Optics
