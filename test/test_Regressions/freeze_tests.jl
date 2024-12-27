using LinearAlgebra

f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

# hdipole=HertzDipole([0.0, 0.0, 0.0], complex.([0.0, 0.0, 1.0]),complex(1.0))
# fdipole=FitzgeraldDipole([0.0, 0.0, 0.0], complex.([0.0, 1.0, 0.0]),complex(1.0))

hdipole = HertzArray([[0.0, 0.0, 0.0]], [complex.([0.0, 0.0, 1.0])], [complex(1.0)], k0)
fdipole =
    FitzgeraldArray([[0.0, 0.0, 0.0]], [complex.([0.0, 1.0, 0.0])], [complex(1.0)], k0)
mixedarray = [hdipole, fdipole]

refvalue = zeros(ComplexF64, 30)
refvalue[4] = -140.54491802563436 - 0.0im
# hspherical=convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [hdipole], k0)
# @test hspherical.coefficients ≈ refvalue
hspherical = SphericalWaveExpansion(Radiated(), refvalue, k0)

refvalue = zeros(ComplexF64, 30)
refvalue[1] = -0.2637968355397918 + 0.0im
refvalue[5] = -0.2637968355397918 + 0.0im
# fspherical=convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [fdipole], k0)
# @test fspherical.coefficients ≈ refvalue
fspherical = SphericalWaveExpansion(Radiated(), refvalue, k0)

R = [1.0, 2.0, 3.0]

refvalue = [
    -53.33289214120721 - 8.34494227296252im
    -106.66578428241442 - 16.68988454592504im
    88.1542208957773 + 18.127788150666984im
]
@test efield(hdipole, R) ≈ refvalue
@test efield(hspherical, R) ≈ refvalue

refvalue = [
    0.35211627648388694 + 0.061245767189904im
    -0.17605813824194347 - 0.030622883594952im
    -0.0 - 0.0im
]
@test hfield(hdipole, R) ≈ refvalue
@test hfield(hspherical, R) ≈ refvalue

refvalue = [
    0.5281744147258304 + 0.09186865078485598im
    0.0 + 0.0im
    -0.17605813824194347 - 0.030622883594952im
]
@test efield(fdipole, R) ≈ refvalue
@test efield(fspherical, R) ≈ refvalue

refvalue = [
    -0.0002505201954500046 - 3.919863493819998e-5im
    0.0012474297325655351 + 0.000225723872935714im
    -0.0007515605863500138 - 0.00011759590481459994im
]
@test hfield(fdipole, R) ≈ refvalue
@test hfield(fspherical, R) ≈ refvalue



dR = [0.1, 0.2, 0.3]
shifted_hdipole =
    HertzArray([[0.0, 0.0, 0.0] + dR], [complex.([0.0, 0.0, 1.0])], [complex(1.0)], k0)
shifted_hspherical = changerepresentation(SphericalWaveExpansion{Radiated}, shifted_hdipole)

# ff = shiftrepresentation(
#     convertrepresentation(FarfieldPattern{ComplexF64}, hspherical, k0),
#     -dR,
#     k0,
# )

# shifted_hspherical2 = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, ff, k0)

# @test shifted_hspherical.coefficients ≈ shifted_hspherical2.coefficients
@test efield(shifted_hspherical, R + dR) ≈
      efield(shifted_hdipole, R + dR) ≈
      efield(hdipole, R)

ϑ = pi / 7
φ = 5pi / 3
refvalue = -207.9842114958308 + 352.0836119998865im
Eθ, Eϕ = farfield(shifted_hdipole, (ϑ, φ))
@test Eθ ≈ refvalue
Eθ, Eϕ = farfield(shifted_hspherical, (ϑ, φ))
@test Eθ ≈ refvalue
