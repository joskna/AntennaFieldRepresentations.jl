ϑlist = [pi / 6, pi / 2, 2 / 3 * pi]
import AntennaFieldRepresentations.legendre_deps_array

for k = 1:length(ϑlist)
    global m, n, P1, P2, P3, ϑ, Nmax, goal
    # Check against Hansen p. 322 ff
    Nmax = 5
    ϑ = ϑlist[k]

    m = 0
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 1
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -√6 / 2 * sin(ϑ)
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = √6 / 2 * cos(ϑ)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 2
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -3 * √10 / 4 * sin(2 * ϑ)
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = √10 / 8 * (3 * cos(2 * ϑ) + 1)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 3
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -3 * √14 / 16 * (5 * sin(3 * ϑ) + sin(ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = √14 / 16 * (5 * cos(3 * ϑ) + 3 * cos(ϑ))
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -15 * √2 / 32 * (7 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = 3 * √2 / 128 * (35 * cos(4 * ϑ) + 20 * cos(2 * ϑ) + 9)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 5
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -15 * √22 / 256 * (21 * sin(5 * ϑ) + 7 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = √22 / 256 * (63 * cos(5 * ϑ) + 35 * cos(3 * ϑ) + 30 * cos(ϑ))
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10

    m = 1
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 1
    goal = √3 / 2
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √3 / 2 * cos(ϑ)
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √3 / 2 * sin(ϑ)
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 2
    goal = √15 / 2 * (cos(ϑ))
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √15 / 2 * cos(2 * ϑ)
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √15 / 4 * (sin(2 * ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 3
    goal = √42 / 16 * (5 * cos(2 * ϑ) + 3)
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √42 / 32 * (15 * cos(3 * ϑ) + cos(ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √42 / 32 * (5 * sin(3 * ϑ) + sin(ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 3 * √10 / 32 * (7 * cos(3 * ϑ) + 9 * cos(ϑ))
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = 3 * √10 / 16 * (7 * cos(4 * ϑ) + cos(2 * ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = 3 * √10 / 64 * (7 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 5
    goal = √165 / 128 * (21 * cos(4 * ϑ) + 28 * cos(2 * ϑ) + 15)
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √165 / 256 * (105 * cos(5 * ϑ) + 21 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √165 / 256 * (21 * sin(5 * ϑ) + 7 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10

    m = 2
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 2
    goal = √15 / 2 * (sin(ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √15 / 4 * (sin(2 * ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√15 / 8 * (cos(2 * ϑ) - 1)
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 3
    goal = √105 / 4 * sin(2 * ϑ)
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √105 / 16 * (3 * sin(3 * ϑ) - sin(ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√105 / 16 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 3 * √5 / 16 * (7 * sin(3 * ϑ) + 3 * sin(ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = 3 * √5 / 16 * (7 * sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -3 * √5 / 64 * (7 * cos(4 * ϑ) - 4 * cos(2 * ϑ) - 3)
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 5
    goal = √1155 / 32 * (3 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √1155 / 128 * (15 * sin(5 * ϑ) - 3 * sin(3 * ϑ) - 2 * sin(ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√1155 / 128 * (3 * cos(5 * ϑ) - cos(3 * ϑ) - 2 * cos(ϑ))
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10

    m = 3
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 3
    goal = -3 * √70 / 16 * (cos(2 * ϑ) - 1)
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √70 / 32 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10
    goal = -√70 / 32 * (sin(3 * ϑ) - 3 * sin(ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10
    n = 4
    goal = -9 * √70 / 32 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √70 / 16 * (cos(4 * ϑ) - cos(2 * ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10 || abs(P2[n-2] - goal) < 1e-14
    goal = -3 * √70 / 64 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10
    n = 5
    goal = -3 * √770 / 256 * (9 * cos(4 * ϑ) - 4 * cos(2 * ϑ) - 5)
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √770 / 512 * (15 * cos(5 * ϑ) - 13 * cos(3 * ϑ) - 2 * cos(ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10
    goal = -√770 / 512 * (9 * sin(5 * ϑ) - 13 * sin(3 * ϑ) - 6 * sin(ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10

    m = 4
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 4
    goal = -3 * √35 / 16 * (sin(3 * ϑ) - 3 * sin(ϑ))
    @test abs(P1[n-3] - goal) / abs(goal) < 1e-10
    goal = -3 * √35 / 32 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P2[n-3] - goal) / abs(goal) < 1e-10
    goal = 3 * √35 / 128 * (cos(4 * ϑ) - 4 * cos(2 * ϑ) + 3)
    @test abs(P3[n-3] - goal) / abs(goal) < 1e-10
    n = 5
    goal = -3 * √385 / 32 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P1[n-3] - goal) / abs(goal) < 1e-10
    goal = -3 * √385 / 256 * (5 * sin(5 * ϑ) - 9 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P2[n-3] - goal) / abs(goal) < 1e-10
    goal = 3 * √385 / 256 * (cos(5 * ϑ) - 3 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P3[n-3] - goal) / abs(goal) < 1e-10

    m = 5
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 5
    goal = 15 * √154 / 256 * (cos(4 * ϑ) - 4 * cos(2 * ϑ) + 3)
    @test abs(P1[n-4] - goal) / abs(goal) < 1e-10
    goal = 15 * √154 / 512 * (cos(5ϑ) - 3 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P2[n-4] - goal) / abs(goal) < 1e-10
    goal = 3 * √154 / 512 * (sin(5 * ϑ) - 5 * sin(3 * ϑ) + 10 * sin(ϑ))
    @test abs(P3[n-4] - goal) / abs(goal) < 1e-10
end

ϑlist = [big(1e-4), big((1 - 1e-4) * pi)]

for k = 1:length(ϑlist)
    global m, n, P1, P2, P3, ϑ, Nmax, goal
    # Check against Hansen p. 322 ff
    Nmax = 5
    ϑ = ϑlist[k]

    m = 0
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 1
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -√6 / 2 * sin(ϑ)
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10
    goal = √6 / 2 * cos(ϑ)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 2
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -3 * √10 / 4 * sin(2 * ϑ)
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10 || abs(P2[n+1] - goal) < 1e-11
    goal = √10 / 8 * (3 * cos(2 * ϑ) + 1)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 3
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -3 * √14 / 16 * (5 * sin(3 * ϑ) + sin(ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10 || abs(P2[n+1] - goal) < 1e-11
    goal = √14 / 16 * (5 * cos(3 * ϑ) + 3 * cos(ϑ))
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -15 * √2 / 32 * (7 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10 || abs(P2[n+1] - goal) < 1e-11
    goal = 3 * √2 / 128 * (35 * cos(4 * ϑ) + 20 * cos(2 * ϑ) + 9)
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10
    n = 5
    goal = 0.0
    @test abs(P1[n+1]) < 1e-10
    goal = -15 * √22 / 256 * (21 * sin(5 * ϑ) + 7 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P2[n+1] - goal) / abs(goal) < 1e-10 || abs(P2[n+1] - goal) < 1e-12
    goal = √22 / 256 * (63 * cos(5 * ϑ) + 35 * cos(3 * ϑ) + 30 * cos(ϑ))
    @test abs(P3[n+1] - goal) / abs(goal) < 1e-10

    m = 1
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 1
    goal = √3 / 2
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √3 / 2 * cos(ϑ)
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √3 / 2 * sin(ϑ)
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 2
    goal = √15 / 2 * (cos(ϑ))
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √15 / 2 * cos(2 * ϑ)
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √15 / 4 * (sin(2 * ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 3
    goal = √42 / 16 * (5 * cos(2 * ϑ) + 3)
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √42 / 32 * (15 * cos(3 * ϑ) + cos(ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √42 / 32 * (5 * sin(3 * ϑ) + sin(ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 3 * √10 / 32 * (7 * cos(3 * ϑ) + 9 * cos(ϑ))
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = 3 * √10 / 16 * (7 * cos(4 * ϑ) + cos(2 * ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = 3 * √10 / 64 * (7 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10
    n = 5
    goal = √165 / 128 * (21 * cos(4 * ϑ) + 28 * cos(2 * ϑ) + 15)
    @test abs(P1[n] - goal) / abs(goal) < 1e-10
    goal = √165 / 256 * (105 * cos(5 * ϑ) + 21 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P2[n] - goal) / abs(goal) < 1e-10
    goal = √165 / 256 * (21 * sin(5 * ϑ) + 7 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P3[n] - goal) / abs(goal) < 1e-10

    m = 2
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 2
    goal = √15 / 2 * (sin(ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √15 / 4 * (sin(2 * ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√15 / 8 * (cos(2 * ϑ) - 1)
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 3
    goal = √105 / 4 * sin(2 * ϑ)
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √105 / 16 * (3 * sin(3 * ϑ) - sin(ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√105 / 16 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 4
    goal = 3 * √5 / 16 * (7 * sin(3 * ϑ) + 3 * sin(ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = 3 * √5 / 16 * (7 * sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -3 * √5 / 64 * (7 * cos(4 * ϑ) - 4 * cos(2 * ϑ) - 3)
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10
    n = 5
    goal = √1155 / 32 * (3 * sin(4 * ϑ) + 2 * sin(2 * ϑ))
    @test abs(P1[n-1] - goal) / abs(goal) < 1e-10
    goal = √1155 / 128 * (15 * sin(5 * ϑ) - 3 * sin(3 * ϑ) - 2 * sin(ϑ))
    @test abs(P2[n-1] - goal) / abs(goal) < 1e-10
    goal = -√1155 / 128 * (3 * cos(5 * ϑ) - cos(3 * ϑ) - 2 * cos(ϑ))
    @test abs(P3[n-1] - goal) / abs(goal) < 1e-10

    m = 3
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 3
    goal = -3 * √70 / 16 * (cos(2 * ϑ) - 1)
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √70 / 32 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10
    goal = -√70 / 32 * (sin(3 * ϑ) - 3 * sin(ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10
    n = 4
    goal = -9 * √70 / 32 * (cos(3 * ϑ) - cos(ϑ))
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √70 / 16 * (cos(4 * ϑ) - cos(2 * ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10 || abs(P2[n-2] - goal) < 1e-14
    goal = -3 * √70 / 64 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10
    n = 5
    goal = -3 * √770 / 256 * (9 * cos(4 * ϑ) - 4 * cos(2 * ϑ) - 5)
    @test abs(P1[n-2] - goal) / abs(goal) < 1e-10
    goal = -3 * √770 / 512 * (15 * cos(5 * ϑ) - 13 * cos(3 * ϑ) - 2 * cos(ϑ))
    @test abs(P2[n-2] - goal) / abs(goal) < 1e-10
    goal = -√770 / 512 * (9 * sin(5 * ϑ) - 13 * sin(3 * ϑ) - 6 * sin(ϑ))
    @test abs(P3[n-2] - goal) / abs(goal) < 1e-10

    m = 4
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= 1 * (sqrt(2 * pi))
    P2 *= 1 * (sqrt(2 * pi))
    P3 *= 1 * (sqrt(2 * pi))

    n = 4
    goal = -3 * √35 / 16 * (sin(3 * ϑ) - 3 * sin(ϑ))
    @test abs(P1[n-3] - goal) / abs(goal) < 1e-10
    goal = -3 * √35 / 32 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P2[n-3] - goal) / abs(goal) < 1e-10
    goal = 3 * √35 / 128 * (cos(4 * ϑ) - 4 * cos(2 * ϑ) + 3)
    @test abs(P3[n-3] - goal) / abs(goal) < 1e-10
    n = 5
    goal = -3 * √385 / 32 * (sin(4 * ϑ) - 2 * sin(2 * ϑ))
    @test abs(P1[n-3] - goal) / abs(goal) < 1e-10
    goal = -3 * √385 / 256 * (5 * sin(5 * ϑ) - 9 * sin(3 * ϑ) + 2 * sin(ϑ))
    @test abs(P2[n-3] - goal) / abs(goal) < 1e-10
    goal = 3 * √385 / 256 * (cos(5 * ϑ) - 3 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P3[n-3] - goal) / abs(goal) < 1e-10

    m = 5
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)
    P1 *= -1 * (sqrt(2 * pi))
    P2 *= -1 * (sqrt(2 * pi))
    P3 *= -1 * (sqrt(2 * pi))

    n = 5
    goal = 15 * √154 / 256 * (cos(4 * ϑ) - 4 * cos(2 * ϑ) + 3)
    @test abs(P1[n-4] - goal) / abs(goal) < 1e-10
    goal = 15 * √154 / 512 * (cos(5ϑ) - 3 * cos(3 * ϑ) + 2 * cos(ϑ))
    @test abs(P2[n-4] - goal) / abs(goal) < 1e-10
    goal = 3 * √154 / 512 * (sin(5 * ϑ) - 5 * sin(3 * ϑ) + 10 * sin(ϑ))
    @test abs(P3[n-4] - goal) / abs(goal) < 1e-10
end






Nmax = 10
ϑ = 0.7im
arr = [5, 6, 7]
for mm in eachindex(arr)
    global m, n, P1, P2, P3, ϑ, Nmax, goal
    m = arr[mm]
    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)

    @test abs(P1[1]) < Inf
    @test abs(P1[2]) < Inf
    @test abs(P1[3]) < Inf

    @test abs(P2[1]) < Inf
    @test abs(P2[2]) < Inf
    @test abs(P2[3]) < Inf


    @test abs(P3[1]) < Inf
    @test abs(P3[2]) < Inf
    @test abs(P3[3]) < Inf

    # Nmax=5
    ϑ = 0.5 + 0.7im

    # m=3  

    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)

    @test abs(P1[1]) < Inf
    @test abs(P1[2]) < Inf
    @test abs(P1[3]) < Inf

    @test abs(P2[1]) < Inf
    @test abs(P2[2]) < Inf
    @test abs(P2[3]) < Inf


    @test abs(P3[1]) < Inf
    @test abs(P3[2]) < Inf
    @test abs(P3[3]) < Inf


    # Nmax=5
    ϑ = pi / 2 + 0.7im

    # m=3  

    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)

    @test abs(P1[1]) < Inf
    @test abs(P1[2]) < Inf
    @test abs(P1[3]) < Inf

    @test abs(P2[1]) < Inf
    @test abs(P2[2]) < Inf
    @test abs(P2[3]) < Inf


    @test abs(P3[1]) < Inf
    @test abs(P3[2]) < Inf
    @test abs(P3[3]) < Inf


    # Nmax=5
    ϑ = pi / 2 + 30.3im

    # m=3  

    P1, P2, P3 = legendre_deps_array(m, Nmax, ϑ)

    @test abs(P1[1]) < Inf
    @test abs(P1[2]) < Inf
    @test abs(P1[3]) < Inf

    @test abs(P2[1]) < Inf
    @test abs(P2[2]) < Inf
    @test abs(P2[3]) < Inf


    @test abs(P3[1]) < Inf
    @test abs(P3[2]) < Inf
    @test abs(P3[3]) < Inf

end


# using GSL


# using GSL
# import ElectromagneticFieldTransformations.sphPlm_deriv_array
# # Check agaainst GNU Scientific Library
# ϑlist=[pi/6, pi/2, 2/3*pi]
# Nmax=5
# for k=1:length(ϑlist)
#     global ϑ, P1, P2
#     ϑ=ϑlist[k]
#     P1,P2= sphPlm_deriv_array(Nmax,0, cos(ϑ))
#     P10,P20= sf_legendre_sphPlm_deriv_array(Nmax,0, cos(ϑ))
#     for ik=1:Nmax+1
#         @test abs(P1[ik]-P10[ik])/abs(P10[ik])< 1e-10
#         @test abs(P2[ik]-P20[ik])< 3e-15
#     end

#     for iik =1:Nmax
#         P1,P2= sphPlm_deriv_array(Nmax,iik, cos(ϑ))
#     P10,P20= sf_legendre_sphPlm_deriv_array(Nmax,iik, cos(ϑ))
#     for ik=1:Nmax-iik+1
#         @test abs(P1[ik]-P10[ik])/abs(P10[ik])< 1e-10
#         @test (abs(P2[ik]-P20[ik])/abs(P20[ik])< 1e-10) || (abs(P2[ik]-P20[ik])< 3e-15 && P20[ik]<1e-9)
#     end
#     end
# end

# ϑlist=[(1e-4), ((1-1e-4)*pi)]
# Nmax=5
# for k=1:length(ϑlist)
#     global ϑ, P1, P2
#     ϑ=ϑlist[k]
#     P1,P2= sphPlm_deriv_array(Nmax,0, cos(ϑ))
#     P10,P20= sf_legendre_sphPlm_deriv_array(Nmax,0, cos(ϑ))
#     for ik=1:Nmax+1
#         @test abs(P1[ik]-P10[ik])/abs(P10[ik])< 1e-8
#         @test (abs(P2[ik]-P20[ik])/abs(P20[ik])< 1e-8) || (abs(P2[ik]-P20[ik])<2e-7)
#     end

#     for iik =1:Nmax
#         P1,P2= sphPlm_deriv_array(Nmax,iik, cos(ϑ))
#     P10,P20= sf_legendre_sphPlm_deriv_array(Nmax,iik, cos(ϑ))
#     for ik=1:Nmax-iik+1
#         @test abs(P1[ik]-P10[ik])/abs(P10[ik])< 1e-8
#         @test (abs(P2[ik]-P20[ik])/abs(P20[ik])< 1e-8) 
#     end
#     end
# end
