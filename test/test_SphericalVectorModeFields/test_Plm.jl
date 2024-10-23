# Check against Wikipedia

import AntennaFieldRepresentations.Plm_deriv_array
Nmax = 4
ϑ = 0.2
ct = cos(ϑ)
st = sin(ϑ)


m = 0
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 1.0
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = cos(ϑ)
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 0.5 * (3 * (ct^2) - 1)
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = 0.5 * (5 * ct^3 - 3 * ct)
@test abs(P1[4] - goal) / abs(goal) < 1e-10

goal = 1 / 8 * (35 * ct^4 - 30ct^2 + 3)
@test abs(P1[5] - goal) / abs(goal) < 1e-10

m = 1
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -st
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = -3 * ct * st
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = -3 / 2 * (5 * ct^2 - 1) * st
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = -5 / 2 * (7ct^3 - 3ct) * st
@test abs(P1[4] - goal) / abs(goal) < 1e-10

m = 2
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 3 * st^2
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = 15 * ct * st^2
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 15 / 2 * (7 * ct^2 - 1) * st^2
@test abs(P1[3] - goal) / abs(goal) < 1e-10

m = 3
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -15 * st^3
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = -105 * ct * st^3
@test abs(P1[2] - goal) / abs(goal) < 1e-10

m = 4
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 105 * st^4
@test abs(P1[1] - goal) / abs(goal) < 1e-10


Nmax = 4
ϑ = pi / 2
ct = cos(ϑ)
st = sin(ϑ)


m = 0
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 1.0
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = cos(ϑ)
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 0.5 * (3 * (ct^2) - 1)
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = 0.5 * (5 * ct^3 - 3 * ct)
@test abs(P1[4] - goal) / abs(goal) < 1e-10

goal = 1 / 8 * (35 * ct^4 - 30ct^2 + 3)
@test abs(P1[5] - goal) / abs(goal) < 1e-10

m = 1
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -st
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = -3 * ct * st
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = -3 / 2 * (5 * ct^2 - 1) * st
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = -5 / 2 * (7ct^3 - 3ct) * st
@test abs(P1[4] - goal) / abs(goal) < 1e-10

m = 2
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 3 * st^2
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = 15 * ct * st^2
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 15 / 2 * (7 * ct^2 - 1) * st^2
@test abs(P1[3] - goal) / abs(goal) < 1e-10

m = 3
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -15 * st^3
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = -105 * ct * st^3
@test abs(P1[2] - goal) / abs(goal) < 1e-10

m = 4
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 105 * st^4
@test abs(P1[1] - goal) / abs(goal) < 1e-10

Nmax = 4
ϑ = 0.0
ct = cos(ϑ)
st = sin(ϑ)


m = 0
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 1.0
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = cos(ϑ)
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 0.5 * (3 * (ct^2) - 1)
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = 0.5 * (5 * ct^3 - 3 * ct)
@test abs(P1[4] - goal) / abs(goal) < 1e-10

goal = 1 / 8 * (35 * ct^4 - 30ct^2 + 3)
@test abs(P1[5] - goal) / abs(goal) < 1e-10

m = 1
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -st
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = -3 * ct * st
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

goal = -3 / 2 * (5 * ct^2 - 1) * st
@test abs(P1[3] - goal) / abs(goal) < 1e-10 || (abs(P1[3] - goal) < 1e-16)

goal = -5 / 2 * (7ct^3 - 3ct) * st
@test abs(P1[4] - goal) / abs(goal) < 1e-10 || (abs(P1[4] - goal) < 1e-16)

m = 2
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 3 * st^2
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = 15 * ct * st^2
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

goal = 15 / 2 * (7 * ct^2 - 1) * st^2
@test abs(P1[3] - goal) / abs(goal) < 1e-10 || (abs(P1[3] - goal) < 1e-16)

m = 3
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -15 * st^3
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = -105 * ct * st^3
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

m = 4
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 105 * st^4
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)


Nmax = 4
ϑ = pi
ct = cos(big(ϑ))
st = sin(big(ϑ))


m = 0
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 1.0
@test abs(P1[1] - goal) / abs(goal) < 1e-10

goal = cos(ϑ)
@test abs(P1[2] - goal) / abs(goal) < 1e-10

goal = 0.5 * (3 * (ct^2) - 1)
@test abs(P1[3] - goal) / abs(goal) < 1e-10

goal = 0.5 * (5 * ct^3 - 3 * ct)
@test abs(P1[4] - goal) / abs(goal) < 1e-10

goal = 1 / 8 * (35 * ct^4 - 30ct^2 + 3)
@test abs(P1[5] - goal) / abs(goal) < 1e-10

m = 1
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -st
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = -3 * ct * st
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

goal = -3 / 2 * (5 * ct^2 - 1) * st
@test abs(P1[3] - goal) / abs(goal) < 1e-10 || (abs(P1[3] - goal) < 1e-16)

goal = -5 / 2 * (7ct^3 - 3ct) * st
@test abs(P1[4] - goal) / abs(goal) < 1e-10 || (abs(P1[4] - goal) < 1e-16)

m = 2
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 3 * st^2
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = 15 * ct * st^2
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

goal = 15 / 2 * (7 * ct^2 - 1) * st^2
@test abs(P1[3] - goal) / abs(goal) < 1e-10 || (abs(P1[3] - goal) < 1e-16)

m = 3
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = -15 * st^3
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)

goal = -105 * ct * st^3
@test abs(P1[2] - goal) / abs(goal) < 1e-10 || (abs(P1[2] - goal) < 1e-16)

m = 4
P1, P2 = Plm_deriv_array(Nmax, m, ct)

goal = 105 * st^4
@test abs(P1[1] - goal) / abs(goal) < 1e-10 || (abs(P1[1] - goal) < 1e-16)
