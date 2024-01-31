# Dipole to Dipole Interactions
The interaction between two elementary dipoles is modeled by letting one of the dipoles generate an electromagnetic field, which is therafter received by the other dipole.

Since the recieved signal of a Hertzian dipole is directly proportional to the parallel component of the incident electric field at the dipole position, the received signals of Hertzian dipoles can be used to ideally probe the electric field. The interaction between two Hertzian dipoles can be used to represent the radiated electric field of one of the dipoles at the location of the other dipole.

Since the recieved signal of a Fitzgerald dipole is directly proportional to the parallel component of the negative incident magnetic field at the dipole position, the received signals of Fitzgerald  dipoles can be used to ideally probe the magnetic field. The interaction between two Fitzgerald dipoles can be used to represent the radiated magnetic field of one of the dipoles at the location of the other dipole.

The detailed derivation can be found in the "Theory" section.
The interaction between two dipoles is implemented in the function
(TODO: insert correct function)
Each elementary dipole is defined by its position ``\bm{r}``, its oriented length ``\bm{\ell}`` and its excitation magnitude ``I``.
Hertzian and Fitzgerald dipoles may be constructed using the constructors as in the following example.


## Interaction between two Hertzian Dipoles
Let ``I_1``, ``\bm{\ell}_1``, and ``\bm{r}_1`` denote the excitation, oriented (equivalent) length, and location of the transmitting Hertzian dipole and let ``I_2``, ``\bm{\ell}_2``, and ``\bm{r}_2`` be the excitation, oriented (equivalent) length, and location of the receiving Hertzian dipole. The received signal which represents the interaction between these two dipoles is given by
```math
b=
\dfrac{-\mathrm{j}}{2} k_0\, Z_\mathrm{F}\, I_1\, I_2\,
\left[
\left(
\dfrac{3}{k^2\,||\bm{r}_2-\bm{r}_1||^2}
+
\dfrac{3\mathrm{j}}{k\,||\bm{r}_2-\bm{r}_1||} 
-1   
\right)
g_0(\bm{r}_2,\bm{r}_1)
\,
\left(\bm{\ell}_2 \cdot \bm{e}_r\right)\, \left(\bm{e}_r \cdot \bm{\ell}_1\right)
-
\left(
\dfrac{\mathrm{j}}{k\,||\bm{r}_2-\bm{r}_1||}
+\dfrac{1}{k^2\,||\bm{r}_2-\bm{r}_1||^2}
-1
\right)
g_0(\bm{r}_2,\bm{r}_1)
\,
\left(\bm{\ell_2} \cdot \bm{\ell_1}
\right)
\right]\, ,
```
where ``\bm{e}_r`` denotes the unit vector pointing in the direction from one dipole location to the other.


## Interaction between two Fitzgerald Dipoles
Let ``I_{\mathrm{m},1}``, ``\bm{\ell}_1``, and ``\bm{r}_1`` denote the excitation, oriented (equivalent) length, and location of the transmitting Fitzgerald dipole and let ``I_{\mathrm{m},2}``, ``\bm{\ell}_2``, and ``\bm{r}_2`` be the excitation, oriented (equivalent) length, and location of the receiving Fitzgerald dipole. The received signal which represents the interaction between these two dipoles is given by
```math
b=
\dfrac{-\mathrm{j}}{2} \dfrac{k_0}{Z_\mathrm{F}}\, I_{\mathrm{m},1}\, I_{\mathrm{m},2}\,
\left[
\left(
\dfrac{3}{k^2\,||\bm{r}_2-\bm{r}_1||^2}
+
\dfrac{3\mathrm{j}}{k\,||\bm{r}_2-\bm{r}_1||} 
-1   
\right)
g_0(\bm{r}_2,\bm{r}_1)
\,
\left(\bm{\ell}_2 \cdot \bm{e}_r\right)\, \left(\bm{e}_r \cdot \bm{\ell}_1\right)
-
\left(
\dfrac{\mathrm{j}}{k\,||\bm{r}_2-\bm{r}_1||}
+\dfrac{1}{k^2\,||\bm{r}_2-\bm{r}_1||^2}
-1
\right)
g_0(\bm{r}_2,\bm{r}_1)
\,
\left(\bm{\ell_2} \cdot \bm{\ell_1}
\right)
\right]\, ,
```
where ``\bm{e}_r`` denotes the unit vector pointing in the direction from one dipole location to the other.

## Interaction between Hertzian and Fitzgerald Dipole
Let ``I_{1}``, ``\bm{\ell}_1``, and ``\bm{r}_1`` denote the excitation, oriented (equivalent) length, and location of the transmitting Hertzian dipole and let ``I_{\mathrm{m},2}``, ``\bm{\ell}_2``, and ``\bm{r}_2`` be the excitation, oriented (equivalent) length, and location of the receiving Fitzgerald dipole. The received signal which represents the interaction between these two dipoles is given by
```math
b= \dfrac{1}{2}
I_1\, I_{\mathrm{m},2}
\left(-\mathrm{j}k-\dfrac{1}{||\bm{r}-\bm{r}'||}\right) \bm{\ell_2} \cdot \left( \bm{e}_r \times \bm{\ell_1}\right)
```
where ``\bm{e}_r`` denotes the unit vector pointing in the direction from one dipole location to the other.


The role of the transmitting and the receiving dipole may be exchanged without changing the result of the calculation (i.e., the electric field of a Fitzgerald dipole matches the electric field of a Fitzgerald dipole).

