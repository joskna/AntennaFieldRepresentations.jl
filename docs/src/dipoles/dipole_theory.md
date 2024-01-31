# Electromagnetic Fields of Elementary Dipoles


## Radiated Fields

### General Current Distribution
For a general source, consisting of an electric volume current distribution ``\bm{J}(\bm{r})`` and a magnetic current distribution ``\bm{M}(\bm{r})``, the radiated electric and magnetic fields can be described by 

```math
\bm{E}(\bm{r})=
\iiint \limits_{V_{\text{src}}}
\left[
\bm{\mathcal{G}}_{JE}
(\bm{r},\bm{r}')
\cdot
\bm{J}(\bm{r}')
+
\bm{\mathcal{G}}_{ME}
(\bm{r},\bm{r}')
\cdot
\bm{M}(\bm{r}')
\right]\, \mathrm{d}v' 
```
and
```math
\bm{H}(\bm{r})=
\iiint \limits_{V_{\text{src}}}
\left[
\bm{\mathcal{G}}_{JH}
(\bm{r},\bm{r}')
\cdot
\bm{J}(\bm{r}')
+
\bm{\mathcal{G}}_{MH}
(\bm{r},\bm{r}')
\cdot
\bm{M}(\bm{r}')
\right]\, \mathrm{d}v'\, , 
```
where 
```math
\bm{\mathcal{G}}_{JE}=
-\mathrm{j} k_0\, Z_\mathrm{F}
\left[
\left(\mathbf{I} +\dfrac{1}{k^2} \nabla \nabla\right) g_0(\bm{r},\bm{r}')
\right]
```
```math
\bm{\mathcal{G}}_{ME}=
-\nabla g_0(\bm{r},\bm{r}') \times \mathbf{I}
```
```math
\bm{\mathcal{G}}_{JH}=
\nabla g_0(\bm{r},\bm{r}') \times \mathbf{I}
```
```math
\bm{\mathcal{G}}_{MH}=
-\mathrm{j} \dfrac{k_0}{Z_\mathrm{F}}
\left[
\left(\mathbf{I} +\dfrac{1}{k^2} \nabla \nabla\right) g_0(\bm{r},\bm{r}')
\right]
```
are the dyadic Green's functions derived from the scalar greens function
```math
 g_0(\bm{r},\bm{r}')
=
\dfrac{\mathrm{e}^{-\mathrm{j}k_0\, ||\bm{r}-\bm{r}'||}}{4\pi\, ||\bm{r}-\bm{r}'||}\, .
``` 

The involved dyadic operators have the explicit expressions
```math
\left(\mathbf{I} +\dfrac{1}{k^2} \nabla \nabla\right) g_0(\bm{r},\bm{r}')
=
\left(
\dfrac{3}{k^2\,||\bm{r}-\bm{r}'||^2}
+
\dfrac{3\mathrm{j}}{k\,||\bm{r}-\bm{r}'||} 
-1   
\right)
g_0(\bm{r},\bm{r}')
\,
\bm{e}_r \bm{e}_r 
-
\left(
\dfrac{\mathrm{j}}{k\,||\bm{r}-\bm{r}'||}
+\dfrac{1}{k^2\,||\bm{r}-\bm{r}'||^2}
-1
\right)
g_0(\bm{r},\bm{r}')
\,
\mathbf{I}
```
and
```math
\nabla g_0(\bm{r},\bm{r}') \times \mathbf{I}
=
\left(-\mathrm{j}k-\dfrac{1}{||\bm{r}-\bm{r}'||}\right) \bm{e}_r \times \mathbf{I}\, .
```
### Hertzian Dipole
A Hertzian dipole is a model for an elementary electric dipole. Formally, it can be regardes as a special electric current distribution of the form
```math
\bm{J}_{\text{Hertz}}(\bm{r})
=
I\bm{\ell}
\,
\delta({\bm{r}-\bm{r}'})\, ,
```
where ``I`` has the unit Ampère and denotes the excitation of the dipole, ``\bm{\ell}`` denotes the oriented dipole length and ``\delta({\bm{r}-\bm{r}'})`` is the Dirac-Delta distribution. In the idealized model, ``||\bm{\ell}||`` tends to zero, while ``I`` increases such that the product ``I\bm{\ell}`` remains constant. Using the dyadic Green's functions from above, the fields of a Hertzian dipole can be calculated as
```math
\bm{E}_{\text{Hertz}}(\bm{r})
=
-\mathrm{j} k_0\, Z_\mathrm{F}\, I\,
\left[
\left(
\dfrac{3}{k^2\,||\bm{r}-\bm{r}'||^2}
+
\dfrac{3\mathrm{j}}{k\,||\bm{r}-\bm{r}'||} 
-1   
\right)
g_0(\bm{r},\bm{r}')
\,
\bm{e}_r \left(\bm{e}_r \cdot \bm{\ell}\right)
-
\left(
\dfrac{\mathrm{j}}{k\,||\bm{r}-\bm{r}'||}
+\dfrac{1}{k^2\,||\bm{r}-\bm{r}'||^2}
-1
\right)
g_0(\bm{r},\bm{r}')
\,
\bm{\ell}
\right]
```
and 
```math
\bm{H}_{\text{Hertz}}(\bm{r})
=I\,
\left(-\mathrm{j}k-\dfrac{1}{||\bm{r}-\bm{r}'||}\right) \bm{e}_r \times \bm{\ell}\, .
```

### Fitzgerald Dipole
A Fitzgerald dipole is a model for an elementary magnetic dipole. Formally, it can be regardes as a special magnetic current distribution of the form
```math
\bm{M}_{\text{Fitzgerald}}(\bm{r})
=
I_{\mathrm{m}}\bm{\ell}
\,
\delta({\bm{r}-\bm{r}'})\, ,
```
where ``I_{m}`` has the unit Volt and denotes the excitation of the dipole, ``\bm{\ell}`` denotes the oriented dipole length and ``\delta({\bm{r}-\bm{r}'})`` is the Dirac-Delta distribution. In the idealized model, ``||\bm{\ell}||`` tends to zero, while ``I_{\mathrm{m}}`` increases such that the product ``I_{\mathrm{m}}\bm{\ell}`` remains constant. Using the dyadic Green's functions from above, the fields of a Fitzgerald dipole can be calculated as
```math
\bm{E}_{\text{Fitzgerald}}(\bm{r})
=-I_{\mathrm{m}}\,
\left(-\mathrm{j}k-\dfrac{1}{||\bm{r}-\bm{r}'||}\right) \bm{e}_r \times \bm{\ell}
``` 
and
```math
\bm{H}_{\text{Fitzgerald}}(\bm{r})
=
-\mathrm{j} \dfrac{k_0}{Z_\mathrm{F}}\, I_{\mathrm{m}}\,
\left[
\left(
\dfrac{3}{k^2\,||\bm{r}-\bm{r}'||^2}
+
\dfrac{3\mathrm{j}}{k\,||\bm{r}-\bm{r}'||} 
-1   
\right)
g_0(\bm{r},\bm{r}')
\,
\bm{e}_r \left(\bm{e}_r \cdot \bm{\ell}\right)
-
\left(
\dfrac{\mathrm{j}}{k\,||\bm{r}-\bm{r}'||}
+\dfrac{1}{k^2\,||\bm{r}-\bm{r}'||^2}
-1
\right)
g_0(\bm{r},\bm{r}')
\,
\bm{\ell}
\right]\, .
``` 

## Far Fields of Dipoles
Often, we are interested in the radiated electric field of a dipole or a collection of dipoles at observation locations ``\bm{r}`` which are far away from the coordinate origin compared to the location ``\bm{r}'`` of the dipoles, i.e., ``r=||\bm{r}||\gg ||\bm{r}'||=r'``. In these cases, the term ``R=||\bm{r}-\bm{r}'||`` can be approximated by
```math
\lim \limits_{r\rightarrow \infty} ||\bm{r}-\bm{r}'|| = r -\bm{e}_r \cdot \bm{r}'\, , 
```
where ``\bm{e}_r=\bm{r}/r`` denotes the unit vector into the direction of the far-field observation location ``\bm{r}``. The dyadic Green's functions for the electric field take the form  
```math
\lim \limits_{r\rightarrow \infty}
\bm{\mathcal{G}}_{JE}(\bm{r},\bm{r}')
=
-\mathrm{j} k_0\, Z_{\mathrm{F}} \, \dfrac{\mathrm{e}^{-\mathrm{j}k_0 r}}{4\pi r}\, 
\left(
    \mathbf{I} -\bm{e}_r\bm{e}_r
\right)
\mathrm{e}^{\,\mathrm{j}k_0 \bm{e}_r\cdot \bm{r}'}
```
```math
\lim \limits_{r\rightarrow \infty}
\bm{\mathcal{G}}_{ME}(\bm{r},\bm{r}')
=
\mathrm{j} k_0 \, \dfrac{\mathrm{e}^{-\mathrm{j}k_0 r}}{4\pi r}\, 
\left(
    \bm{e}_r \times \mathbf{I}
\right)
\mathrm{e}^{\,\mathrm{j}k_0 \bm{e}_r\cdot \bm{r}'}
```

The radiated fields asymptotically far away from the sources have no components on radial direction and depend on the radial variable ``r`` only by the term ``\dfrac{\mathrm{e}^{-\mathrm{j}k_0 r}}{r}``. We have
```math
\lim \limits_{r\rightarrow \infty} 
\bm{E}(r,\vartheta, \varphi)
=
\bm{F}(\vartheta,\varphi)
\dfrac{\mathrm{e}^{-\mathrm{j}k_0 r}}{r}
```
where ``\bm{F}(\vartheta,\varphi)`` denotes the electric far-field pattern (unit: ``V``). 

Obviously, the far-field pattern of a Hertzian dipole can be calculated by
```math
\bm{F}_{\text{Hertz}}(\vartheta,\varphi)
=
I\, \dfrac{-\mathrm{j} k_0\, Z_{\mathrm{F}}}{4\pi}\,
\mathrm{e}^{\,\mathrm{j}k_0 \,\bm{e}_r(\vartheta,\varphi)\cdot \bm{r}'}\, 
\left(
    \mathbf{I} -\bm{e}_r(\vartheta,\varphi)\bm{e}_r(\vartheta,\varphi)
\right)\cdot \ell
\, ,
```
where ``I`` is the excitation of the dipole, ``\bm{\ell}`` is the oriented dipole length and ``\bm{r}'`` is the dipole location. The dyadic operator ``\left(
    \mathbf{I} -\bm{e}_r(\vartheta,\varphi)\bm{e}_r(\vartheta,\varphi)
\right)\cdot`` removes the field component in radial direction.

The far-field pattern of a Fitzgerald dipole can be calculated by
```math
\bm{F}_{\text{Fitzgerald}}(\vartheta,\varphi)
=
I_{\mathrm{m}}\, \dfrac{\mathrm{j} k_0}{4\pi}\,
\mathrm{e}^{\,\mathrm{j}k_0\, \bm{e}_r(\vartheta,\varphi)\cdot \bm{r}'}\, 
\bm{e}_r(\vartheta,\varphi)\times  \ell
\, ,
```
where ``I_{\mathrm{m}}`` is the excitation of the dipole, ``\bm{\ell}`` is the oriented dipole length and ``\bm{r}'`` is the dipole location. Due to the cross product with the radial unit vector, the far field has a zero component in radial direction.

The far fields of a Hertzian or a Fitzgerald dipole can be evaluated by calling the function
`farfield(sourcedipole::AbstractDipole,θ::Number,ϕ::Number,k₀::Flaot64) -> Tuple(Complex64,2)`

## Received signal
It can be shown that an antenna of which the radiated fields can be expressed via the (normalized) equivalent electric and magnetic current distributions ``\hat{\bm{J}}(\bm{r})`` and ``\hat{\bm{M}}(\bm{r})`` will receive the signal 
```math
b=
\dfrac{1}{2}
\iiint \limits_{V_{\text{src}}}
\left[
\bm{E}(\bm{r}) \cdot \hat{\bm{J}}(\bm{r})
-\bm{H}(\bm{r}) \cdot \hat{\bm{M}}(\bm{r})
\right]
\, \mathrm{d}v
``` 
when it is illuminated with an external electromagntic field ``\bm{E}(\bm{r})``, ``\bm{H}(\bm{r})``.

Consequently, a Hertzian dipole which has a (normalized) electric equivalent current distribution of ``\hat{\bm{J}}_{\text{Hertz}}(\bm{r}) = I\bm{\ell}\, \delta({\bm{r}-\bm{r}'})`` receives the signal ``b_{\text{Hertz}}=0.5\,I\, \bm{\ell} \cdot \bm{E}(\bm{r}')`` which is proportional to the parallel component of the incident electric field at the discrete location ``\bm{r}'``.

A Fitzgerald which has a (normalized) magnetic equivalent current distribution of ``\hat{\bm{M}}_{\text{Fitzgerald}}(\bm{r}) = I_{\mathrm{m}}\bm{\ell}\, \delta({\bm{r}-\bm{r}'})`` receives the signal ``b_{\text{Fitzgerald}}=-0.5\,I_{\mathrm{m}}\, \bm{\ell} \cdot \bm{H}(\bm{r}')`` which is proportional to the parallel component of the incident negative magnetic field at the discrete location ``\bm{r}'``.

## Dipole to Dipole Interactions
The interaction between two elementary dipoles is modeled by letting one of the dipoles generate an electromagnetic field, which is therafter received by the other dipole.

Since the recieved signal of a Hertzian dipole is directly proportional to the parallel component of the incident electric field at the dipole position, the received signals of Hertzian dipoles can be used to ideally probe the electric field. The interaction between two Hertzian dipoles can be used to represent the radiated electric field of one of the dipoles at the location of the other dipole.

Since the recieved signal of a Fitzgerald dipole is directly proportional to the parallel component of the negative incident magnetic field at the dipole position, the received signals of Fitzgerald  dipoles can be used to ideally probe the magnetic field. The interaction between two Fitzgerald dipoles can be used to represent the radiated magnetic field of one of the dipoles at the location of the other dipole.

The detailed derivation can be found in the "Theory" section.
The interaction between two dipoles is implemented in the function
(TODO: insert correct function)
Each elementary dipole is defined by its position ``\bm{r}``, its oriented length ``\bm{\ell}`` and its excitation magnitude ``I``.
Hertzian and Fitzgerald dipoles may be constructed using the constructors as in the following example.


### Interaction between two Hertzian Dipoles
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


### Interaction between two Fitzgerald Dipoles
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

### Interaction between Hertzian and Fitzgerald Dipole
Let ``I_{1}``, ``\bm{\ell}_1``, and ``\bm{r}_1`` denote the excitation, oriented (equivalent) length, and location of the transmitting Hertzian dipole and let ``I_{\mathrm{m},2}``, ``\bm{\ell}_2``, and ``\bm{r}_2`` be the excitation, oriented (equivalent) length, and location of the receiving Fitzgerald dipole. The received signal which represents the interaction between these two dipoles is given by
```math
b= \dfrac{1}{2}
I_1\, I_{\mathrm{m},2}
\left(-\mathrm{j}k-\dfrac{1}{||\bm{r}-\bm{r}'||}\right) \bm{\ell_2} \cdot \left( \bm{e}_r \times \bm{\ell_1}\right)
```
where ``\bm{e}_r`` denotes the unit vector pointing in the direction from one dipole location to the other.


The role of the transmitting and the receiving dipole may be exchanged without changing the result of the calculation (i.e., the electric field of a Fitzgerald dipole matches the electric field of a Fitzgerald dipole).