# Electromagnetic Fields of Elementary Dipoles
The radiated fields of an antenna can be approximated by the fields of elementary dipoles, i.e., the fields of very short electrical dipoles (so-called Hertzian dipoles) and very short magnetic dipoles (so-called Fitzgerald dipoles).

Since the equivalent current densities for these very short elementary dipoles are particularly simple, it is possible to find analytical expressions for their respective fields using Green's functions.

## Radiated Fields
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
### Radiated Fields of a Hertzian Dipole
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

### Radiated Fields of a Fitzgerald Dipole
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

## Radiated Far Fields of Elementary Dipoles
Often, we are interested in the radiated electric field of an elementary dipole or a collection of elementary dipoles at observation locations ``\bm{r}`` which are far away from the coordinate origin compared to the location ``\bm{r}'`` of the dipoles, i.e., ``r=||\bm{r}||\gg ||\bm{r}'||=r'``. In these cases, the term ``R=||\bm{r}-\bm{r}'||`` can be approximated by
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

## [Received Signal](@id receivedsignaldipoles)
Let ``\hat{\bm{J}}(\bm{r})`` and ``\hat{\bm{M}}(\bm{r})`` be the normalized equivalent electric and magnetic volume currents for the fields radiated by an antenna such that ``\bm{J}(\bm{r}) = a\,\hat{\bm{J}}(\bm{r})`` and ``\bm{M}(\bm{r}) = a\,\hat{\bm{J}}(\bm{r})`` are the actual equivalent volume current densities when the antenna's transmit port is excited by a signal with wave amplitude ``a \in \mathbb{C} \, \sqrt{\mathrm{W}}``.
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