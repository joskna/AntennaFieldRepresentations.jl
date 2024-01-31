# Spherical Vector-Wave Expansion of General Time Harmonic Electromangentic Fields
Away from sources, the time harmonic electric and magnetic fields must fulfill the homogeneous curl-curl-equation 
```math
\nabla \times \nabla \times \bm{E} - k_0^2 \bm{E} =\bm{0}\, .
```

Arbitrary source-free solutions of Maxwell's equations can be expanded into a series of spherical vector-wave functions as (see pp. 9ff of [^1])
```math
\bm{E}(r, \vartheta, \varphi) =
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(c)}\,
\bm{F_{s\ell m}^{(c)}}(r,\vartheta, \varphi)
```
```math
\bm{H}(r, \vartheta, \varphi) =\mathrm{j}\,
\dfrac{k_0}{\sqrt{Z_{\mathrm{F}}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(c)}\,
\bm{F_{3-s,\ell m}^{(c)}}(r,\vartheta, \varphi)
```  
where ``Z_{\mathrm{F}}\sqrt{\mu_0/ \varepsilon_{0}}\approx 376.730\,313\,669\, \mathrm{\Omega}`` is the impedance of free space and ``\bm{F_{s\ell m}^{(c)}}(r,\vartheta, \varphi)`` are the spherical vector wave mode functions, which are solutions to the homogeneous curl-curl equation (except in the origin of the coordinate system, where some types of vector mode functions become singular). Due to only two of the radial dependencies being linearly independent, it is sufficient to let the sum over ``c`` to be an arbitrary pair of two iindices from the set ``c \in \{1,2,3,4\}``. Usually one will choose the pair ``c=1`` and ``c=4``, where only the ``c=4``-type modes are needed to describe purely radiated fields and only the ``c=1``-type modes are needed to describe purely incident fields.

In the numerical implementation, the coefficients ``\alpha_{s \ell m}`` are stored in a vector, where the triple index
``s \ell m`` is mapped to a single consecutively running index ``j`` by the function `sℓm_to_j(s,ℓ,m)`. The inverse map from a single index to the triple index ``s \ell m`` is implemented in the function `j_to_sℓm(j)`.

The vector of coefficients completely defines the fields of a corresponding radiating or incident spherical wave expansion. 



## Definition of Spherical Vector-Mode Functions

Explicit expressions for the spherical vector-mode functions are (see p. 13 of [^2])
```math
\bm{F}_{1\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}
\left(
\dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \, \dfrac{\mathrm{j}m \,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\sin \vartheta} \bm{e}_\vartheta
-\dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \, \dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta} \bm{e}_\varphi
\right)
```
```math
\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}
\left(
\dfrac{\sqrt{\ell(\ell+1)}}{k_0r}\, z_\ell^{(c)}(k_0r)\, \overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) \bm{e}_r
\right.
\phantom{\dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta}~~~~~}
\\
\phantom{\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}}
+ \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \right\}\, \dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta} \bm{e}_\vartheta
\phantom{~~~~~}
\\ \left.
\phantom{\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}}
+ \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \right\}\, \dfrac{\mathrm{j}m \,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\sin \vartheta} \bm{e}_\varphi\
\right)\, .
```


To ensure that the vector-mode functions ``\bm{F_{s\ell m}^{(c)}}(r,\vartheta, \varphi)`` are solutions to the homogeneous curl-curl equation, they are constructed from solutions of the homogeneous scalar Helmholtz equation as (see, e.g., page 393 of [^3])
```math
\bm{M}_{\ell m}^{(c)}= \nabla f_{\ell m}^{(c)} \times \bm{e}_r \overset{!}{=} \bm{F}_{1\,\ell m}^{(c)}
```
and
```math
\bm{N}_{\ell m}^{(c)}= \dfrac{1}{k_0} \nabla \times \bm{M}_{\ell m}^{(c)} \overset{!}{=} \bm{F}_{2\,\ell m}^{(c)}
``` 
where ``\bm{e}_r`` is the unit vector in radial direction.


General solution to the scalar Helmholtz equation 
```math
\Delta f +k_0^2 f=0\, .
```
can be constructed as superpositions of functions from the set
```math
f_{s\ell m}^{(c)}(r, \vartheta, \varphi)
=
\begin{cases}
\dfrac{\mathrm{j}_{\ell}(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=1\\[1.5em]
\dfrac{\mathrm{\mathrm{n}}_\ell(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=2\\[1.5em]
\dfrac{\mathrm{\mathrm{h}}_\ell^{(1)}(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=3\\[1.5em]
\dfrac{\mathrm{\mathrm{h}}_\ell^{(2)}(k_0 r)}{\sqrt{\ell(\ell+1)}}\,
 \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=4
\end{cases} 
```
where 
- ``\mathrm{j}_{\ell}(k_0 r)`` are spherical Bessel functions which represent the radial behavior standing waves (or *incident* waves)
- ``\mathrm{n}_{\ell}(k_0 r)`` are spherical Neumann functions which usually do not represent a physically meaningful solution
- ``\mathrm{h}_\ell^{(1)}(k_0 r) = \mathrm{j}_{\ell}(k_0 r) +\mathrm{j}\, \mathrm{n}_{\ell}(k_0 r)`` are spherical Hankel functions of first kind which represent the radial behavior of waves coming from infinity being absorbed at the coordinate origin
- ``\mathrm{h}_\ell^{(2)}(k_0 r) = \mathrm{j}_{\ell}(k_0 r) -\mathrm{j}\, \mathrm{n}_{\ell}(k_0 r)`` are spherical Hankel functions of second kind which represent the radial behavior of waves traveling from the coordinate origin toward infinity
- ``\mathrm{Y}_{\ell m}(\vartheta, \varphi)`` are the spherical Harmonics
and ``k_0=2\pi/\lambda`` is the wavenumber defined by the frequency of the electromagnetic wave.



The spherical Harmonics are defined as (see §14.30 of [^1])
```math
\mathrm{Y}_{\ell m}(\vartheta, \varphi) 
=
\begin{cases}
\dfrac{1}{\sqrt{2\pi}}\, \overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) \mathrm{e}^{\, \mathrm{j}m \varphi} & \text{for }|m|\leq \ell\\
0 & \text{esle}
\end{cases}
``` 
where the normalized associated Legendre functions ``\overline{\mathrm{P}}_{\ell}^{m}(x)`` are defined in ordinary Legendre polynomials ``\mathrm{P}_{\ell}(x)`` 
as (compare §14.6 and §14.7 of [^1])
```math
\overline{\mathrm{P}}_{\ell}^{m}(x)
=
\begin{cases}
\sqrt{\dfrac{(2\ell+1)(\ell-m!)}{2(\ell+m)!}}
(-1)^m
\,
(1-x^2)^{m/2}
\, \dfrac{\mathrm{d}^m \mathrm{P}_{\ell}(x)}{\mathrm{d}^m} & m\geq 0\\[2em]
(-1)^m  \overline{\mathrm{P}}_{\ell}^{-m}(x) & m<0\, .
\end{cases}
```

!!! note
    Different sign conventions are used by different authors when the normalized associated Legendre functions are defined. Conventions may differ by a sign-factor of ``(-1)^m`` (this factor is often called *Condon-Shortley phase factor*). We follow the sign convention of [^1], which includes the Condon-Shortley phase factor of ``(-1)^m`` in the definition of the normalized associated Legendre functions for positive ``m``. The main reference [^2] for our implementation does not include the phase factor in the definition for the associated Legendre functions, but includes an explicit ``(-1)^m``-term to the definition of the mode functions (which we don't do), such that in the end, our definition of the final mode functions coincides with the definition used in [^2].
---
!!! warning
    A straightforward implementation of the above definition for the normalized associated Legendre functions will fail already for moderately large mode orders ``\ell`` due to the involved factorials. A more sensible definition is based or recurrence relations[^4]
---

## Far-Field Expressions
Due to the asymptotic behavior of the radial dependencies 
```math
\lim \limits_{k_0r \rightarrow \infty} \mathrm{h}_\ell^{(2)}(k_0r) = \mathrm{j}^{\ell+1} \dfrac{\mathrm{e}^{-\mathrm{j} k_0r}}{k_0r} 
```
```math
\lim \limits_{k_0r \rightarrow \infty}  \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \mathrm{h}_\ell^{(2)}(k_0r) \right\} = \mathrm{j}^{\ell} \dfrac{\mathrm{e}^{-\mathrm{j} k_0r}}{k_0r} 
```
for distances approaching infinity, it is convenient to define far-field expressions of the mode functions as
```math
\bm{K}_{1\,\ell m}^{(4)}(\vartheta, \varphi) =
\lim \limits_{k_0r \rightarrow \infty}
\dfrac{k_0r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{F}_{1\,\ell m}^{(4)}(r,\vartheta, \varphi)
```
```math
\bm{K}_{2\,\ell m}^{(4)}(\vartheta, \varphi) =
\lim \limits_{k_0r \rightarrow \infty}
\dfrac{k_0r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{F}_{2\,\ell m}^{(4)}(r,\vartheta, \varphi)\, .
```

Explicit expressions for the far-field mode functions are given as
```math
\bm{K}_{1\,\ell m}^{(4)}(\vartheta, \varphi) =
\dfrac{\mathrm{j}^{\ell+1}}{\sqrt{\ell(\ell+1)}}
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2\pi}}
\left(
\dfrac{\mathrm{j}m\, \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\sin \vartheta} \bm{e}_\vartheta
-
\dfrac{\mathrm{d} \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\mathrm{d} \vartheta} \bm{e}_\varphi  
\right)
```
```math
\bm{K}_{2\,\ell m}^{(4)}(\vartheta, \varphi) =
\dfrac{\mathrm{j}^{\ell}}{\sqrt{\ell(\ell+1)}}
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2\pi}}
\left(
\dfrac{\mathrm{d} \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\mathrm{d} \vartheta} \bm{e}_\vartheta
+
\dfrac{\mathrm{j}m\, \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\sin \vartheta}
 \bm{e}_\varphi  
\right)\, .
```
!!! note
    This definition of the far-field mode functions lacks the arbitrary factor of ``\sqrt{4\pi}`` which is introduced in [^2] .
---


Naturally, the far field ``\bm{F}(\vartheta, \varphi) = \lim \limits_{k_0r \rightarrow \infty} \dfrac{r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{E}(r,\vartheta, \varphi) `` can be calculated from a given set of expansion coefficients ``\alpha_{s\ell m}^{(4)}`` as
```math
\bm{F}(\vartheta, \varphi) =
\sqrt{Z_{\mathrm{F}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(4)}\,
\bm{K_{s\ell m}^{(c)}}(\vartheta, \varphi)\, .
```
---
## [References](@id refs)
[^1]:  [F. W. J. Olver and National Institute of Standards and Technology (U.S.), eds., "*NIST Handbook of Mathematical Functions*", Cambridge ; New York: Cambridge University Press : NIST, 2010](https://dlmf.nist.gov/).
[^2]:  J. E. Hansen, ed., *Spherical Near-Field Antenna Measurements*, The Institution of Engineering and Technology, Michael Faraday House, Six Hills Way, Stevenage SG1 2AY, UK: IET, 1988.
[^3]:  J. A. Stratton, "*Electromagnetic Theory*", McGraw-Hill, 1st ed. International series in pure and applied physics, OCLC: 536704, New York, 1941.
[^4]: [T. Limpanuparb, J. Milthorpe: "Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications", arXiv 2014](https://arxiv.org/abs/1410.1748).
