# Fast Algorithms in Spherical Measurement Setups 

Very efficient post-processing algorithms exist for measurements of the ``S_{12}``-parameter between a radiating probe antenna and a rotating receiving antenna under test (AUT). The AUT is rotated around the coordinate origin by the angles ``\vartheta`` and ``\varphi``. The probe antenna remains at a fixed position in the positive ``z``-diretion but may be rotated around the ``z``-axis by the angle ``\chi`` to  realize different polarizations of the incident field at the AUT. The measurement setup is depicted in the figure below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/geom_measure_dark.svg" width="200">
  <source media="(prefers-color-scheme: light)" srcset="../assets/geom_measure.svg" width="200" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Setup for the antenna field measurement.
  </figcaption>
</figure>
<br/>
```

## Fast Evaluation of ``S_{12}``-Parameters in Spherical Measurement Setups 

The  ``S_{12}``-parameter can be calculated by (cf. [1, eq. (4.40)])


```math
S_{12}(\vartheta, \varphi, \chi)=\sum \limits_{s=1}^2 \sum \limits_{\ell=1}^L \sum \limits_{m=-\ell}^\ell \sum \limits_{\mu=-\ell}^\ell \beta_{s\ell m}^{\mathrm{aut}}\,\mathrm{e}^{-\mathrm{j} m \varphi}\, \mathrm{d}_{m,\mu}^{\ell}(\vartheta)\, \mathrm{e}^{-\mathrm{j} \mu \chi}\,\alpha_{s\ell \mu}^{(1),\mathrm{pro}}\,.
```

Here. ``\alpha_{s\ell \mu}^{(1),\mathrm{pro}}`` are the spherical mode coefficients of the incident probe antenna field (normalized to a unit excitation of the probe), ``\beta_{s\ell m}^{\mathrm{aut}}`` are the spherical receive coefficients of the AUT and ``\mathrm{d}_{m,\mu}^{\ell}(\vartheta)`` is the so-called Wigner-d-matrix which has a finite Fourier-series expansion

```math
\mathrm{d}_{m,\mu}^{\ell}(\vartheta)= \mathrm{j}^{m-\mu} \sum \limits_{m_\vartheta=-\ell}^\ell \, \Delta_{m_\vartheta,\mu}^\ell\, \Delta_{m_\vartheta,m}^\ell\, \mathrm{e}^{-\mathrm{j} m_\vartheta \, \vartheta}\, ,
```

where the deltas

```math
\Delta_{m_\vartheta,m}^\ell= \mathrm{d}_{m,\mu}^{\ell}(\mathrm{\pi}/2)
```
correspond to the calue of the Wigner-d-matrix at ``\vartheta=\mathrm{\pi/2}``.

Substituting the series expansion into the original equation, one obtains

```math
S_{12}(\vartheta, \varphi, \chi)=
\sum \limits_{\mu=-L}^L 
\mathrm{e}^{-\mathrm{j} \mu \chi}
\underbrace{
\sum \limits_{m=-L}^L
\mathrm{e}^{-\mathrm{j} m \varphi}\,
\underbrace{
\sum \limits_{m_\vartheta=-L}^L
\mathrm{e}^{-\mathrm{j} m_\vartheta \, \vartheta}\, 
\underbrace{
\mathrm{j}^{m-\mu}
\sum \limits_{\ell=\max(|m|,|m_\vartheta|,1)}^L
\Delta_{m_\vartheta,\mu}^\ell\, \Delta_{m_\vartheta,m}^\ell\,
 \underbrace{
 \sum \limits_{s=1}^2 
 \alpha_{s\ell \mu}^{(1),\mathrm{pro}}\,
\beta_{s\ell m}^{\mathrm{aut}}
}_{u_{\mu,m,\ell}}
}_{\tilde{v}_{\mu,m,m_\vartheta}}
}_{v_{\mu,m}(\vartheta)}
}_{w_\mu(\vartheta,\varphi)}\,.
```

The form of the above equation gives rise to an efficient algorithm to evaluate the ``S_{12}``-parameter on a spherical measurement surface by successively calculating 
`` \alpha_{s\ell \mu}^{(1),\mathrm{pro}}\,
\beta_{s\ell m}^{\mathrm{aut}} \rightarrow u_{\mu,m,\ell} \rightarrow \tilde{v}_{\mu,m,m_\vartheta} \rightarrow v_{\mu,m}(\vartheta) \rightarrow w_\mu(\vartheta,\varphi)``. The occurring Fourier series are evaluated using a fast Fourier transform giving rise to the resulting data available at a regular spherical sampling grid.

### Radiating Hertzian Dipole Probe and Receiving AUT
If we want to calculate the transmission between a dipole probe and the AUT, we need to know the spherical wave expansion of the probe field at the AUT location. 
The incident field at the coordinate origin of a Hertzian dipole, oriented in ``x``- direction with a given dipole moment ``I\ell`` at a given location ``\bm{r}= r\, \bm{e}_z`` on the ``z``-axis, is characterized by the spherical wave coefficients 
```math
\alpha_{1,\ell,1}^{(1)}=\dfrac{-\mathrm{j}}{4} \,k\, \sqrt{Z_{\mathrm{F}}}\, I\ell
\, \sqrt{\dfrac{2\ell+1}{\pi}}\, h_\ell^{(2)}(kr)
```
```math
\alpha_{1,\ell,-1}^{(1)}=\dfrac{-\mathrm{j}}{4} \,k\, \sqrt{Z_{\mathrm{F}}}\, I\ell
\, \sqrt{\dfrac{2\ell+1}{\pi}}\, h_\ell^{(2)}(kr)
```
```math
\alpha_{2,\ell,1}^{(1)}=\dfrac{1}{4} \,k\, \sqrt{Z_{\mathrm{F}}}\, I\ell
\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{1}{kr} \, \dfrac{\mathrm{d}}{\mathrm{d} kr} \left\{kr\,h_\ell^{(2)}(kr)\right\}
```
```math
\alpha_{2,\ell,-1}^{(1)}=\dfrac{-1}{4} \,k\, \sqrt{Z_{\mathrm{F}}}\, I\ell
\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{1}{kr} \, \dfrac{\mathrm{d}}{\mathrm{d} kr} \left\{kr\,h_\ell^{(2)}(kr)\right\}\,,
```

where ``h_\ell^{(2)}(kr)`` are the spherical Hankel functions of second kind.
All incident field coefficients with ``|m|\neq 1`` are equal to zero. 

### Regularly Sampled Far Field of AUT
The far field pattern of the AUT can be calculated by determining the received signal of a plane wave illumination traveling into the desired propagation directions. To determine the correctly normalized far-field pattern (with unit V), we utilize the spherical expansion of a ``x``-polarized plane wave traveling into negativ ``z``-direction with the coefficients

```math
\alpha_{1,\ell,1}^{(1)}=- \mathrm{j}^{\ell}\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{\sqrt{Z_{\mathrm{F}}}}{2}
```

```math
\alpha_{2,\ell,1}^{(1)}=- \mathrm{j}^{\ell}\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{\sqrt{Z_{\mathrm{F}}}}{2}
```

```math
\alpha_{1,\ell,-1}^{(1)}=- \mathrm{j}^{\ell}\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{\sqrt{Z_{\mathrm{F}}}}{2}
```

```math
\alpha_{2,\ell,-1}^{(1)}= \mathrm{j}^{\ell}\, \sqrt{\dfrac{2\ell+1}{\pi}}\, \dfrac{\sqrt{Z_{\mathrm{F}}}}{2}\, .
```
All incident field coefficients with ``|m|\neq 1`` are equal to zero. 
Using the reciprocity relation ``\beta_{s\ell m}^{\mathrm{aut}}=(-1)^{m}\, \alpha_{s,\ell,- m}^{(4),\mathrm{aut}}``, we can efficiently calculate the desired far fields for regularly sampled ``\vartheta`` and ``\varphi``.

# Fast Reconstruction of Expansion Coefficients from ``S_{12}``-Parameters in Spherical Measurement Setups

The  ``S_{12}``-parameter can be calculated by (cf. [1, eq. (4.40)])

```math
S_{12}(\vartheta, \varphi, \chi)=
\sum \limits_{\mu=\pm 1} 
\mathrm{e}^{-\mathrm{j} \mu \chi}
\underbrace{
\sum \limits_{m=-L}^L
\mathrm{e}^{-\mathrm{j} m \varphi}\,
\underbrace{ 
\sum \limits_{\ell=\max(|m|,|m_\vartheta|,1)}^L
\mathrm{d}_{\mu,m}^\ell(\vartheta)\,
 \underbrace{
 \sum \limits_{s=1}^2 
 \alpha_{s\ell \mu}^{(1),\mathrm{pro}}\,
\beta_{s\ell m}^{\mathrm{aut}}
}_{u_{\mu,m,\ell}}
}_{v_{\mu,m}(\vartheta)}
}_{w_\mu(\vartheta,\varphi)}\,.
```

Here, ``\mathrm{d}_{m,\mu}^{\ell}(\vartheta)`` is the so-called Wigner-d-matrix, ``\alpha_{s\ell \mu}^{(1),\mathrm{pro}}`` are the spherical mode coefficients of the incident probe antenna field (normalized to a unit excitation of the probe), and ``\beta_{s\ell m}^{\mathrm{aut}}`` are the spherical receive coefficients of the AUT which we aim to reconstruct from the measured data. Note, that a first-order probe is assumed here (i.e., all incident modes with ``\mu\neq \pm 1`` are zero), reducing the number of summation terms in the ``\mu``-summation to two.  

Starting at ``S_{12} (\vartheta, \varphi, \chi)``, we can successively calculate ``S_{12}(\vartheta, \varphi, \chi) \rightarrow w_\mu(\vartheta,\varphi) \rightarrow v_{\mu,m}(\vartheta) \rightarrow u_{\mu,m,\ell}``
using the orthogonality relations
```math
\int \limits_{0}^{2\pi}
\mathrm{e}^{\, \mathrm{j} (m-m') \varphi} \, \mathrm{d} \varphi = 2\pi \, \delta_{m,m'}
```
and
```math
\int \limits_{0}^{\pi}
\mathrm{d}_{m,m'}^\ell(\vartheta)
\,
\mathrm{d}_{m,m'}^{\ell'}(\vartheta)
\,
\sin(\vartheta) 
\, \mathrm{d} \vartheta = \dfrac{2}{2\ell+1} \, \delta_{\ell,\ell'}\,
```
from where the ``\beta_{s\ell m}^{\mathrm{aut}}`` can be found by solving a ``2\times 2`` system of linear equations
```math
	\begin{bmatrix}
		{u_{+1,m,\ell}}\\[0.5em]
		{u_{-1,m,\ell}}
	\end{bmatrix}
=
\begin{bmatrix}
\alpha_{1,\ell ,+1}^{(1),\,\mathrm{pro}} 
&
\alpha_{2,\ell ,+1}^{(1),\,\mathrm{pro}} 
\\[0.5em]
\alpha_{1,\ell ,-1}^{(1),\,\mathrm{pro}} 
&
\alpha_{2,\ell ,-1}^{(1),\,\mathrm{pro}}	
\end{bmatrix}
\,
	\begin{bmatrix}
	{\beta}_{1, \ell, m }^{\mathrm{aut}}\\[0.5em]
	{\beta}_{2, \ell, m }^{\mathrm{aut}}
\end{bmatrix}\, .
```

The individual transformation steps are as follows. Since we only have to consider ``\mu``-values for ``\mu = \pm 1``, we can calculate the ``\rightarrow w_\mu(\vartheta,\varphi)`` from the ``S_{12} (\vartheta, \varphi, \chi)``-values measured at ``\chi \in \{0, \pi/2\}``
by
```math
w_{+1}(\vartheta,\varphi)= 
\dfrac{1}{2}
\left[
S_{12} (\vartheta, \varphi, \chi=0) + \mathrm{j}\, S_{12} (\vartheta, \varphi, \chi=\pi/2)
\right]
```
and
```math
w_{-1}(\vartheta,\varphi)= 
\dfrac{1}{2}
\left[
S_{12} (\vartheta, \varphi, \chi=0) - \mathrm{j}\, S_{12} (\vartheta, \varphi, \chi=\pi/2)
\right]
``` 

---
## [References](@id refs)
- [1] J. E. Hansen, ed., *Spherical Near-Field Antenna Measurements*, The Institution of Engineering and Technology, Michael Faraday House, Six Hills Way, Stevenage SG1 2AY, UK: IET, 1988.