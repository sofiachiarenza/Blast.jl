# The algorithm
**Blast.jl** computes integrals of the form: 

```math
C_{ij}^{\mathrm{AB}}(\ell)= N(\ell) \int_0^{\infty} \mathrm{d} \chi_1 W_i^{\mathrm{A}}(\chi_1) \int_0^{\infty} \mathrm{d} \chi_2 W_j^{\mathrm{B}}(\chi_2) \int_0^{\infty} \mathrm{d} k\, k^2 P_{\mathrm{AB}}(k, \chi_1, \chi_2) \frac{j_\ell(k \chi_1) j_\ell(k \chi_2)}{(k \chi_1)^\alpha (k \chi_2)^\beta}
```
Currently, **Blast.jl** supports two kinds of kernels:

- Clustering: ``W_i^g(\chi)=H(z)n_i(z)b_g(z)``, ``\alpha = 0``

- Weak Lensing: ``W_i^s(\chi)= \frac{3H_0^2\Omega_m}{2a}\chi \int_z^{\infty} dz' n_i(z')\frac{\chi(z')-\chi}{\chi(z')}``, ``\beta = 2``

The algorithm takes advantage of the useful properties of Chebyshev polynomials, using them as a basis for efficiently decomposing the 3D matter power spectrum $P(k,\chi_1,\chi_2)$:

```math
P\left(k, \chi_1, \chi_2 ; \theta\right) \approx \sum_{n=0}^{n_{\max }} c_n\left(\chi_1, \chi_2 ; \theta\right) T_n(k)
```

This decomposition is advantageous for multiple reasons: first of all, the coefficients $c_n(\chi_1, \chi_2)$ can be quickly obtained through a DFT. Moreover, the approximation of a function on the basis of the Chebyshev polynomials is the most accurate compared to any other method for fixed number of interpolation points (see e.g. Trefethen (2019)). Finally, this decomposition enables a clear separation of geometric and cosmological components in the integrals, simplifying the overall calculation. We can write:

```math
w_{\ell}^{\mathrm{AB}}\left(\chi_1, \chi_2 ; \theta\right) \equiv \sum_{n=0}^{n_{\max }} c_n\left(\chi_1, \chi_2 ; \theta\right) \tilde{T}_{n ; \ell}^{\mathrm{AB}}\left(\chi_1, \chi_2\right),
```
where we defined: 
```math
\tilde{T}_{n ; \ell}^{\mathrm{AB}}\left(\chi_1, \chi_2\right) \equiv \int_{k_{\min }}^{k_{\max }} \mathrm{d} k f^{\mathrm{AB}}(k) T_n(k) j_{\ell}\left(k \chi_1\right) j_{\ell}\left(k\chi_2\right),
```
with:
```math
f^{\mathrm{AB}}(k)= \begin{cases}k^2 & \mathrm{AB}=g g, \\ 1 / k^2 & \mathrm{AB}=s s, \\ 1 & \mathrm{AB}=g s.\end{cases}
```
The dependence on the cosmological parameters is only present in the coefficients of the Chebyshev expansion of the power spectrum $c_n(\chi_1,\chi_2;\theta)$, while $\tilde{T}^{\mathrm{AB}}_{n;\ell}(\chi_1,\chi_2)$ is cosmology-independent as it is the integral of the two Bessel functions against the Chebyshev polynomials. This is the key idea of the algorithm: the $\tilde{T}^{\mathrm{AB}}_{n;\ell}(\chi_1,\chi_2)$ integrals are still challenging to compute for the presence of the Bessel functions, but they can be computed once-for-all.

The last ingredient for a successful computation of the integral is a change of variable: introducing $R \equiv \chi_2/\chi_1$, we can switch from the $\chi_1$-$\chi_2$ to the $\chi$-$R$ basis, which allows for a better sampling of the regions that most contribute to the integral, i.e., when $\chi_1 \approx \chi_2$ (or, equivalently,$R \approx 1$). In these new variables, the integral becomes:

```math
C_{ij}^{\mathrm{AB}}(\ell) = \int_0^{\infty} \mathrm{d} \chi \int_0^1 \mathrm{d} R \, \chi \left[ \mathcal{K}_i^{\mathrm{A}}(\chi) \mathcal{K}_j^{\mathrm{B}}(R \chi) + \mathcal{K}_j^{\mathrm{B}}(\chi) \mathcal{K}_i^{\mathrm{A}}(R \chi) \right] w_{\ell}^{\mathrm{AB}}(\chi, R \chi)

```
with:
```math
\mathcal{K}_i^{\mathrm{A}}(\chi)= \begin{cases}K_i^{\mathrm{A}}(\chi) & \text { for clustering, } \\ K_i^{\mathrm{A}}(\chi) / \chi^2 & \text { for lensing. }\end{cases}
```




