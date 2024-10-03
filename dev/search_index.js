var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API","title":"API reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Blast.compute_T̃\nBlast.bessel_cheb_eval\nBlast.get_clencurt_weights\nBlast.get_clencurt_grid\nBlast.w_ell_tullio","category":"page"},{"location":"api/#Blast.compute_T̃","page":"API","title":"Blast.compute_T̃","text":"compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb = 119, N=2^(15)+1)\n\nCompute integrals of the Bessels function and the Chebyshev polynomials. This is the precomputation part of the code. The parameters are:\n\n- ℓ: Multipole order\n\n- χ: Array containing values of the comoving distance. \n\n- R: Array containing values for the R=χ₁/χ₂ variable.\n\n- kmin-kmax: Integration range in k.\n\n- β: Exponent of the k dependence of the integral. This parameter depends on the combination of tracers: β=2,-2,0 for clustering, cosmic shear and the cross-correlation respectively.\n\n- n_cheb: Number of chebyshev polynomials used in the approximation of the power spectra.\n\n- N: Number of integration points in k.\n\n\n\n\n\n","category":"function"},{"location":"api/#Blast.bessel_cheb_eval","page":"API","title":"Blast.bessel_cheb_eval","text":"bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)\n\nReturn the Chebyshev polynomials up to order 'n_cheb+1' and the Bessel function of order 'ℓ' evaluated on the grid of 'N' Chebyshev points in the interval ['kmin', 'kmax'] and on the specified 'χ' points. \n\n\n\n\n\n","category":"function"},{"location":"api/#Blast.get_clencurt_weights","page":"API","title":"Blast.get_clencurt_weights","text":"get_clencurt_weights(kmin::Number, kmax::Number, N::Number)\n\nReturn the set of 'N' weights needed to perform the integration with the Clenshaw-Curtis quadrature rule. The weights are rescaled between 'kmin' and 'kmax'.  \n\n\n\n\n\n","category":"function"},{"location":"api/#Blast.get_clencurt_grid","page":"API","title":"Blast.get_clencurt_grid","text":"get_clencurt_grid(kmin::Number, kmax::Number, N::Number)\n\nReturn the integration points in k. They are a set of 'N' Chebyshev points rescaled between 'kmin' and 'kmax'.\n\n\n\n\n\n","category":"function"},{"location":"api/#Blast.w_ell_tullio","page":"API","title":"Blast.w_ell_tullio","text":"w_ell_tullio(c::AbstractArray,T::AbstractArray)\n\nCompute the tensor contraction of the chebyshev coefficients of the power spectrum 'c' and the precomputed integrals 'T' to obtain the projected matter densities.\n\n\n\n\n\n","category":"function"},{"location":"#Blast.jl","page":"Home","title":"Blast.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Hello, welcome to Blast.jl documentation! ","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Sofia Chiarenza, PhD student at Waterloo Centre for Astrophysics.","category":"page"}]
}
