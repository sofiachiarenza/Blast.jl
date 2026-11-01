# Fixed non-Limber tensor banks.
#
# A bank is read-only data: the eight `T̃` tensors (one per geometry) together
# with the spectral sites they were generated for and a precomputed
# compatibility signature. Banks are never mutated after construction, so they
# may be shared by every setup, workspace and chain.

"""
    BankSignature

Immutable, precomputed identity of a tensor bank and the spectral contract it
implies. `digest` is a versioned SHA-256 over a canonical byte encoding of: the
format version, basis, coordinate and (direct) representation conventions, the
dataset identity string (the artifact git-tree id for published data), and the
exact bits and order of the ℓ, χ, R axes and the k sites. It is deterministic
across processes and Julia versions.
It is attached to plans, power spectra and projected-matter components, and
compared at public boundaries so that mismatched configurations fail before
any contraction. It is computed once at bank construction, never in hot loops.
"""
struct BankSignature
    basis::Symbol
    n_k::Int
    digest::NTuple{4, UInt64}
    geometry::NTuple{4,UInt64}
end
BankSignature(basis,n_k,digest)=BankSignature(basis,n_k,digest,PRODUCTION_GEOMETRY_SIGNATURE)

const BankAxes=NamedTuple{(:ell,:chi,:R),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
_production_bank_axes()=(ell=copy(ℓ_nonlimber),chi=copy(χ),R=copy(R))
function _copy_bank_axes(axes)
    a=(ell=collect(Float64,axes.ell),chi=collect(Float64,axes.chi),R=collect(Float64,axes.R))
    !isempty(a.ell) && all(isfinite,a.ell) && all(>(0),a.ell) ||
        throw(ArgumentError("bank multipoles must be finite and positive"))
    all(v->!isempty(v) && all(isfinite,v) && first(v)>0 && all(>(0),diff(v)),(a.chi,a.R)) ||
        throw(ArgumentError("bank chi/R axes must be finite, positive and strictly increasing"))
    return a
end

"""
    TensorBank

Eight `T̃` tensors of shape `(ℓ, χ, R, k)`, the log₁₀ k `sites` they were
generated at, and their [`BankSignature`](@ref). Field names follow the
geometry labels used by `_PROJECTED_MATTER_COMPONENTS`.
"""
struct TensorBank{A<:AbstractArray{Float64,4}}
    T_2_00::A
    T_0_00::A
    T_minus2_00::A
    T_0_02::A
    T_0_20::A
    T_2_02::A
    T_2_20::A
    T_2_22::A
    sites::Vector{Float64}
    signature::BankSignature
    axes::BankAxes
end
TensorBank(a::A,b::A,c::A,d::A,e::A,f::A,g::A,h::A,sites,signature) where {A<:AbstractArray{Float64,4}} =
    TensorBank(a,b,c,d,e,f,g,h,sites,signature,_production_bank_axes())
bank_axes(bank::TensorBank)=bank.axes

const _BANK_GEOMETRIES = (:T_2_00, :T_0_00, :T_minus2_00, :T_0_02, :T_0_20, :T_2_02, :T_2_20, :T_2_22)

const _SIGNATURE_FORMAT = "blast-tensor-bank-signature-v1"

function _digest_text!(ctx, x::AbstractString)
    bytes = codeunits(x)
    SHA.update!(ctx, reinterpret(UInt8, [htol(UInt64(length(bytes)))]))
    SHA.update!(ctx, bytes)
end
function _digest_floats!(ctx, v::AbstractVector{Float64})
    SHA.update!(ctx, reinterpret(UInt8, [htol(UInt64(length(v)))]))
    SHA.update!(ctx, reinterpret(UInt8, htol.(collect(v))))
end

"""
    bank_signature_digest(basis, dataset_id, sites)

Canonical signature for a direct-representation bank on the production ℓ, χ, R
axes. `dataset_id` identifies the tensor data: `"artifact:<name>:<git-tree-sha1>"`
for a published artifact, `"private:<label>"` for an artificially constructed
bank (so the two can never collide). Computed once at construction, never in
differentiated code.
"""
function bank_signature_digest(basis::Symbol, dataset_id::AbstractString, sites::AbstractVector{<:Real}; axes=_production_bank_axes())
    ctx = SHA.SHA256_CTX()
    for text in (_SIGNATURE_FORMAT, String(basis), "coordinate=log10(k/Mpc^-1)",
                 "representation=direct;tensor_axes=(ell,chi,R,k)", dataset_id)
        _digest_text!(ctx, text)
    end
    for axis in (axes.ell, axes.chi, axes.R, sites)
        _digest_floats!(ctx, Float64.(axis))
    end
    d = SHA.digest!(ctx)
    words = ltoh.(reinterpret(UInt64, d))
    return BankSignature(basis, length(sites), (words[1], words[2], words[3], words[4]),
        _geometry_signature(axes.ell,axes.chi,axes.R))
end

"""
    TensorBank(tensors, sites; basis, dataset_id)

Validate shapes and build the signature. Use a `dataset_id` starting with
`"private:"` for artificial/test banks; published data is built only by the
loader below, with the artifact's tree id.
"""
# Shared cold construction: geometry set, shapes, signature.
function _build_tensor_bank(tensors, sites, basis::Symbol, dataset_id::AbstractString; axes=_production_bank_axes())
    axes=_copy_bank_axes(axes)
    Set(keys(tensors)) == Set(_BANK_GEOMETRIES) ||
        throw(ArgumentError("a tensor bank needs exactly the geometries $(_BANK_GEOMETRIES)"))
    expected = (length(axes.ell), length(axes.chi), length(axes.R), length(sites))
    for name in _BANK_GEOMETRIES
        size(tensors[name]) == expected || throw(DimensionMismatch(
            "tensor $name has size $(size(tensors[name])); expected (ℓ, χ, R, k) = $expected"))
    end
    signature = bank_signature_digest(basis, dataset_id, sites; axes)
    return TensorBank(tensors.T_2_00, tensors.T_0_00, tensors.T_minus2_00, tensors.T_0_02,
                      tensors.T_0_20, tensors.T_2_02, tensors.T_2_20, tensors.T_2_22,
                      collect(Float64, sites), signature, axes)
end

"""
    TensorBank(tensors, sites; basis, dataset_id)

Artificial/test bank. `dataset_id` must start with `"private:"`; the
`"artifact:"` namespace is reserved for the trusted loader of published data
and is rejected here.
"""
function TensorBank(tensors::NamedTuple{<:Any, <:NTuple{8, AbstractArray{Float64,4}}},
                    sites::AbstractVector{<:Real}; basis::Symbol, dataset_id::AbstractString, axes=_production_bank_axes())
    startswith(dataset_id, "private:") || throw(ArgumentError(
        "artificial banks need a dataset_id starting with \"private:\"; \"artifact:\" is reserved for published data"))
    return _build_tensor_bank(tensors, sites, basis, dataset_id; axes)
end

# The legacy Chebyshev bank (artifact `T_tildes`, 128 χ × 64 R × 161 k). Its
# signature is a compile-time constant, so tagging legacy components never
# touches the filesystem.
const _LEGACY_BANK_NAME = "T_tildes_128_64_161"
# Authoritative dataset identity: the git tree id pinned in Artifacts.toml.
const _LEGACY_DATASET_ID = "artifact:T_tildes:" *
    Artifacts.artifact_meta("T_tildes", joinpath(@__DIR__, "..", "Artifacts.toml"))["git-tree-sha1"]
const LEGACY_BANK_SIGNATURE = bank_signature_digest(:chebyshev, _LEGACY_DATASET_ID, k_cheb)

function _load_legacy_bank()
    dir = joinpath(artifact"T_tildes", _LEGACY_BANK_NAME)
    load(name) = load_Ts(joinpath(dir, name), 128, 64, 161)[:, :, _R_KEEP_IDX, :]
    tensors = (T_2_00 = load("T_tilde_2_00"), T_0_00 = load("T_tilde_0_00"),
               T_minus2_00 = load("T_tilde_-2_00"), T_0_02 = load("T_tilde_0_02"),
               T_0_20 = load("T_tilde_0_20"), T_2_02 = load("T_tilde_2_02"),
               T_2_20 = load("T_tilde_2_20"), T_2_22 = load("T_tilde_2_22"))
    bank = _build_tensor_bank(tensors, k_cheb, :chebyshev, _LEGACY_DATASET_ID)
    bank.signature == LEGACY_BANK_SIGNATURE || error("legacy bank signature mismatch")
    return bank
end

const _LEGACY_BANK = Ref{TensorBank{Array{Float64,4}}}()
const _LEGACY_BANK_LOCK = ReentrantLock()

"""
    legacy_tensor_bank()

The validated fixed legacy Chebyshev bank, read from the `T_tildes` artifact on
first use only (never at module import) and shared by every later call. Loading
is serialized by a lock; once loaded, access is a lock-free read of immutable
data. The result does not depend on any differentiated quantity.
"""
function legacy_tensor_bank()
    isassigned(_LEGACY_BANK) && return _LEGACY_BANK[]
    lock(_LEGACY_BANK_LOCK) do
        isassigned(_LEGACY_BANK) || (_LEGACY_BANK[] = _load_legacy_bank())
    end
    return _LEGACY_BANK[]
end
