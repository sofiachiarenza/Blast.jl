# Temporary bootstrap: Blast's B-spline backend needs the unreleased
# AbstractCosmologicalEmulators 0.12 (branch 96-b-spline), which is not in the
# General registry. Run it with the project to be prepared, before any
# instantiate/build/test step:
#
#     julia --project=.     scripts/bootstrap_ace_dev.jl   # Blast's own project
#     julia --project=docs/ scripts/bootstrap_ace_dev.jl   # docs project
#
# The pin is explicit and does not rely on Project.toml `[sources]`, which
# Julia 1.10 ignores. The docs project declares the unregistered Blast itself,
# so there both unregistered sources must be given to a single resolution:
# the local checkout and a clone of ACE at the pinned revision (under
# `<project>/dev/AbstractCosmologicalEmulators-<rev>`, since `Pkg.develop` has no `rev` option) are developed
# together in one call. Requires the `git` command line tool.
#
# Remove this script and its CI steps once a registered ACE release containing
# these APIs replaces it.
using Pkg, TOML

const ACE_URL = "https://github.com/CosmologicalEmulators/AbstractCosmologicalEmulators.jl"
const ACE_REV = "af7944975f422af0271ebdfb03a6ba2395a5e2c8"
const BLAST_UUID = "aee4e916-1e40-4aff-8d35-b84b50e5d5e3"

deps = get(TOML.parsefile(Base.active_project()), "deps", Dict{String,Any}())

if get(deps, "Blast", nothing) == BLAST_UUID
    clone = joinpath(dirname(Base.active_project()), "dev", "AbstractCosmologicalEmulators-" * ACE_REV[1:7])
    if ispath(clone)
        # Reuse only an untouched clone of exactly the pinned revision; never
        # reset, clean, stash or check out anything in an existing directory.
        git(args...) = readchomp(`git -C $clone $args`)
        problems = String[]
        try
            realpath(git("rev-parse", "--show-toplevel")) == realpath(clone) ||
                push!(problems, "not the root of its own git repository")
            git("rev-parse", "HEAD") == ACE_REV || push!(problems, "HEAD is not $ACE_REV")
            isempty(git("status", "--porcelain", "--untracked-files=all")) ||
                push!(problems, "working tree has tracked or untracked changes")
        catch err
            push!(problems, "not a usable git clone ($(sprint(showerror, err)))")
        end
        isempty(problems) || error("Refusing to reuse $clone: ", join(problems, "; "),
                                   ". Move it away or fix it manually; nothing was modified.")
    else
        run(`git clone --quiet $ACE_URL $clone`)
        run(`git -C $clone checkout --quiet $ACE_REV`)
    end
    Pkg.develop([PackageSpec(path=dirname(@__DIR__)), PackageSpec(path=clone)])
else
    Pkg.add(PackageSpec(name="AbstractCosmologicalEmulators", url=ACE_URL, rev=ACE_REV))
    Pkg.resolve()
end
