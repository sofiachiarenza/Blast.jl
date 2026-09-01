using DelimitedFiles
using Printf
using TOML

include(joinpath(@__DIR__, "..", "..", "characterization_workload.jl"))

function write_floats(path, values)
    open(path, "w") do io
        for value in vec(values)
            @printf(io, "%.17g\n", value)
        end
    end
end

function main()
    output_directory = @__DIR__
    result = run_characterization()
    arrays = characterization_arrays(result)
    shapes = Dict{String, Any}()
    for (name, values) in sort!(collect(arrays); by=first)
        write_floats(joinpath(output_directory, name * ".txt"), values)
        shapes[name] = collect(size(values))
    end

    metadata = Dict(
        "description" => "Legacy Blast marco_ad characterization before ACE/little-omega/state fixes",
        "git_commit" => readchomp(`git -C $(normpath(joinpath(@__DIR__, "..", "..", ".."))) rev-parse HEAD`),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blast_version" => string(Base.pkgversion(Blast)),
        "cosmology" => Dict(
            "ln10As" => result.cosmo.ln10Aₛ,
            "ns" => result.cosmo.nₛ,
            "h" => result.cosmo.h,
            "omegab" => result.cosmo.ωb,
            "omegac" => result.cosmo.ωc,
            "omegak" => result.cosmo.ωk,
            "mnu" => result.cosmo.mν,
            "w0" => result.cosmo.w0,
            "wa" => result.cosmo.wa,
        ),
        "pmm_contract" => "All supplied power spectra are total-matter Pmm in Mpc^3.",
        "shapes" => shapes,
    )
    open(joinpath(output_directory, "metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Wrote $(length(arrays)) text fixtures to $output_directory")
end

main()
