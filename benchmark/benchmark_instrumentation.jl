# Progress logging and an optional Linux CPU-temperature gate shared by the
# controlled and realistic benchmark runners. Minimal and specific on purpose.
#
# Thermal gate: set BLAST_BENCH_THERMAL_MAX_C (e.g. 60) to wait, before each
# timed section, until the hottest CPU sensor is below that threshold and has
# been stable (default 3 consecutive readings 2 s apart, max-min <= 2 C).
# BLAST_BENCH_THERMAL_TIMEOUT_S (default 600) bounds the wait. If requested but
# no sensor exists, or the timeout expires, the run fails. Unset = disabled and
# no behaviour change. The wait happens outside @benchmark.
using Dates

progress(message) = (println(Dates.format(now(), "HH:MM:SS"), " ", message); flush(stdout))

thermal_threshold() = (v = get(ENV, "BLAST_BENCH_THERMAL_MAX_C", ""); isempty(v) ? nothing : parse(Float64, v))

"Hottest CPU temperature in C from hwmon coretemp/k10temp sensors, or `nothing`."
function cpu_max_temp_c(root::AbstractString="/sys/class/hwmon")
    hottest = nothing
    isdir(root) || return nothing
    for entry in readdir(root; join=true)
        namefile = joinpath(entry, "name")
        isfile(namefile) && strip(read(namefile, String)) in ("coretemp", "k10temp") || continue
        for file in readdir(entry; join=true)
            occursin(r"/temp\d+_input$", file) || continue
            value = parse(Float64, strip(read(file, String))) / 1000
            hottest = isnothing(hottest) ? value : max(hottest, value)
        end
    end
    return hottest
end

function thermal_log(label; sampler=cpu_max_temp_c)
    isnothing(thermal_threshold()) && return nothing
    temperature = sampler()
    progress("thermal $label: cpu_max_c=$temperature wall=$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
    return temperature
end

function wait_for_thermal_stability(label; threshold=thermal_threshold(), sampler=cpu_max_temp_c,
                                    sleeper=sleep, clock=time, interval=2.0, required=3, drift=2.0,
                                    timeout=parse(Float64, get(ENV, "BLAST_BENCH_THERMAL_TIMEOUT_S", "600")))
    isnothing(threshold) && return nothing
    started = clock()
    window = Float64[]
    while true
        temperature = sampler()
        isnothing(temperature) && error("BLAST_BENCH_THERMAL_MAX_C is set but no coretemp/k10temp sensor was found")
        push!(window, temperature)
        length(window) > required && popfirst!(window)
        if length(window) == required && maximum(window) < threshold && maximum(window) - minimum(window) <= drift
            progress("thermal ready before $label: cpu_max_c=$temperature window=$window waited_s=$(round(clock() - started; digits=1)) wall=$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
            return temperature
        end
        clock() - started > timeout && error("thermal timeout before $label: cpu_max_c=$temperature threshold=$threshold after $(timeout) s")
        sleeper(interval)
    end
end

"Wait (if enabled) then log; call before each timed section."
function thermal_gate(label)
    wait_for_thermal_stability(label)
    return nothing
end
