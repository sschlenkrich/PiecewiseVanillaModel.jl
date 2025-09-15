
function model_string(m::Model)
    function slopes(ds, dv)
        ds0 = vcat(0.0, ds)
        dv0 = vcat(0.0, dv)
        return (dv0[2:end] .- dv0[1:end-1]) ./ (ds0[2:end] .- ds0[1:end-1])
    end
    rl = -slopes(m.dsl, m.dvl)
    ru = slopes(m.dsu, m.dvu)
    indent = "  "
    s = "PiecewiseVanillaModel:\n"
    s = s * indent
    # scalar parameters
    s = s * "s0: " * (@sprintf "%.4f" m.s0) * ", "
    s = s * "v0: " * (@sprintf "%.4f" m.v0) * ", "
    s = s * "w0: " * (@sprintf "%.4f" m.w0) * ", "
    s = s * "T: "  * (@sprintf "%.4f" m.T) * "\n"
    # header line
    s = s * indent * "        ds"
    s = s *  "                dv"
    s = s *  "                dw"
    s = s *  "                r\n"
    # lower wing table
    for (ds, dv, dw, r) in zip(reverse(m.dsl), reverse(m.dvl), reverse(m.dwl), reverse(rl))
        s = s * indent 
        s = s * (@sprintf "%16.8f" ds) * "  "
        s = s * (@sprintf "%16.8f" dv) * "  "
        s = s * (@sprintf "%16.8f" dw) * "  "
        s = s * (@sprintf "%16.8f" r) * "\n"
    end
    # delimiter
    s = s * indent * "----------------"
    s = s * "  ----------------"
    s = s * "  ----------------"
    s = s * "  ----------------\n"
    # upper wing table
    for (ds, dv, dw, r) in zip(m.dsu, m.dvu, m.dwu, ru)
        s = s * indent
        s = s * (@sprintf "%16.8f" ds) * "  "
        s = s * (@sprintf "%16.8f" dv) * "  "
        s = s * (@sprintf "%16.8f" dw) * "  "
        s = s * (@sprintf "%16.8f" r) * "\n"
    end
    # derived model parameters
    s = s * "Model parameters:\n"
    # header line
    s = s * indent * "        s"
    s = s * "                 σ"
    s = s * "                 λ" 
    s = s * "                 f\n"
    # lower wing
    sl = m.s0 .- vcat(0.0, m.dsl)
    for strike in reverse(sl)
        p = model_parameters(m, strike, -1)
        s = s * indent
        s = s * (@sprintf "%16.8f" strike) * "  "
        s = s * (@sprintf "%16.8f" p.σ) * "  "
        s = s * (@sprintf "%16.8f" p.λ) * "  "
        s = s * (@sprintf "%16.8f" p.f) * "\n"
    end
    # delimiter
    s = s * indent * "----------------"
    s = s * "  ----------------"
    s = s * "  ----------------"
    s = s * "  ----------------\n"
    # upper wing
    su = m.s0 .+ vcat(0.0, m.dsu)
    for strike in su
        p = model_parameters(m, strike, 1)
        s = s * indent
        s = s * (@sprintf "%16.8f" strike) * "  "
        s = s * (@sprintf "%16.8f" p.σ) * "  "
        s = s * (@sprintf "%16.8f" p.λ) * "  "
        s = s * (@sprintf "%16.8f" p.f) * "\n"
    end
    return s
end
