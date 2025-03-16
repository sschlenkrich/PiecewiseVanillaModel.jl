
const Model = NamedTuple  # an alias

const r_ε = 1.0e-8  # r -> 0, avoid division by zero


function brownian_grid(v0, ds, dv)
    @assert length(ds) == length(dv)
    dw = zeros(0)
    for k=1:length(ds)
        ds0_ = (k==1) ? 0.0 : ds[k-1]
        dv0_ = (k==1) ? 0.0 : dv[k-1]
        dw0_ = (k==1) ? 0.0 : dw[k-1]
        #
        ds1_ = ds[k]
        dv1_ = dv[k]
        #
        r1 = (dv1_ - dv0_) / (ds1_ - ds0_)
        if abs(r1) < r_ε
            step = 1.0 / (v0 + dv0_)
        else
            step = (log(v0 + dv1_) - log(v0 + dv0_)) / (dv1_ - dv0_)
        end
        dw_k = dw0_ + (ds1_ - ds0_) * step
        dw = vcat(dw, dw_k)
    end
    return dw
end


function model(s0, v0, w0, T, dsl, dsu, dvl, dvu; rexl = 0.0, rexu = 0.0)
    @assert v0 > 0.0
    @assert T > 0.0
    @assert length(dsl) == length(dvl)
    @assert length(dsu) == length(dvu)
    # TODO: test monotonicity of S-grids, positivity of v
    dwl = brownian_grid(v0, dsl, dvl)
    dwu = brownian_grid(v0, dsu, dvu)
    return (
        s0  = s0,
        v0  = v0,
        w0  = w0,
        T   = T,
        dsl = dsl,
        dsu = dsu,
        dvl = dvl,
        dvu = dvu,
        dwl = dwl,
        dwu = dwu,
        rexl = rexl,
        rexu = rexu,
    )
end
