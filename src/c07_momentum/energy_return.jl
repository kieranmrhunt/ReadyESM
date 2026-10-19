module C07EnergyReturn

include("flux_form_cmt.jl")
const C = FluxFormCMT

"""One convex momentum substep with local kinetic-to-thermal energy return.

Transport the pre-mixing kinetic energy with the same linear mixing map.
The excess over resolved post-mixing kinetic energy becomes internal heat.
Evaluate in the lowest-level translating frame to avoid background-wind
cancellation. No clipping, uniform redistribution or energy correction.
"""
function exchange(u, v, mass, M, E, D, mass_flux, dt; strength=one(mass_flux), cp=1004)
    isfinite(cp) && cp > 0 || throw(ArgumentError("positive finite heat capacity required"))
    momentum = C.transport_momentum(u, v, mass, M, E, D, mass_flux, dt; strength)
    F = eltype(momentum.du)
    n = length(u)
    if iszero(momentum.effective_mass_flux)
        return merge(momentum, (; dT=zeros(F,n), transported_ke_tendency=zeros(F,n)))
    end
    a = F.(u) .- F(u[end])
    b = F.(v) .- F(v[end])
    kinetic = (a .^ 2 .+ b .^ 2) ./ F(2)
    plume_ke = C._plume_momentum(kinetic, F.(M), F.(E), F.(D))
    dk = C._environment_tendency(kinetic, plume_ke, F.(mass), F.(M), F.(E),
        F.(D), momentum.effective_mass_flux, F(dt))
    # Diagnose the finite-step work from the tendencies actually supplied.
    # Float64 arithmetic avoids subtracting rounded old/new total energies.
    dT = F[(Float64(dk[k]) - Float64(a[k])*Float64(momentum.du[k]) -
        Float64(b[k])*Float64(momentum.dv[k]) - Float64(dt)/2 *
        (Float64(momentum.du[k])^2 + Float64(momentum.dv[k])^2)) / Float64(cp)
        for k in 1:n]
    return merge(momentum, (; dT, transported_ke_tendency=dk))
end

function ledger(x, u, v, mass, dt, cp)
    m, du, dv = Float64.(mass), Float64.(x.du), Float64.(x.dv)
    work = Float64.(u).*du .+ Float64.(v).*dv .+ Float64(dt)/2 .* (du.^2 .+ dv.^2)
    heating = Float64(cp) .* Float64.(x.dT)
    residual = sum(m .* (work .+ heating))
    gross = sum(m .* (abs.(work) .+ abs.(heating)))
    return (; residual, gross, relative=abs(residual)/max(gross,floatmin(Float64)),
        heating_W_m2=sum(m.*heating), kinetic_loss_W_m2=-sum(m.*work),
        minimum_heating_K_s=minimum(x.dT))
end
end
