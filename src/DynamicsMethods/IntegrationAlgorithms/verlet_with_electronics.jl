
struct VerletwithElectronics <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

OrdinaryDiffEq.isfsal(::VerletwithElectronics) = false

mutable struct VerletwithElectronicsCache{uType,vType,rateType} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
end

function OrdinaryDiffEq.alg_cache(::VerletwithElectronics,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    VerletwithElectronicsCache(u, uprev, tmp, vtmp, k)
end

function OrdinaryDiffEq.initialize!(integrator, cache::VerletwithElectronicsCache)
    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    Calculators.update_electronics!(integrator.p.calculator, r)
    DynamicsMethods.SurfaceHoppingMethods.acceleration!(cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, cache::VerletwithElectronicsCache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, vtmp = cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)
    σfinal = DynamicsUtils.get_quantum_subsystem(u)

    step_B!(vtmp, vprev, dt/2, k)
    step_A!(rfinal, rprev, dt, vtmp)

    Calculators.update_electronics!(p.calculator, rfinal)
    DynamicsMethods.SurfaceHoppingMethods.acceleration!(k, vtmp, rfinal, p, t, p.method.state)

    step_B!(vfinal, vtmp, dt/2, k)

    DynamicsMethods.SurfaceHoppingMethods.propagate_wavefunction!(σfinal, σprev, vfinal, p, dt)

end

DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.SurfaceHoppingMethods.DiabaticIESH}) = VerletwithElectronics()
