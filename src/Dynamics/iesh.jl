export IESHPhasespace
export IESH
export IESH_callback
using Revise
using DiffEqCallbacks
using Combinatorics

"""This module controles how IESH is executed. For a description of IESH, see e.g.
Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009) and 
Shenvi, Roy,  Tully,J. Chem. Phys. 130, 174107 (2009).
It first needs to be initialized, e.g.
dynam = Dynamics.IESH{Float64}(DoFs, atoms, n_states)
After setting up the simulation (e.g. sim = Simulation(atoms, model, dynam; DoFs=1)),
for a single trajectory, the Phasespace is set up inside this module:
z = SurfaceHoppingPhasespace(r,p, n_states+2, k)  
and the trajectory can be run with:
solution = Dynamics.run_trajectory(z, (0.0, 15000.0), sim)
The latter calls back to DifferentialEquations, but the callback of DifferentialEquations
with the hopping between surfaces is handled here.
"""

struct IESH{T} <: Method
    nonadiabatic_coupling::Matrix{Matrix{T}}
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    momentum_rescale::Vector{T}
    function IESH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
        nonadiabatic_coupling = [zeros(states, states) for i=1:DoFs, j=1:atoms]
        density_propagator = zeros(states, states)
        #hopping_probability = zeros(states)
        hopping_probability = zeros(3)
        momentum_rescale = zeros(atoms)
        new{T}(nonadiabatic_coupling, density_propagator, hopping_probability, momentum_rescale)
    end
end

mutable struct IESHPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{Matrix{T}, Matrix{T}, Matrix{Complex{T}}}}
    state::Vector{Int}
end

function IESHPhasespace(x::ArrayPartition{T}, state::Vector{Int}) where {T<:AbstractFloat}
    IESHPhasespace{T}(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function IESHPhasespace(R::Matrix{T}, P::Matrix{T}, σ::Matrix{Complex{T}}, state::Vector{Int}) where {T}
    IESHPhasespace{T}(ArrayPartition(R, P, σ), state)
end

# construct the density matrix
# Adapted from James's definition. We now need states as the Vector |k>
# See: ShenviRoyTully_JChemPhys_130_174107_2009
function IESHPhasespace(R::Matrix{T}, P::Matrix{T}, n_states::Integer, state::Vector{Int}) where {T}
    σ = zeros(Complex{T}, n_states, n_states)
    # 1 needs to be replaced by the info from state, once I've gotten things to run.
    c = 0
    for i=1:length(state)/2
        c = c + 1
        σ[c, c] = 1
    end
    IESHPhasespace(R, P, σ, state)
end

get_density_matrix(z::IESHPhasespace) = z.x.x[3]

# This belongs to run_trajectory.
# It first follows motion! and then iesh_callback
# Doesn't work, yet. HERE!
function create_problem(u0::IESHPhasespace, tspan::Tuple, sim::AbstractSimulation{<:IESH})
    #IESH_callback=create_energy_saving_callback(u0,integrator::DiffEqBase.DEIntegrator)
    # IESH_callback=create_saving_callback()
    cb2 = DiscreteCallback(condition, affect!; save_positions=(false, false))
    ODEProblem(motion!, u0, tspan, sim; callback=cb2)
 end


function motion!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH}, t)
    println("ping1 ", time())
    #println(get_positions(u))
    set_velocity!(du, u, sim)  # classical.jl momenta/(atom_mass) v = p/m 
    println("ping2 ", time())
    update_electronics!(sim, u) # The next three routines 
    println("ping3 ", time())
    set_force!(du, u, sim) # 4th routine
    println("ping4 ", time())
    set_density_matrix_derivative!(du, u, sim)# 5th routine
    println("ping5 ", time())
end

# Get the current eneries, the current derivate and the current eigenvectors
# see ../Calculators/Calculators.jl for the functions used here.
# Then get the eigenvectors of the energy
# then use them to transform the derivatives into the adiabatic form.
function update_electronics!(sim::Simulation{<:IESH}, u::IESHPhasespace)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    Calculators.transform_derivative!(sim.calculator)
    evaluate_nonadiabatic_coupling!(sim)
end

#These routines point to the previous routine and hands down an empty array,
# the adiabatic_derivatives calculated in transform_drivatives and the 
# energy eigenvectors
function evaluate_nonadiabatic_coupling!(sim::Simulation{<:IESH})
    evaluate_nonadiabatic_coupling!.(
        sim.method.nonadiabatic_coupling,
        sim.calculator.adiabatic_derivative,
        Ref(sim.calculator.eigenvalues))
end

# Calculates the nonadiabatic coupling. This is probably also true for IESH
function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
end

# gets the momenta corresponding to the current state from the adiabatic derivates 
# (transformation of the diabatic ones, see 3 routines above)
# Here, instead of just calculating the moment, I assume I will need to some up (?)
# the different contributions to the force according to Eq.(12) in the Tully paper.
function set_force!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH})
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            get_momenta(du)[j,i] = 0.0
            for n=1:length(u.state)
                get_momenta(du)[j,i] = get_momenta(du)[j,i] 
                                       #-sim.calculator.adiabatic_derivative[j,i][u.state, u.state]
                                    -sim.calculator.adiabatic_derivative[j,i][n, n]
            end
        end
    end
end

# Calculate the _time_-derivate  of the  density matrix
# Goes back to motion!
# This one should be the same as for FSSH
function set_density_matrix_derivative!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH})
    σ = get_density_matrix(u)
    velocity = get_positions(du)
    V = sim.method.density_propagator # this is the potential energy surface

    V .= diagm(sim.calculator.eigenvalues)# Creates a diagonal matrix from eigenvalues

    # Calculation going on here from: Martens_JPhysChemA_123_1110_2019, eq. 6
    # d is the nonadiabatic coupling matrix
    #i ħ dσ/dt = iħ sum_l [(V_{m,l} - i ħ velocity d_{m,l})*σ_{l,n} - &
    #                      σ_{m,l}*(V_{l,n} - iħ velocity d_{l,n})]
    # vgl. also SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, eq. 5 and 6
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            V .-= im*velocity[j,i].*sim.method.nonadiabatic_coupling[j,i]
        end
    end
    get_density_matrix(du) .= -im*(V*σ - σ*V)
 end

# This sets when the condition will be true. In this case, always.
condition(u, t, integrator::DiffEqBase.DEIntegrator) = true

function affect!(integrator::DiffEqBase.DEIntegrator)
    println("ping6 ", time())
    update_hopping_probability!(integrator)
    println("ping7 ", time())
    
    # not necessary anymore
    #new_state = select_new_state(integrator.p.method.hopping_probability, integrator.u.state)
    
    #if new_state != 0
    if integrator.p.method.hopping_probability[1] !=0
        # Set new state population
        new_state = copy(integrator.u.state)
        new_state[Int(integrator.p.method.hopping_probability[2])] = 0
        new_state[Int(integrator.p.method.hopping_probability[3])] = 1

        # needs probably a vector as input for new_state (i.e. the state distribution)
        if calculate_rescaling_constant!(integrator, new_state)
            execute_hop!(integrator, new_state)
        end
    end
end

function update_hopping_probability!(integrator::DiffEqBase.DEIntegrator)
    sim = integrator.p
    coupling = sim.method.nonadiabatic_coupling
    velocity = get_positions(get_du(integrator))
    s = integrator.u.state
    σ = get_density_matrix(integrator.u)
    dt = get_proposed_dt(integrator)
    
    sim.method.hopping_probability .= 0 # Set all entries to 0
    hop_mat = zeros(length(s),length(s))
    sumer=0
    sum_before = 0
    random_number = rand()
    first = true
    #random_number = 1e-8
    # Calculate matrix with hopping probs
    # This should probably be optimized at some point.
    for l = 1:length(s)
        # If occupied
        if(integrator.u.state[l]==1)
            for m = 1:length(s)
                # if unoccupied
                if(integrator.u.state[m]==0)
                    for i=1:length(sim.atoms)
                        for j=1:sim.DoFs
                            #sim.method.hopping_probability[m] += 2*velocity[j,i]*real(σ[m,s]/σ[s,s])*coupling[j,i][s,m] * dt
                            hop_mat[l,m] += 2*velocity[j,i]*real(σ[m,l]/σ[l,l])*coupling[j,i][l,m] * dt
                        end
                    end
                end    
                clamp(hop_mat[l,m], 0, 1)
                # Calculate the hopping probability. Hopping occures for
                # the transition that's first above the random number.
                # See: Tully_JChemPhys_93_1061_1990
                sumer = sumer+abs(hop_mat[l,m]) # cumulative sum.
                if (random_number > sumer)
                    sum_before = sumer
                elseif (random_number < sumer && random_number > sum_before && first)
                    sim.method.hopping_probability[1] = sumer
                    sim.method.hopping_probability[2] = l
                    sim.method.hopping_probability[3] = m
                    first = false
                elseif (sumer > 1)
                    println("Error: Sum of hopping probability above 1!")
                    println("Sum: ", sumer, " Individ. hopping probability: ", hop_mat[l,m])
                    println("l = ", l, " m = ", m)
                    exit()
                end
            end
        end
    end
    
    #a=findmax(hop_mat)
    # Write the hopping probability and the array positions into array
    # This one just extracts the maximum hopping probability. 
    # May be alternative to above
    #sim.method.hopping_probability[1] = maximum(hop_mat)
    #sim.method.hopping_probability[2] = a[2][1]
    #sim.method.hopping_probability[3] = a[2][2]

    #@time a=combinations(integrator.u.state)
    #b=collect(permutations(integrator.u.state))
    #println(b)

    # Old version from FSSH
    # for m=1:integrator.p.calculator.model.n_states
    #     if m != s
    #         for i=1:length(sim.atoms)
    #             for j=1:sim.DoFs
    #                 sim.method.hopping_probability[m] += 2*velocity[j,i]*real(σ[m,s]/σ[s,s])*coupling[j,i][s,m] * dt
    #             end
    #         end
    #     end
    # end
    #clamp!(sim.method.hopping_probability[1], 0, 1) # Restrict probabilities between 0 and 1
    # A cumulative sum. No idea why.
    #cumsum!(sim.method.hopping_probability[1], sim.method.hopping_probability[1]) # Cumulative sum
end

# Not necessary anymore, bcs hopping probabilty selected inside loop.
# function select_new_state(probability::Vector{T}, current_state::Integer)::UInt where {T<:AbstractFloat}
#     random_number = rand()
#     for (i, prob) in enumerate(probability)
#         if i != current_state # Avoid self-hops
#             if prob > random_number
#                 return i # Return index of selected state
#             end
#         end
#     end
#     0 # Return 0 if no hop is desired
# end

# Calculate if hop frustrated
# This includes a momentum rescaling (IESH paper of Shenvi et al. expresses it in terms of velocity)
# I believe (also according to Subotnic&Miao JCP 150 2019) that this is the same.
# In any case, it's used to conserve energy. Reini remarked that we might not want it 
# eventually, but I'm leaving it for now.
# It should be related to: HammesSchifferTully_JChemPhys_101_4657_1994
function calculate_rescaling_constant!(integrator::DiffEqBase.DEIntegrator, new_state)::Bool
    sim = integrator.p
    old_state = integrator.u.state
    state_diff = integrator.p.method.hopping_probability
    velocity = get_positions(get_du(integrator))

    # Loop over and sum over optential energies, according to Eq. 12 Shenvi, Roy,  Tully,J. Chem. Phys. 130, 174107 (2009)
    c = 0
    for i=1:length(old_state)
        c = c + calculate_potential_energy_change(sim.calculator.eigenvalues, 
                                                  new_state[i], old_state[i],i)
    end
    a = zeros(length(sim.atoms))
    b = zero(a)
    # view: treats data structure from array as another array
    #': conjucated transposition (adjoint)
    # Might be that old and new state need to be other way around; changes sign.
    # (but probably  not; looks like any sign vanishes if real valued)
    @views for i in range(sim.atoms)
        coupling = [sim.method.nonadiabatic_coupling[j,i][Int(state_diff[3]), Int(state_diff[2])] for j=1:sim.DoFs]
        #coupling = [sim.method.nonadiabatic_coupling[j,i][Int(state_diff[2]), Int(state_diff[3])] for j=1:sim.DoFs]
        a[i] = coupling'coupling / sim.atoms.masses[i]
        b[i] = velocity[:,i]'coupling
    end
    
    discriminant = b.^2 .- 2a.*c
    if any(discriminant .< 0)
        return false
    else
        root = sqrt.(discriminant)
        integrator.p.method.momentum_rescale .= min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
    return true
    end
end

# This does to update_electronics above.
function calculate_potential_energy_change(eigenvalues::Vector, new_state::Integer, 
    current_state::Integer, counter::Integer)
    new_state*eigenvalues[counter] - current_state*eigenvalues[counter]
end

function execute_hop!(integrator::DiffEqBase.DEIntegrator, new_state::Vector)
    state_diff = integrator.p.method.hopping_probability
    # For momentum rescaling, see eq. 7 and 8 SubotnikBellonzi_AnnuRevPhyschem_67_387_2016
    for i in range(integrator.p.atoms)
        coupling = [integrator.p.method.nonadiabatic_coupling[j,i][Int(state_diff[3]), 
                   Int(state_diff[2])] for j=1:integrator.p.DoFs]
        #for n=1:length(new_state)
        #coupling = [integrator.p.method.nonadiabatic_coupling[j,i][new_state, integrator.u.state] for j=1:integrator.p.DoFs]
        #get_momenta(integrator.u)[:,i] .-= integrator.p.method.momentum_rescale[i] .* coupling
        # No sum over states necessary, because momenta just rescaled & 
        # not calculated from scratch
        get_momenta(integrator.u)[:,i] .-= integrator.p.method.momentum_rescale[i] .* coupling
        #end
    end
    integrator.u.state = new_state
end