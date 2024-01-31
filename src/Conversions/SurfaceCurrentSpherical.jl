# function convertrepresentation(T::Type{<:AbstractSphericalExpansion},currents::AbstractSurfaceCurrentDensity, L::Integer, k₀::Number)
#     C=SphericalVectorModeFields.elementtype(T)
#     αvec=zeros(C, sℓm_to_j(2,L,L))

#     Threads.@threads for j in eachindex(αvec)
#         s, ℓ, m = j_to_sℓm(j)
#         αtmp=zeros(C, sℓm_to_j(2,L,L))
#         αtmp[sℓm_to_j(s, ℓ, -m)]=(-1)^(m+1)
#         αvec[j]=2*transmission(reciprocaltype(T)(αtmp), currents, k₀)
#     end
#     return T(αvec)
# end

function convertrepresentation(T::Type{<:RadiatingSphericalExpansion}, currents::AbstractSurfaceCurrentDensity, k₀::Number)
    ϵ = 1e-7

    rmax = 0.0
    for pos in currents.functionspace.geo.vertices
        if norm(pos) > rmax
            rmax = norm(pos)
        end
    end
    rmax *= 2
    L = maximum([3, Int(ceil(k₀ * rmax + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k₀ * rmax)^(1 / 3)))])
    return convertrepresentation(T, currents, L, k₀)
end
function convertrepresentation(T::Type{<:IncidentSphericalExpansion}, currents::AbstractSurfaceCurrentDensity, k₀::Number)
    ϵ = 1e-7

    rmin = Inf
    for dipole in dipoles
        if norm(dipole.pos) < rmin
            rmin = norm(dipole.pos)
        end
    end
    L = Int(ceil((k₀ * rmin + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k₀ * rmin)^(1 / 3))))
    return convertrepresentation(T, currents, L, k₀)
end


function convertrepresentation(T::Type{<:AbstractSphericalExpansion}, currents::ElmagSurfaceCurrentDensity, L::Integer, k₀::Number)
    electric_currents = ElectricSurfaceCurrentDensity(currents.functionspace, currents.electricexcitations)
    magnetic_currents = MagneticSurfaceCurrentDensity(currents.functionspace, currents.magneticexcitations)
    α1 = convertrepresentation(T, electric_currents, L, k₀)
    α2 = convertrepresentation(T, magnetic_currents, L, k₀)
    return T(α1.coefficients + α2.coefficients)
end


#######################################
# Fast implementation but relies on internals of BEAST
#  |
#  V
######################################
function convertrepresentation(T::Type{<:AbstractSphericalExpansion}, currents::AbstractSurfaceCurrentDensity, L::Integer, k₀::Number)
    C = SphericalVectorModeFields.elementtype(T)
    Jmax = sℓm_to_j(2, L, L)
    αvec = zeros(C, Jmax)
    currenttype = typeof(currents)

    fieldop = _sphericalwavefieldoperator(currenttype, reciprocaltype(T), Jmax, k₀)
    multifield = MultiFunctional(fieldop, Jmax)

    tested_incident_field = BEAST.assemble(multifield, currents.functionspace)
    return T(k₀ * _fieldfactor(currenttype) * sum(tested_incident_field .* currents.excitations))
end

function _sphericalwavefieldoperator(
    currenttype::Type{MagneticSurfaceCurrentDensity{R}}, sphericaltype::Type{<:AbstractSphericalExpansion}, Jmax::Integer, k₀::Number
) where {R}
    return (x -> curlF_sℓm_cartesian_array(Jmax, sphericaltype, x, k₀))
end
function _sphericalwavefieldoperator(
    currenttype::Type{ElectricSurfaceCurrentDensity{R}}, sphericaltype::Type{<:AbstractSphericalExpansion}, Jmax::Integer, k₀::Number
) where {R}
    return (x -> F_sℓm_cartesian_array(Jmax, sphericaltype, x, k₀))
end
function _fieldfactor(currenttype::Type{MagneticSurfaceCurrentDensity{R}}) where {R}
    return (1 / sqrt(Z₀))
end
function _fieldfactor(currenttype::Type{ElectricSurfaceCurrentDensity{R}}) where {R}
    return (sqrt(Z₀))
end

mutable struct MultiFunctional <: BEAST.Functional
    field
    valuecount::Integer
end

function (MF::MultiFunctional)(p)
    F = MF.field
    x = cartesian(p)
    return F(x)
end

function BEAST.assemble(multifield::MultiFunctional, tfs; quaddata=BEAST.quaddata, quadrule=BEAST.quadrule)

    R = scalartype(tfs)
    b = fill(zeros(Complex{R}, multifield.valuecount), numfunctions(tfs))
    store(v, m) = (b[m] .+= v)
    BEAST.assemble!(multifield, tfs, store; quaddata=quaddata, quadrule=quadrule)
    return b
end


function BEAST.celltestvalues(tshs::BEAST.RefSpace{T,NF}, tcell, field::MultiFunctional, qr) where {T,NF}

    num_tshs = BEAST.numfunctions(tshs)
    interactions = fill(zeros(Complex{T}, field.valuecount), num_tshs)
    storeinteractions(w, m) = (interactions[m] .+= w)

    num_oqp = length(qr)

    for p in 1:num_oqp
        mp = BEAST.cartesian(qr[p].point)

        dx = qr[p].weight

        Fx, Fy, Fz = field.field(mp)
        Fvec = [Fx Fy Fz]
        tvals = qr[p].value

        for m in 1:num_tshs
            tval = tvals[m]
            # interactions[:,m] += Fvec* tval[1] *dx
            storeinteractions(Fvec * tval[1] * dx, m)
        end
    end

    return interactions
end

# function assemble!(field::MultiFunctional, tfs::subdBasis, store;
#     quaddata=quaddata, quadrule=quadrule)

#     tels, tad = assemblydata(tfs)

#     trefs = refspace(tfs)
#     qd = quaddata(field, trefs, tels)

#     for (t, tcell) in enumerate(tels)

#         # compute the testing with the reference elements
#         qr = quadrule(field, trefs, t, tcell, qd)
#         blocal = celltestvalues(trefs, tcell, field, qr)

#         for i in 1 : length(tad[t])
#             for (m,a) in tad[t][i]
#                 store(a*blocal[:,i], m)
#             end
#         end

#     end

# end

# function convertrepresentation(T::Type{<:AbstractSphericalExpansion},receivecurrents::ElectricSurfaceCurrentDensity, L::Integer, k₀::Number)
#     C=SphericalVectorModeFields.elementtype(T)
#     Jmax=sℓm_to_j(2,L,L)
#     efieldfun=x->efield(reciprocaltype(T)(zeros(C, Jmax)), x, k₀)
#     tangential_trace=((BEAST.n × (efieldfun)) × BEAST.n)

#     tels, tad = assemblydata(receivecurrents.functionspace)
#     trefs = refspace(receivecurrents.functionspace)
#     qd = quaddata(tangential_trace, trefs, tels)

#     storevec = zeros(C, Jmax, numfunctions(receivecurrents.functionspace))
#     store(v,m)= (storevec[:,m]+=v)

#     for (t, tcell) in enumerate(tels)
#         qr = quadrule(tangential_trace, trefs, t, tcell, qd)
#         num_tshs = numfunctions(trefs)
#         interactions = zeros(Complex{Float64},Jmax, num_tshs)
#         storeinteraction(w,m)=(interactions[:,m]+=w)

#         num_oqp = length(qr)

#         for p in 1 : num_oqp
#             mp = BEAST.cartesian(qr[p].point)

#             dx =qr[p].weight

#             Fx, Fy, Fz = F_sℓm_cartesian_array(Jmax, SphericalVectorModeFields.IncidentSphericalExpansion{ComplexF64}, mp, k₀)

#             Fvec = [Fx Fy Fz]  * k₀ * sqrt(Z₀)    
#             tvals = qr[p].value

#             for m in 1 : num_tshs
#                 tval = tvals[m]
#                 # interactions[:,m]+= Fvec* tval[1] *dx
#                 storeinteraction(Fvec* tval[1] *dx, m)
#             end
#         end


#         for i in 1 : length(tad[t])
#             for (m,a) in tad[t][i]
#                 store(a*interactions[:,i],m)
#             end
#         end
#     end
#     α_hlp=storevec* receivecurrents.excitations
#     α=zeros(ComplexF64, size(α_hlp))
#     for j = 1:length(α_hlp)
#         s,ℓ, m = j_to_sℓm(j)
#         jj = sℓm_to_j(s,ℓ, -m)
#         α[jj] = (-1)^(m + 1) * α_hlp[j] 
#     end
#     return α
# end

# function convertrepresentation(T::Type{<:AbstractSphericalExpansion},receivecurrents::MagneticSurfaceCurrentDensity, L::Integer, k₀::Number)
#     C=SphericalVectorModeFields.elementtype(T)
#     Jmax=sℓm_to_j(2,L,L)
#     # fieldfun=x->efield(reciprocaltype(T)(zeros(C, Jmax)), x, k₀)
#     fieldfun=x->curlF_sℓm_cartesian_array(Jmax, SphericalVectorModeFields.IncidentSphericalExpansion{C}, x, k₀)
#     tangential_trace=((BEAST.n × (fieldfun)) × BEAST.n)

#     tels, tad = assemblydata(receivecurrents.functionspace)
#     trefs = refspace(receivecurrents.functionspace)
#     qd = quaddata(tangential_trace, trefs, tels)

#     storevec = zeros(C,Jmax, numfunctions(receivecurrents.functionspace))
#     store(v,m)= (storevec[:,m]+=v)

#     for (t, tcell) in enumerate(tels)
#         qr = quadrule(tangential_trace, trefs, t, tcell, qd)
#         num_tshs = numfunctions(trefs)
#         interactions = zeros(Complex{Float64}, Jmax, num_tshs)
#         storeinteraction(w,m)=(interactions[:,m]+=w)

#         num_oqp = length(qr)

#         for p in 1 : num_oqp
#             mp = BEAST.cartesian(qr[p].point)

#             dx =qr[p].weight

#             Fx, Fy, Fz = curlF_sℓm_cartesian_array(Jmax, SphericalVectorModeFields.IncidentSphericalExpansion{C}, mp, k₀)

#             Fvec = [Fx Fy Fz]  *complex(0.0,-k₀) / sqrt(Z₀)    
#             tvals = qr[p].value

#             for m in 1 : num_tshs
#                 tval = tvals[m]
#                 # interactions[:,m]+= Fvec* tval[1] *dx
#                 storeinteraction(Fvec* tval[1] *dx, m)
#             end
#         end


#         for i in 1 : length(tad[t])
#             for (m,a) in tad[t][i]
#                 # storevec[:,m]+=a*interactions[:,i]
#                 store(a*interactions[:,i],m)
#             end
#         end
#     end
#     α_hlp=storevec* receivecurrents.excitations
#     α=zeros(ComplexF64, size(α_hlp))
#     for j = 1:length(α_hlp)
#         s,ℓ, m = j_to_sℓm(j)
#         jj = sℓm_to_j(s,ℓ, -m)
#         α[jj] = (-1)^(m + 1) * α_hlp[j] 
#     end
#     return α
# end

#######################################
#  ^
#  |
#  Fast implementation but relies on internals of BEAST
######################################


