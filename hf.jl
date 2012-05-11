#!/usr/bin/env julia
abstract BasisFunction

#1s Gaussian-type function
type Gaussian1s <: BasisFunction
    #Gaussian orbital exponent
    alpha::Float64
    #nuclear coordinates
    center::Float64
end

#Slater Type Orbital fit with N primative gausians (STO-NG) type basis
type STONG <: BasisFunction
    n::Integer
    #contraction coeffiecents
    d::Array{Float64} 
    #primative gaussians
    g::Array{Gaussian1s}
end

#STO-3G basis for hydrogen
function sto3g_hydrogen(center)
    sto3g(center, 1.24)
end

#STO-3G basis for helium
function sto3g_helium(center)
    sto3g(center, 2.0925)
end

#Builds a STO-3G basis that best approximates a single slater type
#orbital with Slater orbital exponent zeta
function sto3g(center, zeta)
    scaling = zeta^2
    STONG(3,[0.444635, 0.535328, 0.154329], 
            [Gaussian1s(scaling*.109818, center),
             Gaussian1s(scaling*.405771, center),
             Gaussian1s(scaling*2.22766, center)])
end

#The overlap integrals describe how the basis functions overlap
#as the atom centered gaussian basis functions are non-orthognal
#they have a non-zero overlap. The integral has the following form:
#S_{ij} = \int \phi_i(r-R_a) \phi_j(r-R_b) \mathrm{d}r
function overlap_integral(b1::STONG, b2::STONG)
    two_center_contraction(b1, b2, overlap_integral)
end

#This function calculates the overlap integral for 1s Gaussian
#orbitals using a closed-form expression.
function overlap_integral(g1::Gaussian1s, g2::Gaussian1s)
    alpha = g1.alpha
    beta = g2.alpha
    Ra = g1.center
    Rb = g2.center

    #normalization constant
    n = (2*alpha/pi)^(3/4) * (2*beta/pi)^(3/4)

    I  = n * (pi/(alpha+beta))^(3/2) 
    I *= exp(-alpha*beta/(alpha+beta) * abs(Ra-Rb)^2)

    return I
end

#
function nuclear_attraction_integral(Zc::Int, Rc::Float64, b1::STONG, 
                                     b2::STONG)
    integral(g1, g2) = nuclear_attraction_integral(Zc, Rc, g1, g2)
    two_center_contraction(b1, b2, integral)
end

function nuclear_attraction_integral(Zc::Int, Rc::Float64, g1::Gaussian1s, g2::Gaussian1s)
    alpha = g1.alpha
    beta  = g2.alpha
    Ra = g1.center
    Rb = g2.center
    Rp = (alpha*Ra + beta*Rb)/(alpha + beta)

    n = (2*alpha/pi)^(3/4) * (2*beta/pi)^(3/4)
    matrix_element  = n*-2*pi/(alpha+beta)*Zc
    matrix_element *= exp(-alpha*beta/(alpha+beta)*abs(Ra-Rb)^2)

    t = (alpha+beta)*abs(Rp-Rc)^2
    if abs(t) < 1e-8
        return matrix_element
    end

    matrix_element *= 0.5 * sqrt(pi/t) * erf(sqrt(t))
    return matrix_element
end

function kinetic_energy_integral(b1::STONG, b2::STONG)
    two_center_contraction(b1, b2, kinetic_energy_integral)
end

function kinetic_energy_integral(g1::Gaussian1s, g2::Gaussian1s)
    alpha = g1.alpha
    beta = g2.alpha
    Ra = g1.center
    Rb = g2.center

    n = (2*alpha/pi)^(3/4) * (2*beta/pi)^(3/4)

    matrix_element  = n * alpha*beta/(alpha+beta)
    matrix_element *= (3-2*alpha*beta/((alpha+beta)/abs(Ra-Rb)^2 )) 
    matrix_element *= (pi/(alpha+beta))^(3/2)
    matrix_element *= exp(-alpha*beta/(alpha+beta) * abs(Ra-Rb)^2)
end

function two_electron_integral(g1::STONG, g2::STONG, g3::STONG, g4::STONG)
    four_center_contraction(g1, g2, g3, g4, two_electron_integral)
end

function two_electron_integral(g1::Gaussian1s, g2::Gaussian1s, g3::Gaussian1s, 
                               g4::Gaussian1s)
    alpha = g1.alpha
    beta  = g2.alpha
    gamma = g3.alpha
    delta = g4.alpha
    Ra = g1.center
    Rb = g2.center
    Rc = g3.center
    Rd = g4.center
    Rp = (alpha*Ra + beta*Rb)/(alpha + beta)
    Rq = (gamma*Rc + delta*Rd)/(gamma + delta)

    n  = (2*alpha/pi)^(3/4) * (2*beta/pi)^(3/4)
    n *= (2*gamma/pi)^(3/4) * (2*delta/pi)^(3/4)

    matrix_element  = n*2*pi^(5/2)
    matrix_element /= ((alpha+beta)*(gamma+delta)*sqrt(alpha+beta+gamma+delta))
    matrix_element *= exp(-alpha*beta/(alpha+beta)*abs(Ra-Rb)^2 - gamma*delta/(gamma+delta)*abs(Rc-Rd)^2)
    t = (alpha+beta)*(gamma+delta)/(alpha+beta+gamma+delta)*abs(Rp-Rq)^2
    if abs(t) < 1e-8
        return matrix_element
    end

    matrix_element *= 0.5 * sqrt(pi/t) * erf(sqrt(t))
    return matrix_element
end

function two_center_contraction(b1::STONG, b2::STONG, integral::Function)
    total = 0.0
    for p = 1:b1.n
        for q = 1:b2.n
            d1 = b1.d[p]
            d2 = b2.d[q]
            total += d1*d2*integral(b1.g[p], b2.g[q]) 
        end
    end
    return total
end

function four_center_contraction(b1::STONG, b2::STONG, b3::STONG, b4::STONG,
                                 integral::Function)
    total = 0.0
    for p = 1:b1.n
        for q = 1:b2.n
            for r = 1:b3.n
                for s = 1:b4.n
                    dp = b1.d[p]
                    dq = b2.d[q]
                    dr = b3.d[r]
                    ds = b4.d[s]
                    total += dp*dq*dr*ds*integral(b1.g[p], b2.g[q], b3.g[r], b4.g[s])
                end
            end
        end
    end
    total
 end


function hartree_fock(R2, P)
    R1 = 0.0
    R = [R1, R2]
    R12 = abs(R2-R1)

    Z1 = 2
    Z2 = 1
    Z = [Z1, Z2]

    #println("constructing basis set")
    phi1 = sto3g_helium(R1)
    phi2 = sto3g_hydrogen(R2)
    phi = [phi1, phi2]
    #calculate the overlap matrix S
    #the matrix should be symmetric with diagonal entries equal to one
    #println("building overlap matrix")
    S = eye(2)
    S[1,2] = overlap_integral(phi1, phi2)

    #calculate the overlap matrix S
    #the matrix should be symmetric with diagonal entries equal to one
    #println("building overlap matrix")
    S = eye(2)
    for i = 1:length(phi)
        for j = (i+1):length(phi)
            S[i,j] = S[j,i] = overlap_integral(phi[i], phi[j])
        end
    end

    #println("S: ", S)


    #calculate the kinetic energy matrix T
    #println("building kinetic energy matrix")
    T = zeros(2,2) 
    for i = 1:length(phi)
        for j = i:length(phi)
            T[i,j] = T[j,i] = kinetic_energy_integral(phi[i], phi[j])
        end
    end

    #println("T: ", T)

    #calculate nuclear attraction matrices V_i
    #println("building nuclear attraction matrices")
    V = zeros(length(Z), length(Z))
    for A = 1:length(Z)
        for i = 1:length(phi)
            for j = i:length(phi)
                v = nuclear_attraction_integral(Z[A], R[A], phi[i], phi[j])
                V[i,j] += v
                if i != j
                    V[j,i] += v
                end
            end
        end
    end

    #println("V: ", V)

    #build core-Hamiltonian matrix
    #println("building core-Hamiltonian matrix")
    Hcore = T + V

    #println("Hcore: ", Hcore)

    #diagonalize overlap matrix to get transformation matrix X
    #println("diagonalizing overlap matrix")
    s, U = eig(S)
    #println("building transformation matrix")
    X = U*diagm(s.^(-1/2))*U'
    #println("X: ", X)


    #calculate all of the two-electron integrals
    K = length(phi)
    two_electron = zeros(K,K,K,K)
    for mu = 1:2
        for v = 1:2
            for lambda = 1:2
                for sigma = 1:2
                    coulomb  = two_electron_integral(phi[mu], phi[v], 
                                                     phi[sigma], phi[lambda])
                    two_electron[mu,v,sigma,lambda] = coulomb
                    #println("coulomb  ($mu $v | $sigma $lambda): $coulomb")
                    exchange = two_electron_integral(phi[mu], phi[lambda], 
                                                     phi[sigma], phi[v])
                    #println("exchange ($mu $lambda | $sigma $v): $exchange")
                    two_electron[mu,lambda,sigma,v] = exchange
                end
            end
        end
    end

    total_energy = 0.0
    old_energy = 0.0
    electronic_energy = 0.0

    for scf_iter = 1:100
        #calculate the two electron part of the Fock matrix
        G = zeros(size(Hcore))
        for mu = 1:2
            for v = 1:2
                for lambda = 1:2
                    for sigma = 1:2
                        coulomb  = two_electron[mu,v,sigma,lambda]
                        exchange = two_electron[mu,lambda,sigma,v]
                        G[mu,v] += P[lambda,sigma]*(coulomb - 0.5*exchange)
                    end
                end
            end
        end

        #println("G: ", G)
        F = Hcore + G

        nuclear_energy = 0.0
        for A = 1:length(Z)
            for B = (A+1):length(Z)
                nuclear_energy += Z[A]*Z[B]/abs(R[A]-R[B])
            end
        end
        #println("E_nclr: $nuclear_energy")

        electronic_energy = 0.0
        for mu = 1:length(phi)
            for v = 1:length(phi)
                electronic_energy += P[v,mu]*(Hcore[mu,v]+F[mu,v])
            end
        end
        electronic_energy *= 0.5
        #println("E_elec: $electronic_energy")
        total_energy = electronic_energy + nuclear_energy
        printf("SCF i: %3i E_elec: %12.8f E_total: %12.8f de: %12.4e\n", 
               scf_iter, electronic_energy, total_energy, 
               total_energy - old_energy)

        if scf_iter > 2 && abs(old_energy - total_energy) < 1e-6
            break 
        end

        #println("F: ", F)
        Fprime = X' * F * X
        #println("F': $Fprime")
        epsilon, Cprime = eig(Fprime)
        #println("epsilon: ", epsilon)
        #println("C': ", Cprime)
        C = real(X*Cprime)
        #println("C: ", C)

        P = zeros(size(Hcore))
        for mu = 1:length(phi)
            for v = 1:length(phi)
                P[mu,v] = 2*C[mu,1]*C[v,1]
            end
        end
        #println("P: ", P)


        old_energy = total_energy
    end

    printf("r12: %.4f e_tot: %12.8f e_elec: %12.8f\n", R12, total_energy, 
           electronic_energy)

    return total_energy, electronic_energy, P
end

P = zeros(2,2)
total_energy, electronic_energy, P = hartree_fock(1.4632, P)
szabo_energy = -4.227529
if abs(electronic_energy - szabo_energy) > 1e-6
    println("TEST FAILED")
else
    println("TEST PASSED")
end
