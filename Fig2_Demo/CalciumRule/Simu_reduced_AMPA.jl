# Defines output directory
directory_name = "/Users/kathleen/Documents/GitHub/Thesis/Julia/CalciumRule"

# Loads packages
using Plots
using DelimitedFiles
using Statistics
using Images, ImageView
using DifferentialEquations
using DataFrames
using Printf
using CSV
using LinearAlgebra
using Distributions
using DataStructures

## Include model

include("model_reduced_AMPA.jl")
const SD = "no"
const new_network=1


if(SD=="yes")
    experiment_name = "AMPA_SD"
else
    experiment_name = "AMPA"
end
idx_ntk = 1

# Network parameters
const ncellsI = 1
const nPre = 2
const nPost = 1
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC


## Simulation parameters
const dt = 0.01

const N_cycles=2
const Duration_cycle = 30000
const Tdt_cycle = convert(Int64, Duration_cycle/dt)

const N_patterns = 1
const Duration_state = convert(Int64, Duration_cycle/2)
const Tdt_state = convert(Int64, Duration_state/dt)

const N_samples = 1
const Duration_sample = convert(Int64, Duration_state/N_samples)
const Tdt_sample = convert(Int64, Duration_sample / dt)

const Duration_testing = 1
const Tdt_testing = convert(Int64, Duration_testing/dt)

const Duration_transient_testing = 1
const Tdt_transient_testing = convert(Int64, Duration_transient_testing/dt)

const T = N_cycles*Duration_cycle + Duration_testing*N_patterns*N_cycles + Duration_transient_testing*N_cycles
const Tdt = convert(Int64, T / dt)
const t = range(dt, T, length = Tdt)

const Duration_set = Duration_cycle + Duration_testing*N_patterns + Duration_transient_testing
const Tdt_set = convert(Int64, Duration_set/dt)

const state = zeros(Tdt,1)
for idx_cycle=1:1:N_cycles
    # wake
    state[(idx_cycle-1)*Tdt_set+1: (idx_cycle-1)*Tdt_set+Tdt_state] .= 0
    # sleep
    if(SD=="yes")
        state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= -1
    else
        state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= 1
    end
    # transient testing
    TA = (idx_cycle-1)*Tdt_set+Tdt_cycle
    state[TA+1: TA+Tdt_transient_testing] .= 2
    # test
    TB= TA+ Tdt_transient_testing
    state[TB+1: TB+Tdt_testing*N_patterns] .= 3
end


## Neurons' model parameters

# Global parameters
const C = 1
const VNa = 50
const VK = -85
const VCa = 120
const Vl = -55
const VH = -20
const Kd = 170

# Cells parameters
const gl = 0.055
const gNa = 170.0
const gKd = 40
const k1 = 1.e-1
const k2 = 0.1e-1
const gH = 0.01
const gKCa = 4
const gCaT = 0.55

const gamma = 0.0 #1  #10% of variability in the network put 0.1 to get 10% etc
const gl_cells = rand((gl-(gamma*gl):0.001:gl+(gamma*gl)),ncells)
const gNa_cells = rand((gNa-(gamma*gNa):0.001:gNa+(gamma*gNa)),ncells)
const gKd_cells = rand((gKd-(gamma*gKd):0.001:gKd+(gamma*gKd)),ncells)
const k1_cells = rand((k1-(gamma*k1):0.001:k1+(gamma*k1)),ncells)
const k2_cells = rand((k2-(gamma*k2):0.001:k2+(gamma*k2)),ncells)
const gH_cells = rand((gH-(gamma*gH):0.001:gH+(gamma*gH)),ncells)
const gKCa_cells = rand((gKCa-(gamma*gKCa):0.001:gKCa+(gamma*gKCa)),ncells)
const gCaT_cells = rand((gCaT-(gamma*gCaT):0.001:gCaT+(gamma*gCaT)),ncells)


## Current parameters
const IappI = 3.
const IappC = 0.

const spike_duration = 3
const IstepI = -1.2-IappI
const IstepC = 50.

IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
Istep_cell = [IstepI_cell; IstepC_cell]


if(new_network==1)
    neurons_freq  = zeros(N_samples*N_cycles, ncells)   
    neurons_freq[:,1] .= 1 # inhibitory cell

    neurons_freq[1,2] = 60 # associative input
    neurons_freq[1,3] = 1
    neurons_freq[1,4] = 10

    neurons_freq[2,2] = 55 # associative input
    neurons_freq[2,3] = 0.5
    neurons_freq[2,4] = 30
    #writedlm(@sprintf("%s/data/_neurons_freq%d.dat", directory_name, idx_ntk), neurons_freq, header=false)
else
    #neurons_freq= readdlm(@sprintf("%s/data/_neurons_freq%d.dat",directory_name, idx_ntk))
end

if(SD=="yes")
    if(new_network==1)
        neurons_freq_noise = zeros(N_cycles, ncells)
        neurons_freq_noise[:,1].=1
        neurons_freq_noise[:,2:end] = rand([0.1 0.5 1], size(neurons_freq_noise[:,2:end]))
        #writedlm(@sprintf("%s/data/_noise_neurons_freq%d.dat", directory_name, idx_ntk), neurons_freq_noise, header=false)
    else
        #neurons_freq_noise= readdlm(@sprintf("%s/data/_noise_neurons_freq%d.dat", directory_name, idx_ntk))
    end
end

Iapp_cell = zeros(ncells, Tdt)

for idx_cycle= 1:1:N_cycles
    for idx=1:1:N_samples
        T1 = (idx_cycle-1)*Tdt_set +(idx-1)*Tdt_sample+1
        T2 = (idx_cycle-1)*Tdt_set + idx*Tdt_sample
        Iapp_cell[:, T1:T2] = get_Iapp(Duration_sample, dt, neurons_freq[(idx_cycle-1)*N_samples+idx,:], spike_duration)
    end
    if(SD=="yes")
        TA = (idx_cycle-1)*Tdt_set+Tdt_state+1
        TB = (idx_cycle-1)*Tdt_set+Tdt_cycle
        Iapp_cell[:,TA:TB] = get_Iapp(Duration_state, dt, neurons_freq_noise[idx_cycle,:], spike_duration)
    end
end

idxTesting1 = [1 2].+1
for idx_cycle=1:1:N_cycles
    T1 = (idx_cycle-1)*Tdt_set + Tdt_cycle +  Tdt_transient_testing +1
    T2 = (idx_cycle-1)*Tdt_set + Tdt_cycle + Tdt_transient_testing +Tdt_testing
    for idx=1:1:length(idxTesting1)
        Iapp_cell[idxTesting1[idx], T1:T2] = get_Iapp(Duration_testing, dt, [50], spike_duration)
    end
end
Iapp_cell[1:ncellsI,:] .= IappI


const  BurstTime = zeros(1,N_cycles)
for idx=1:1:N_cycles
     BurstTime[idx] = Duration_state + (idx-1)*Duration_set
end
const BurstDuration = Duration_state #20000



## synaptic plasticity
const expm = "Control"
const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410
const C_Post = 1.62138
const D_pre  = 9.53709 #[ms]

const tau_w    = 520.76129e3 #[s>ms]
const gamma_p  = 597.08922
const gamma_d  = 137.7586
const n_lin = 1.
const theta_p = 2.009289
const theta_d = 1.0

const Omega_p = gamma_p/(gamma_p+gamma_d)
const Omega_sleep = 0.75
const Omega_d = 0.
const Omega_0 = 0.

const tauw_p = tau_w /(gamma_p + gamma_d)
const tauw_d = tau_w / gamma_d
const tauw_0 = 0.
const zeta = tauw_d / tauw_p
const tau_x = 1


## AMPAFICATION
#const slope = 20
#const tauG = 5e5
#const tauINCR = 100

##   CONNECTIVITY

const gIGABAA = (2. * ones(ncellsC))./ ncellsI
const gIGABAB = (1.5 * ones(ncellsC))./ ncellsI


# Pre to post

idxPrePost = [CartesianIndex(1,ncells-2)]
for idx_q = 2:1:nPre
  push!(idxPrePost, CartesianIndex(idx_q, ncells-2))
end
for idx_q = 1:1:nPre
    push!(idxPrePost, CartesianIndex(idx_q, ncells-1))
end


#=
# to connect neuron manually
idxPrePost = [CartesianIndex(1, 2);
              CartesianIndex(3, 4)]
=#





gCAMPA = 0.1*ones(nPre,nPost)
w_init = 0.5*ones(nPre,nPost)


@time (V, CPre_plot, CPost_plot, Ca_plot, w, gCAMPA_plot) = simulateTOY_ncellsScenarioNMOD(
    ncells,
    ncellsI,
    ncellsC,
    Iapp_cell,
    Istep_cell,
    gCAMPA,
    idxPrePost
)


pw =plot(w[1:end,:]', legend= (frameon=false))
pg = plot(gCAMPA_plot[1:end,:]')
plot(pw,pg, layout=(2,1))



pCpre1 =plot(CPre_plot[1,:], legend= (frameon=false))
pCpre2 =plot(CPre_plot[2,:], legend= (frameon=false))
pCpost1 =plot(CPost_plot[1,:], legend= (frameon=false))
pCpost2 =plot(CPost_plot[2,:], legend= (frameon=false))
pCa1 =plot(Ca_plot[1,:], legend= (frameon=false))
pCa2 =plot(Ca_plot[2,:], legend= (frameon=false))

sum1 = CPre_plot[1,:]+ CPost_plot[1,:]
sum2 = CPre_plot[2,:]+ CPost_plot[2,:]
psum1 =plot(sum1, legend= (frameon=false))
psum2 =plot(sum2, legend= (frameon=false))


plot(pCpre1,pCpre2, layout=(2,1))
plot(pCpost1,pCpost2, layout=(2,1))
plot(pCa1,pCa2, layout=(2,1))
plot(psum1,psum2, layout=(2,1))



##
pV1 =plot(V[1,12000:18000], legend= (frameon=false))
pV2 =plot(V[2,:], legend= (frameon=false))
pV3 =plot(V[3,:], legend= (frameon=false))
pV4 =plot(V[4,:], legend= (frameon=false))
plot(pV1, pV2, pV3, pV4, layout=(4,1))


writedlm(@sprintf("%s/V.dat",directory_name), V[:,12000:18000], header=false)
writedlm(@sprintf("%s/CPre.dat",directory_name), CPre_plot[:,12000:18000], header=false)
writedlm(@sprintf("%s/CPost.dat",directory_name), CPost_plot[:,12000:18000], header=false)
writedlm(@sprintf("%s/w.dat",directory_name), w, header=false)
writedlm(@sprintf("%s/w_frame.dat",directory_name), w[:,12000:18000], header=false)