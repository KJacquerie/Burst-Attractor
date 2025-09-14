# Defines output directory
directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig4_Model"

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
using DSP

## Include model

include("model_PairBased.jl")
include("PARAMS_cycle.jl")
const new_network=0

const bound_type = "HB"
const region = "HPC"
const SD = "no"
const experiment_name = "PairBased"
const wMAX = 1
const wMIN = 0

idx_ntk = 1

# Network parameters
const ncellsI = 1
const nPre = 50
const nPost = 50
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC


## Simulation parameters
const dt = 0.01

const N_cycles=2
const Duration_cycle = 40000
const Tdt_cycle = convert(Int64, Duration_cycle/dt)

const N_patterns = 1
const Duration_state = convert(Int64, Duration_cycle/2)
const Tdt_state = convert(Int64, Duration_state/dt)

const N_samples = 1
const Duration_sample = convert(Int64, Duration_state/N_samples)
const Tdt_sample = convert(Int64, Duration_sample / dt)


const T = N_cycles*Duration_cycle 
const Tdt = convert(Int64, T / dt)
const t = range(dt, T, length = Tdt)

const Duration_set = Duration_cycle 
const Tdt_set = convert(Int64, Duration_set/dt)

const state = zeros(Tdt,1)
for idx_cycle=1:1:N_cycles
    # wake
    state[(idx_cycle-1)*Tdt_set+1: (idx_cycle-1)*Tdt_set+Tdt_state] .= 0
    # burst
    if(SD=="yes")
        state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= -1
    else
        state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= 1
    end
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




    gion = readdlm(@sprintf("%s/data/gion.dat",directory_name))
    const gl_cells = gion[:,6]
    const gNa_cells = gion[:,1]
    const gKd_cells = gion[:,2]
    const k1_cells = gion[:,end-1]
    const k2_cells = gion[:,end]
    const gH_cells = gion[:,4]
    const gKCa_cells = gion[:,5]
    const gCaT_cells = gion[:,3]




const IappI = 3.
const IappC = 0.

const spike_duration = 3 
const IstepI = -1.2-IappI
const IstepC = 50.0 

IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
Istep_cell = [IstepI_cell; IstepC_cell]


if(new_network==1)
    neurons_freq  = zeros(N_samples*N_cycles, ncells)   
    neurons_freq[:,1] .= 1 # inhibitory cell

    neurons_freq[:,2:end] = round.(rand(Uniform(0.1,50),N_samples*N_cycles, ncells-1))
    writedlm(@sprintf("%s/data/neurons_freq_PB.dat", directory_name), neurons_freq, header=false)
else
    neurons_freq= readdlm(@sprintf("%s/data/neurons_freq_PB.dat",directory_name))
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
Iapp_cell[1:ncellsI,:] .= IappI
LL = copy(Iapp_cell)
for idw=1:1:nPost
    Duration_space = rand(-10:10) #max 100ms 
    Tdt_space = convert(Int64, Duration_space/dt)
    LL[nPre+1+idw,1:1:Tdt_state] = circshift(Iapp_cell[1+idw,1:1:Tdt_state], Tdt_space)
end

for idw=1:1:nPost
    Duration_space = rand(-10:10) #max 100ms 
    Tdt_space = convert(Int64, Duration_space/dt)
    LL[nPre+1+idw,Tdt_cycle:1:Tdt_cycle+Tdt_state] = circshift(Iapp_cell[1+idw,Tdt_cycle:1:Tdt_cycle+Tdt_state], Tdt_space)
end
Iapp_cell = copy(LL)
Iapp_cell[1:ncellsI,:] .= IappI

const  BurstTime = zeros(1,N_cycles)
for idx=1:1:N_cycles
     BurstTime[idx] = Duration_state + (idx-1)*Duration_set
end
const BurstDuration = Duration_state #20000

const StateTime = zeros(1,N_cycles*2)
let idx_count
    idx_count=1
    for idx_cycle=1:1:N_cycles
        StateTime[idx_count]   = (idx_cycle-1)*Duration_set+Duration_state
        StateTime[idx_count+1] = (idx_cycle-1)*Duration_set+Duration_cycle
        idx_count = idx_count+2
    end
end

## synaptic plasticity
const expm = "Control"

const Am = 0.0053
const Ap = 0.0096
const tau_m = 33.7
const tau_p = 16.8

const wMax = 1.
const wMin = 0.
    
if(bound_type=="SB")
    const SB = 1
else
    const SB = 0
end



##   CONNECTIVITY

const gIGABAA_unit = 2.0
const gIGABAB_unit = 1.5

    gsyn = readdlm(@sprintf("%s/data/gsyn.dat",directory_name))
    const gIGABAA = gsyn[:,1]
    const gIGABAB = gsyn[:,2]


gCAMPA = 0.001*ones(nPre,nPost)
w_init = 0.5*ones(nPre,nPost)

@time (Vspk, Vspk_bin,w) = simulateTOY_ncellsScenarioNMOD(
    ncells,
    ncellsI,
    ncellsC,
    Iapp_cell,
    Istep_cell,
    gCAMPA
)



SPB_depol1 = zeros(ncells,N_cycles)
SPB_hyperpol1 = zeros(ncells,N_cycles)
PER_depol1 = zeros(ncells,N_cycles)
PER_hyperpol1 = zeros(ncells,N_cycles)
DC_depol1 = zeros(ncells,N_cycles)
DC_hyperpol1 = zeros(ncells,N_cycles)
IBF_depol1 = zeros(ncells,N_cycles)
IBF_hyperpol1 = zeros(ncells,N_cycles)
freq_depol1 = zeros(ncells,N_cycles)
freq_hyperpol1 = zeros(ncells,N_cycles)

for idx_cycle=1:1:N_cycles
    Duration_CompParam= 1000
    T0 = (idx_cycle-1)*Duration_set + Duration_CompParam
    println("T0",T0)
    TA = (idx_cycle-1)*Duration_set+Duration_state+1
    println("TA",TA)
    TB = (idx_cycle-1)*Duration_set+Duration_cycle
    println("TB",TB)
    (ISIs_depol, ISIs_hyperpol, PARAMS_depol, PARAMS_hyperpol, freq_depol, freq_hyperpol) = compute_params(Vspk[:,:],T0,TA,TA,TB)
    SPB_depol1[:,idx_cycle] = PARAMS_depol[:,1]
    SPB_hyperpol1[:,idx_cycle] = PARAMS_hyperpol[:,1]
    PER_depol1[:,idx_cycle] = PARAMS_depol[:,2]
    PER_hyperpol1[:,idx_cycle] = PARAMS_hyperpol[:,2]
    DC_depol1[:,idx_cycle] = PARAMS_depol[:,3]
    DC_hyperpol1[:,idx_cycle] = PARAMS_hyperpol[:,3]
    IBF_depol1[:,idx_cycle] = PARAMS_depol[:,4]
    IBF_hyperpol1[:,idx_cycle] = PARAMS_hyperpol[:,4]
    freq_depol1[:,idx_cycle] = freq_depol
    freq_hyperpol1[:,idx_cycle] = freq_hyperpol
end


writedlm(@sprintf("%s/data/%s/%s/%s/Vspk.dat",directory_name, experiment_name, region, bound_type), Vspk, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/Vbin.dat",directory_name, experiment_name, region, bound_type), Vspk_bin, header=false)
#if(bound_type=="HB")
    writedlm(@sprintf("%s/data/%s/%s/%s/w.dat",directory_name, experiment_name, region, bound_type), w, header=false)
#end
writedlm(@sprintf("%s/data/%s/%s/%s/SPB_tonic.dat",directory_name, experiment_name, region, bound_type), SPB_depol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/SPB_burst.dat",directory_name, experiment_name, region, bound_type), SPB_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/PER_tonic.dat",directory_name, experiment_name, region, bound_type), PER_depol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/PER_burst.dat",directory_name, experiment_name, region, bound_type), PER_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/DC_tonic.dat",directory_name, experiment_name, region, bound_type), DC_depol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/DC_burst.dat",directory_name, experiment_name, region, bound_type), DC_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/IBF_tonic.dat",directory_name, experiment_name, region, bound_type), IBF_depol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/IBF_burst.dat",directory_name, experiment_name, region, bound_type), IBF_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/freq_tonic.dat",directory_name, experiment_name, region, bound_type), freq_depol1, header=false)
writedlm(@sprintf("%s/data/%s/%s/%s/freq_burst.dat",directory_name, experiment_name, region, bound_type), freq_hyperpol1, header=false)





