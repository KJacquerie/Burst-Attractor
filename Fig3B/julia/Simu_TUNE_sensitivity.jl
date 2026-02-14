# Defines output directory
directory_name = "/Users/kathleen/Documents/PhD/2025-Paper1/Fig3C"

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

include("model_TUNE_sensitivity.jl")
include("PARAMS_cycle.jl")
const new_network=0
const fMax = 50

const bound_type = "SB"
const SD = "no"
if(SD=="yes")
    experiment_name = "Graupner2012_TUNE_SD"
else
    experiment_name = "Graupner2012_TUNE"
end
idx_ntk = 1

# Network parameters
const ncellsI = 1
const nPre = 50
const nPost = 50
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC


## Simulation parameters
const dt = 0.01

const N_states= 10

const N_cycles=1
const Duration_cycle = 20000*N_states
const Tdt_cycle = convert(Int64, Duration_cycle/dt)


const Duration_state = convert(Int64, Duration_cycle/N_states)
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
#=
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
=#
state[1:Tdt_state] .=0
state[Tdt_state+1:end] .=1


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



if(new_network==1 && fMax==70)
    const gamma = 0.20 #1  #10% of variability in the network put 0.1 to get 10% etc
    const gl_cells = rand(Uniform(gl*(1-gamma),gl*(1+gamma)),ncells)
    const gNa_cells = rand(Uniform(gNa*(1-gamma),gNa*(1+gamma)),ncells)
    const gKd_cells = rand(Uniform(gKd*(1-gamma),gKd*(1+gamma)),ncells)
    const k1_cells = rand(Uniform(k1*(1-gamma),k1*(1+gamma)),ncells) #k1*ones(ncells) 
    const k2_cells = rand(Uniform(k2*(1-gamma),k2*(1+gamma)),ncells) #k2*ones(ncells) 
    const gH_cells = rand(Uniform(gH*(1-gamma),gH+(gamma*gH)),ncells)
    const gKCa_cells = rand(Uniform(gKCa*(1-gamma),gKCa*(1+gamma)),ncells)
    const gCaT_cells = rand(Uniform(gCaT*(1-gamma),gCaT*(1+gamma)),ncells)
    const g_cond = [gNa_cells'; gKd_cells'; gCaT_cells'; gH_cells'; gKCa_cells'; gl_cells'; k1_cells'; k2_cells']'
    writedlm(@sprintf("%s/data/%s/gion.dat",directory_name, experiment_name), g_cond, header=false)
else 
    gion = readdlm(@sprintf("%s/data/%s/gion.dat",directory_name, experiment_name))
    const gl_cells = gion[:,6]
    const gNa_cells = gion[:,1]
    const gKd_cells = gion[:,2]
    const k1_cells = gion[:,end-1]
    const k2_cells = gion[:,end]
    const gH_cells = gion[:,4]
    const gKCa_cells = gion[:,5]
    const gCaT_cells = gion[:,3]
end



const IappI = 3.
const IappC = 0.

const spike_duration = 3 


const IstepI_sens = -1.7-IappI:0.10:-0.9-IappI

##

const IstepC = 50.0 

#IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
#Istep_cell = [IstepI_cell; IstepC_cell]


if(new_network==1)
    neurons_freq  = zeros(N_samples*N_cycles, ncells)   
    neurons_freq[:,1] .= 1 # inhibitory cell

    neurons_freq[:,2:end] = round.(rand(Uniform(0.1,fMax),N_samples*N_cycles, ncells-1))
    writedlm(@sprintf("%s/data/%s/neurons_freq.dat", directory_name, experiment_name), neurons_freq, header=false)
else
    neurons_freq= readdlm(@sprintf("%s/data/%s/neurons_freq.dat",directory_name, experiment_name))
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
Iapp_cell[1:ncellsI,:] .= IappI


const  BurstTime = zeros(1,N_states)
const StateTime = zeros(1,N_states)
for idx_state=1:1:N_states
     BurstTime[idx_state] = Duration_state*idx_state
     StateTime[idx_state]   = idx_state*Duration_state
end
const BurstDuration = Duration_state #20000


## synaptic plasticity
const expm = "Control"
const tau_Ca = 22.6936 #[ms]
const C_Pre  = 0.56#17539
const C_Post = 1.24#23964
const D_pre  = 4.60#98 #[ms]

const tau_w    = 346.3615e3 #[s>ms]
const gamma_p  = 725.085*1.1
const gamma_p_burst  = 725.085*0.95
const gamma_d  = 331.909
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5


##   CONNECTIVITY

const gIGABAA_unit = 2.0
const gIGABAB_unit = 1.5
if(new_network==1 && fMax==70)
    const gIGABAA = rand(Uniform(gIGABAA_unit*(1-gamma),gIGABAA_unit*(1+gamma)),ncellsC)./ ncellsI
    const gIGABAB = rand(Uniform(gIGABAB_unit*(1-gamma),gIGABAB_unit*(1+gamma)),ncellsC)./ ncellsI
    const g_syn = [gIGABAA'; gIGABAB']'
    writedlm(@sprintf("%s/data/%s/gsyn.dat",directory_name, experiment_name), g_syn, header=false)
else
    gsyn = readdlm(@sprintf("%s/data/%s/gsyn.dat",directory_name, experiment_name))
    const gIGABAA = gsyn[:,1]
    const gIGABAB = gsyn[:,2]
end

gCAMPA = 0.001*ones(nPre,nPost)
w_init = 0.5*ones(nPre,nPost)

idx_wl = convert(Matrix{Int64}, readdlm("idx_wl.dat"))

@time (Vspk,LFP_C, Vbin) = simulateTOY_ncellsScenarioNMOD(
    ncells,
    ncellsI,
    ncellsC,
    Iapp_cell,
    gCAMPA
)



SPB_depol1 = zeros(ncells,N_states)
SPB_hyperpol1 = zeros(ncells,N_states)
PER_depol1 = zeros(ncells,N_states)
PER_hyperpol1 = zeros(ncells,N_states)
DC_depol1 = zeros(ncells,N_states)
DC_hyperpol1 = zeros(ncells,N_states)
IBF_depol1 = zeros(ncells,N_states)
IBF_hyperpol1 = zeros(ncells,N_states)
freq_depol1 = zeros(ncells,N_states)
freq_hyperpol1 = zeros(ncells,N_states)

for idx_state=1:1:N_states
    Duration_CompParam= 1000
    T0 = Duration_CompParam
    TA = Duration_state
    TB = idx_state*Duration_state
    TC = (idx_state+1)*Duration_state
    (ISIs_depol, ISIs_hyperpol, PARAMS_depol, PARAMS_hyperpol, freq_depol, freq_hyperpol) = compute_params(Vspk[:,:],T0,TA,TB,TC)
    SPB_depol1[:,idx_state] = PARAMS_depol[:,1]
    SPB_hyperpol1[:,idx_state] = PARAMS_hyperpol[:,1]
    PER_depol1[:,idx_state] = PARAMS_depol[:,2]
    PER_hyperpol1[:,idx_state] = PARAMS_hyperpol[:,2]
    DC_depol1[:,idx_state] = PARAMS_depol[:,3]
    DC_hyperpol1[:,idx_state] = PARAMS_hyperpol[:,3]
    IBF_depol1[:,idx_state] = PARAMS_depol[:,4]
    IBF_hyperpol1[:,idx_state] = PARAMS_hyperpol[:,4]
    freq_depol1[:,idx_state] = freq_depol
    freq_hyperpol1[:,idx_state] = freq_hyperpol
end


const fcutoff = 100
responsetype = Lowpass(fcutoff; fs=1000/(dt))
designmethod = Butterworth(4)
LFP_C_filt = filt(digitalfilter(responsetype, designmethod), LFP_C)

writedlm(@sprintf("%s/data/%s/Vspk.dat",directory_name, experiment_name), Vspk, header=false)
writedlm(@sprintf("%s/data/%s/Vbin.dat",directory_name, experiment_name), Vbin, header=false)
#writedlm(@sprintf("%s/data/%s/w.dat",directory_name, experiment_name), w, header=false)
writedlm(@sprintf("%s/data/%s/LFP.dat",directory_name, experiment_name), LFP_C_filt, header=false)

writedlm(@sprintf("%s/data/%s/SPB_tonic.dat",directory_name, experiment_name), SPB_depol1, header=false)
writedlm(@sprintf("%s/data/%s/SPB_burst.dat",directory_name, experiment_name), SPB_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/PER_tonic.dat",directory_name, experiment_name), PER_depol1, header=false)
writedlm(@sprintf("%s/data/%s/PER_burst.dat",directory_name, experiment_name), PER_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/DC_tonic.dat",directory_name, experiment_name), DC_depol1, header=false)
writedlm(@sprintf("%s/data/%s/DC_burst.dat",directory_name, experiment_name), DC_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/IBF_tonic.dat",directory_name, experiment_name), IBF_depol1, header=false)
writedlm(@sprintf("%s/data/%s/IBF_burst.dat",directory_name, experiment_name), IBF_hyperpol1, header=false)
writedlm(@sprintf("%s/data/%s/freq_tonic.dat",directory_name, experiment_name), freq_depol1, header=false)
writedlm(@sprintf("%s/data/%s/freq_burst.dat",directory_name, experiment_name), freq_hyperpol1, header=false)




