# Defines output directory

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

include("model_scenario_NMOD.jl")

## Simulation parameters
const N_cycles=2
const Duration_cycle = 30000
const Duration_state = convert(Int64, Duration_cycle/2)
const T = N_cycles*Duration_cycle
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Ttransient = convert(Int64, 1 / dt)
const t = range(dt, T, length = Tdt)

## Neurons' model parameters

# Global parameters
const C = 1
const VNa = 50
const VK = -85
const VCa = 120
const Vl = -55
const VH = -20
const Kd = 170

# Network parameters
const ncellsI = 1
const ncellsC = 12
const ncells = ncellsI+ncellsC

# /!\ For all vectors, the first rows correspond to the inhibitory cells (ncellsI)

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


#possible_freq =  [ 0.1, 0.5, 1, 5,  10, 15, 20, 25, 30, 35, 40]#, 45, 50, 55, 60]
#rand(possible_freq, ncells)# Neurons spiking frequencies in Hz

neurons_freq = [1  30 25 28   10 5 10   30 60 40  10 5 1;
                1  35 50 60   0.5 10 1  47 28 45  0.1 10 3]
                #spike_duration = 3 # spike duration in ms
Iapp_cell = zeros(ncells, Tdt)
Tdt_cycle = convert(Int64, Duration_cycle / dt)
for idx=1:1:N_cycles
    Iapp_cell[:, (idx-1)*Tdt_cycle+1:idx*Tdt_cycle] = get_Iapp(Duration_cycle, dt, neurons_freq[idx,:], spike_duration)
    #Iapp_cell[:, (idx-1)*Tdt_cycle+1:idx*Tdt_cycle] = get_Iapp(Duration_cycle, dt, neurons_freq[idx,:], spike_duration)
end
Iapp_cell[1:ncellsI,:] .= IappI

const  BurstTime = zeros(1,N_cycles)
for idx=1:1:N_cycles
     BurstTime[idx] = Duration_state + (idx-1)*Duration_cycle
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
const Omega_p_sleep = 0.65#gamma_p/(gamma_p+gamma_d)
const Omega_d = 0.
const Omega_0 = 0.

const tauw_p = tau_w /(gamma_p + gamma_d)
const tauw_d = tau_w / gamma_d
const tauw_0 = 0.
const zeta = tauw_d / tauw_p
const tau_x = 1

##   CONNECTIVITY

const gIGABAA = (2. * ones(ncellsC))./ ncellsI
const gIGABAB = (1.5 * ones(ncellsC))./ ncellsI


# Pre to post
idxPrePost = [CartesianIndex(1, 7);
              CartesianIndex(2, 8);
              CartesianIndex(3, 9);
              CartesianIndex(4, 10);
              CartesianIndex(5, 11);
              CartesianIndex(6, 12)]

IstepPrePost_cell = zeros(ncellsC,ncellsC)
IstepPrePost_cell[idxPrePost] .= 50.

gCAMPA = zeros(ncellsC,ncellsC)
gCAMPA[idxPrePost] .= 0.1
const w_init = zeros(ncellsC,ncellsC)
w_init[idxPrePost] .= 0.5



@time (w,V) = simulateTOY_ncellsScenarioNMOD(
    ncells,
    ncellsI,
    ncellsC,
    Iapp_cell,
    Istep_cell,
    gCAMPA,
    idxPrePost
)

##

T1 = convert(Int64, 1)
T2 = convert(Int64, T)


## W ensemble
colorss = ["blue","red","cyan","black","gold","green","orchid","teal","indigo"]
#const idx_color = 1

plot()
#idx_color=copy(1)
for idx=1:1:size(w,1)
    if(CartesianIndex(rowColFromIdx(idx, ncellsC)) in idxPrePost)
        plot!(
            t[T1:T2],
            w[idx, T1:T2],
            #color = color_choice,
            ylabel="w",
            legend= (frameon=false)
        )

    end

end
p1 = plot!()





writedlm("/Users/kathleen/Documents/PhD/2022-reset/NMOD-SP/scenario/data/w_NMOD.dat", w, header = false)
#writedlm("/Users/kathleen/Documents/PhD/2022-reset/NMOD-SP/scenario/data/V.dat", V, header = false)
