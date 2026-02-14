Each Figure from the article is associated with a corresponding folder. 
Each folder is divided into the same structure:  
/data - contains the parameter used to generate the simulation 
/julia - main scripts to run the simulations.  
All computational experiments are launched thanks to a *Simu_... .jl* calling a function with the differential equations of the conductance-based models (neuronal model) and the synaptic weight change (plasticity) found in *model_... .jl*.  
/matlab - codes to analyze and plot the network simulations.  

Note: 
l_ij is indicated at g_AMPA_ij 
GB means Graupner and Brunel 's model

### Fig1
Simulate a network of neurons with / /Simu_scenario_GB2012.jl 
Inputs
gion.dat: intrinsic parameters of the neurons 
gsyn.dat: initial parameters of the late weights 
neurons_freq.dat: frequencies associated with each neuron at the different states.

Outputs
- Voltage recordings
- Firing pattern properties
- LFP recordings
- Synaptic weight recordings

MATLAB: Plot traces and quantify firing pattern properties with /matlab/Fig1.m 

### Fig2_Demo
Simulate a tonic or a burst phase for the model Graupner2016.
Figure2B is made using the code Simu_reduced_AMPA.jl in CalciumRule 

Inputs
gion.dat: intrinsic parameters of the neurons 
gsyn.dat: initial parameters of the synaptic weights 
neurons_freq.dat: frequencies associated with each neuron at the different states.

Outputs
- synaptic weight recordings
- histogram of alpha_d, alpha_p, alpha_0
- prediction vs simulation


### Fig3 
Simu_TUNE.jl is associated with Fig SI for one tonic firing state followed by three burst firing states with different hyperpolarizing currents.  
 
MATLAB 
Fig3.m: code associated with SIMU_TUNE.jl 

### Fig3B
Simu_TUNE_sensitivity.jl is associated with Fig3B. 
 

### Fig4
Simu_NMOD.jl associated with Fig 4B
Simu_NMOD_TAG.jl is associated with Fig4D.

MATLAB 
Fig4_scenario_NMOD.m: plots the evolution of the weight as function of time. Change the expm_name for plotting NMOD or TAG. 
 
  

###  FigS2
Simu_XX.jl runs the different simulations for the model  XX     (can be Graupner, Deperrois,  etc).  
 
MATLAB 
Fig3_scatter_SB.m retrieves the scatter plot for soft-bound models 
Fig3_scatter_HB.m retrieves the scatter plot for hard-bound models 
 
 

}

