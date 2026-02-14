# Functions
function ISIfunc(Spkt::Array{Float64})   # Computes the interspike intervals out of spike times
  ISI = zeros(200000)
  l::Int64 = 1
  for i = 2:length(Spkt)
      ISI[l] = Spkt[i] - Spkt[i-1]
      l += 1
  end
  return ISI
end

function remove0(ISI::Array{Float64})    # Removes 0's in spike times/interspike intervals vectors
  f = findall(ISI .== 0)
    if f[1] > 1
      ISI = ISI[1:f[1]-1]
    else
      ISI = [0.]
    end
  return ISI
end

function isbursting(ISI::Array{Float64})   # Check if the neuron is bursting using the rule 3*minISI < maxISI
  bursting::Int64 = 0
  if minimum(ISI)*4 < maximum(ISI)
      bursting = 1
  end
  return bursting
end

function SPB_PER_DC_IBFfunc(ISI::Array{Float64}) # Computes spike per burst (SPB), period (PER), duty cycle (DC) and mean intraburst frequency of bursting neurons.
  minISI = minimum(ISI)
  maxISI = maximum(ISI)
  interburst = findall(ISI .>= maxISI/3)
  intraburst = findall(ISI .<= maxISI/3)
  SPB = round(length(intraburst)/length(interburst))+1
  IBP = mean(ISI[intraburst])
  Burstdur = IBP*(SPB-1)
  PER = Burstdur+mean(ISI[interburst])
  DC = Burstdur/PER
  IBF = 1000/IBP
  return (SPB, PER, DC, IBF)
end

function compute_params(spk_ncells::Array{Float64},Ttransient::Int64,tstepstart::Int64,tstepinit::Int64,tstepstop::Int64) # Extracts frequency or (SPB, PER, DC and IBF) out of the spiketimes saved in the Spkt$n.dat files. The number of cells is given by the number of lines in gs.
  ncells = length(spk_ncells[:,1])

  ISIs_depol = zeros(200000)
  ISIs_hyperpol = zeros(200000)
  PARAMS_depol = zeros(ncells,4) #[SPB;PER;DC;IBF]
  PARAMS_hyperpol = zeros(ncells,4) #[SPB;PER;DC;IBF]
  freq_depol = zeros(ncells)
  freq_hyperpol = zeros(ncells)

  for n = 1:ncells
    f1 = findall(spk_ncells[n,:] .>= Ttransient) #Removes the transient of depol period
    f2 = findall(spk_ncells[n,:] .>= tstepstart) #Finds end of depol period
    f3 = findall(spk_ncells[n,:] .>= tstepinit) #Removes the transient of hyperpol period
    f4 = findall(spk_ncells[n,:] .>= tstepstop) #Finds end of hyperpol period
    f5 = findall(spk_ncells[n,:] .== 0) #Removes the zeros at the end of the vector

    if length(f1) == 0
      Spkt_depol = [0.]
      Spkt_hyperpol = [0.]
    elseif length(f1) == length(f2)
      Spkt_depol = [0.]
      if length(f3) == length(f4)
        Spkt_hyperpol = [0.]
      elseif length(f4) > 0
        if(isempty(f3))
          Spkt_hyperpol = [0.]
        else
          Spkt_hyperpol = spk_ncells[n,f3[1]:f4[1]-1]
        end
      else
        if(isempty(f5)||isempty(f3))
          Spkt_hyperpol = [0.]
        else
          Spkt_hyperpol = spk_ncells[n,f3[1]:f5[1]-1]
        end
      end
    else
      if(isempty(f2))
        Spkt_depol = [0.]
      else
        Spkt_depol = spk_ncells[n,f1[1]:f2[1]-1]
      end
      if length(f3) == length(f4)
        Spkt_hyperpol = [0.]
      elseif length(f4) > 0
        if(isempty(f3))
          Spkt_hyperpol = [0.]
        else
          Spkt_hyperpol = spk_ncells[n,f3[1]:f4[1]-1]
        end
      else
        if(isempty(f5)||isempty(f3))
          Spkt_hyperpol = [0.]
        else
          Spkt_hyperpol = spk_ncells[n,f3[1]:f5[1]-1]
        end
      end
    end

    ISI_depol = ISIfunc(Spkt_depol)
    ISI_hyperpol = ISIfunc(Spkt_hyperpol)
    ISItemp_depol::Array{Float64} = remove0(ISI_depol)
    ISItemp_hyperpol::Array{Float64} = remove0(ISI_hyperpol)
    if isbursting(ISItemp_depol) == 1
      (PARAMS_depol[n,1],PARAMS_depol[n,2],PARAMS_depol[n,3],PARAMS_depol[n,4]) = SPB_PER_DC_IBFfunc(ISItemp_depol)
      #freq_depol[n] = 0.
    else
      PARAMS_depol[n,:] .= 0.
      freq_depol[n] = 1000/mean(ISItemp_depol)
    end

    if isbursting(ISItemp_hyperpol) == 1
      (PARAMS_hyperpol[n,1],PARAMS_hyperpol[n,2],PARAMS_hyperpol[n,3],PARAMS_hyperpol[n,4]) = SPB_PER_DC_IBFfunc(ISItemp_hyperpol)
      #freq_hyperpol[n] = 0.
    else
      PARAMS_hyperpol[n,:] .= 0.
      freq_hyperpol[n] = 1000/mean(ISItemp_hyperpol)
    end
    ISIs_depol = [ISIs_depol ISI_depol]
    ISIs_hyperpol = [ISIs_hyperpol ISI_hyperpol]
  end
  ISIs_depol = ISIs_depol[:,2:end]
  ISIs_hyperpol = ISIs_hyperpol[:,2:end]

  return ISIs_depol, ISIs_hyperpol, PARAMS_depol, PARAMS_hyperpol, freq_depol, freq_hyperpol
end
