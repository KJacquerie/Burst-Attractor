boltz(V::Float64,A::Float64,B::Float64) = 1/(1 + exp((V+A)/B))
tauX(V::Float64,A::Float64,B::Float64,D::Float64,E::Float64) = A - B/(1+exp((V+D)/E))
mNainf(V::Float64) = boltz(V,35.5,-5.29)
taumNa(V::Float64) = tauX(V,1.32,1.26,120.,-25.)
hNainf(V::Float64) = boltz(V,48.9,5.18)
tauhNa(V::Float64) = (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)))
mKdinf(V::Float64) = boltz(V,12.3,-11.8)
taumKd(V::Float64) = tauX(V,7.2,6.4,28.3,-19.2)
mCaTinf(V::Float64) = boltz(V,67.1,-7.2)
taumCaT(V::Float64) = tauX(V,21.7,21.3,68.1,-20.5)
hCaTinf(V::Float64) = boltz(V,80.1,5.5)
tauhCaT(V::Float64) = 2*tauX(V,205.,89.8,55.,-16.9)
mHinf(V::Float64) = boltz(V,80.,6.)
taumH(V::Float64) = tauX(V,272.,-1149.,42.2,-8.73)
mKCainf(Ca::Float64) = (Ca/(Ca+Kd))^2
Tm(V::Float64) = 1/(1+exp(-(V-2)/5))

ICaT(V::Vector{Float64},mCaT::Vector{Float64},hCaT::Vector{Float64},gCaT::Vector{Float64}) = gCaT .* mCaT.^3 .* hCaT .* (V .- VCa)
ICaT(V::Float64,mCaT::Float64,hCaT::Float64,gCaT::Float64) = gCaT * mCaT^3 * hCaT * (V - VCa)

function dV(gNa::Float64, gKd::Float64, gl::Float64, V::Float64, mNa::Float64, hNa::Float64, mKd::Float64, mCaT::Float64, hCaT::Float64, mH::Float64, Ca::Float64, Iapp::Float64, Istep::Float64, gCaT::Float64, gH::Float64, gKCa::Float64)
  (dt)*(1/C)*(-gNa*mNa^3*hNa*(V-VNa) -gKd*mKd^4*(V-VK) -gl*(V-Vl) -ICaT(V,mCaT,hCaT,gCaT) -gKCa*mKCainf(Ca)*(V-VK) -gH*mH*(V-VH) +Iapp +Istep)
end


dmNa(V::Float64,mNa::Float64) = (dt)*((1/taumNa(V))*(mNainf(V) - mNa))
dhNa(V::Float64,hNa::Float64) = (dt)*((1/tauhNa(V))*(hNainf(V) - hNa))
dmKd(V::Float64,mKd::Float64) = (dt)*((1/taumKd(V))*(mKdinf(V) - mKd))
dmCaT(V::Float64,mCaT::Float64) = (dt)*((1/taumCaT(V))*(mCaTinf(V) - mCaT))
dhCaT(V::Float64,hCaT::Float64) = (dt)*((1/tauhCaT(V))*(hCaTinf(V) - hCaT))
dmH(V::Float64,mH::Float64) = (dt)*((1/taumH(V))*(mHinf(V) - mH))
dCa(V::Float64,mCaT::Float64,hCaT::Float64,Ca::Float64,gCaT::Float64,k1::Float64,k2::Float64) = (dt)*(-k1*ICaT(V,mCaT,hCaT,gCaT) -k2*Ca)
dAMPA(V::Float64,AMPA::Float64) = (dt)*(1.1*Tm(V)*(1-AMPA)-0.19*AMPA)
dGABAA(V::Float64,GABAA::Float64) = (dt)*(0.53*Tm(V)*(1-GABAA)-0.18*GABAA)
dGABAB(V::Float64,GABAB::Float64) = (dt)*(0.016*Tm(V)*(1-GABAB)-0.0047*GABAB)



function rowColFromIdx(index,max_column)
    i = ((index-1)%max_column)+1
    j = trunc(Int64,((index-1)/max_column))+1
    return i,j
end

function IdxFromRowCol(row,col,max_column)
    idx = col*max_column + row - max_column
    return idx
end

function dxSpike(D_pre::Float64)
   return D_pre
 end

function dx(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_p + 1*(Vprev<= 0 && V >0)
end

function dx2(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_x + 1*(Vprev<= 0 && V >0)
end

function dy(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_m + 1*(Vprev<= 0 && V >0)
end

function dy2(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_y + 1*(Vprev<= 0 && V >0)
end

function dw(new::Float64, w::Float64, tau::Float64)
  (dt)*(new - w)*1/tau
end

function get_Iapp(T, dt, neurons_freq, spike_duration)

    nb_neurons = size(neurons_freq)[1]
    nb_steps = convert(Int64, round(T/dt))
    Iapp = zeros(nb_neurons, nb_steps)

     # Convert spike duration from ms to index
    spike_duration = convert(Int64, round(spike_duration/dt))

     # Generate Iapp
    for (neuron_idx, freq) in enumerate(neurons_freq)
      # Compute spike period in index
      if(freq>0)
          spike_period = convert(Int64, round(1000/(freq*dt)))
      
          if(spike_period<= convert(Int64, round(T/dt)))
            # Compute indices at which neurons spikes
            spike_indices = collect(rand(1:spike_period):spike_period:nb_steps)
            spike_std = 0.1 * spike_period
            spike_indices_err =   convert(
                                        Array{Int64,1},
                                        broadcast(round, spike_std * randn(
                                            size(spike_indices))))
            spike_indices += spike_indices_err
            spike_indices = broadcast(abs, spike_indices)
            # Check if first index is not 0
            if spike_indices[1] == 0
                spike_indices[1] += 1
            end
          end
      end
      

        # Set Iapp to 50 when spiking
      if(freq>0)
          for spike_idx in spike_indices
            up_val = min(nb_steps, spike_idx + spike_duration - 1)
            interval = spike_idx:up_val

            Iapp[neuron_idx, interval] = IstepC * ones(size(interval))
          end
      end

           
    end


    return Iapp
end



###
function simulateTOY_ncellsScenarioNMOD(ncells::Int64,ncellsI::Int64,ncellsC::Int64,Iapp_cell,Istep_cell::Vector{Float64},gCAMPA::Matrix{Float64})
  # Initial conditions
  V=-60*ones(ncells)
  Vprev=-60*ones(ncells)

  mNa=mNainf(V[1])*ones(ncells)
  hNa=hNainf(V[1])*ones(ncells)
  mKd=mKdinf(V[1])*ones(ncells)
  mCaT=mCaTinf(V[1])*ones(ncells)
  hCaT=hCaTinf(V[1])*ones(ncells)
  mH=mHinf(V[1])*ones(ncells)
  Ca= (-k1_cells ./ k2_cells) .* ICaT(V,mCaT,hCaT,gCaT_cells)

  AMPA=zeros(nPre)
  GABAA = zeros(ncellsI)
  GABAB = zeros(ncellsI)
  GABAA_prev = zeros(ncellsI)
  GABAB_prev = zeros(ncellsI)

  ww_var = ones(nPre, nPost)*0.5
  ww=zeros(length(ww_var),convert(Int64,T))
  gCAMPA_var = gCAMPA
  w_BACK = zeros(length(ww_var), N_cycles*2)

  x = zeros(nPre)
  x2 = zeros(nPre)
  x_prev = zeros(nPre)
  x2_prev = zeros(nPre)
  y = zeros(nPost)
  y2 = zeros(nPost)
  y_prev = zeros(nPost)
  y2_prev = zeros(nPost)
  
  pre_spikes = zeros(nPre)
  post_spikes = zeros(nPost)

  idx_w = 1
  Spkt = zeros(ncells,T)
  Spkt_binary = zeros(ncells,T)
  lncells::Array{Int64} = ones(ncells,1)

  que = []
  for idx_q = 1:1:ncellsC
    push!(que, Queue{Int}())
  end

  # Burst start and stop values
  TstartBurstTime = convert(Array{Int64,2}, BurstTime./dt)
  TBurstDuration::Int64 = convert(Int64,BurstDuration/dt)
  Tdt_StateTime = convert(Array{Int64,2}, StateTime./dt)
  idx_burst = 1

  idx_set= 1
  println("cycle=",idx_set)
  idx_state=1


  for z = 1:Tdt
    Isyn=zeros(ncells)

    if (z==TstartBurstTime[idx_burst]+TBurstDuration+1)
      idx_burst+=1
      if(idx_burst>N_cycles) # Testing phase !
        idx_burst=N_cycles
      end
    end

    if(z==Tdt_StateTime[idx_state])
      w_BACK[:,idx_state] = copy(reshape(ww_var, length(ww_var), 1))
      #g_BACK[:,idx_state] = copy(reshape(gCAMPA_var, length(gCAMPA_var), 1)) 
    
      writedlm(@sprintf("%s/data/%s/%s/%s/w_BACK.dat",directory_name, experiment_name, region, bound_type), w_BACK, header=false)
      #writedlm(@sprintf("%s/data/task_%s/%s/g_BACK.dat",directory_name), g_BACK, header=false)
      
      idx_state = idx_state+1
      if(idx_state >=N_cycles*2)
        idx_state = N_cycles*2
      end
      println("state=",idx_state)
      println("time=", t[z])
    end


    # ---- INHIBITORY CELLS -----
    for j = 1:ncellsI
      Iapp = Iapp_cell[j,z]

      if state[z] == 1 # burst
        Iappstep = Istep_cell[j]
      else
        Iappstep = 0.
      end

      V[j] += dV(gNa_cells[j], gKd_cells[j], gl_cells[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaT_cells[j], gH_cells[j], gKCa_cells[j])
      Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaT_cells[j],k1_cells[j], k2_cells[j])

      if(V[j]>=0. && Vprev[j]<0.)
        Spkt[j,lncells[j]] = t[z]
        Spkt_binary[j,idx_w] = 1
        lncells[j] +=1
      end

      mNa[j]  += dmNa(Vprev[j],mNa[j])
      hNa[j]  += dhNa(Vprev[j],hNa[j])
      mKd[j]  += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j]   += dmH(Vprev[j],mH[j])

      GABAA_prev[j] = GABAA[j]
      GABAB_prev[j] = GABAB[j]

      GABAA[j] += dGABAA(Vprev[j],GABAA[j])
      GABAB[j] += dGABAB(Vprev[j],GABAB[j])
      Vprev[j] = copy(V[j])



    end

    # ---- CORTICAL CELLS ------
    for j=ncellsI+1:ncellsC+ncellsI
        jj  = j-ncellsI

        Iapp = 0.
        Iappstep = 0.

        if state[z] == 1 # burst
            Iappstep = rand()
        else  # wake / test
            Iappstep = Iapp_cell[j,z]#+10*rand()
        end

        V[j] += dV(gNa_cells[j], gKd_cells[j], gl_cells[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaT_cells[j], gH_cells[j], gKCa_cells[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaT_cells[j],k1_cells[j], k2_cells[j])

        # Inhibitory synaptic currents from the I cell
        for k=1:1:ncellsI
            V[j] += (dt)*(1/C)*(-gIGABAA[jj]*GABAA_prev[k]*(Vprev[j]+70))
            V[j] += (dt)*(1/C)*(-gIGABAB[jj]*GABAB_prev[k]*(Vprev[j]+85))
        end

        if(j>ncellsI+nPre) # for all postsynaptic cells
          aa=j-ncellsI-nPre
          for k=1:1:nPre # check each presynaptic neuron
            wnorm = gCAMPA_var[k,aa]*ww_var[k,aa]
            V[j] += (dt)*(1/C)*(-wnorm*AMPA[k]*(Vprev[j]-0))
          end
        end


        # Check the state of the PRESYNAPTIC NEURONS
        if(j<=ncellsI+nPre)
          aa = j-ncellsI

          pre_spikes[aa] = 0
          x_prev[aa] = x[aa]
          x[aa] += dx(x[aa], V[j], Vprev[j])

          x2_prev[aa] = x2[aa]
          x2[aa] += dx2(x2[aa], V[j], Vprev[j])

          if(V[j]>=0. && Vprev[j]<0.)
            pre_spikes[aa] = 1 # Pré en dernier
          end
          AMPA[aa] += dAMPA(Vprev[j], AMPA[aa])
          
        # Check the state of  the POSTSYNAPTIC NEURONS
        else
          aa = j-ncellsI-nPre
          post_spikes[aa] = 0

          y_prev[aa] = y[aa]
          y[aa] += dy(y[aa], V[j], Vprev[j])

          y2_prev[aa] = y2[aa]
          y2[aa] += dy2(y2[aa], V[j], Vprev[j])

          if (V[j] >= 0 && Vprev[j]< 0)
            post_spikes[aa] = 1 #Post en dernier
          end
        end

        # We check if the presynaptic neuron has spiked recently
        # to update the weight in depression 
        if(j>ncellsI+nPre)
          ll = j-ncellsI-nPre # (ll) POST 
          for k=1:1:nPre      # (k)  PRE
            if pre_spikes[k] ==1
              # -DEPRESSION- 
              # the presynaptic neuron has spiked last
              if SB == 0 
                new = -y_prev[ll]*A2_m  + ww_var[k,ll]
              else
                new = -y_prev[ll]*A2_m* ww_var[k,ll] + ww_var[k,ll]
              end
              ww_var[k,ll] = new
              if(ww_var[k,ll]> wMax)
                ww_var[k,ll] = wMax
              end
              if(ww_var[k,ll]< wMin)
                ww_var[k,ll] = wMin
              end
            end
          end
        end

        # neuron pre
        if(j>ncellsI && j<=ncellsI+nPre)
          ll = j-ncellsI #pre
          for k=1:1:nPost # post
            # __ POTENTIATION ___
            if post_spikes[k] ==1
              if SB==0
                new = x_prev[ll]*(A2_p + A3_p*y2_prev[k])+ww_var[ll,k] 
              else
                new = x_prev[ll]*(A2_p + A3_p*y2_prev[k]) *(1-ww_var[ll,k]) + ww_var[ll,k]
              end
              ww_var[ll,k] = new
              if ww_var[ll,k] > wMax
                ww_var[ll,k] = wMax
              end
              if ww_var[ll,k] < wMin
                ww_var[ll,k] = wMin
              end
            end
          end
        end

        mNa[j]  += dmNa(Vprev[j],mNa[j])
        hNa[j]  += dhNa(Vprev[j],hNa[j])
        mKd[j]  += dmKd(Vprev[j],mKd[j])
        mCaT[j] += dmCaT(Vprev[j],mCaT[j])
        hCaT[j] += dhCaT(Vprev[j],hCaT[j])
        mH[j]   += dmH(Vprev[j],mH[j])
        
        if(V[j]>=0. && Vprev[j]<0.)
          Spkt[j,lncells[j]] = t[z]
          Spkt_binary[j,idx_w] = 1
          lncells[j] +=1
        end

        Vprev[j] = copy(V[j])

    end # end of j loop for cortical cells

    if(mod(z,1/dt)==0)
      #VV[:,idx_w] = copy(V)
      ww[:,idx_w] = copy(reshape(ww_var, length(ww_var), 1))


      
      idx_w+=1
    end
    #LFP_C[z] = (1/nPre).*sum(Isyn[ncellsI+nPre:end])
    

  end # end of loop on z
  


  return  Spkt, Spkt_binary,ww
end
