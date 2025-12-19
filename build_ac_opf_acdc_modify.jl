function build_ac_opf_acdc_modify!(m::Model)
    # This function builds the polar form of a nonlinear AC power flow formulation

    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    # AC network
    N = m.ext[:sets][:N]
    N_sl = m.ext[:sets][:N_sl]
    B = m.ext[:sets][:B]
    B_ac_fr = m.ext[:sets][:B_ac_fr]
    B_ac_to = m.ext[:sets][:B_ac_to]
    G = m.ext[:sets][:G]
    G_ac = m.ext[:sets][:G_ac]
    L = m.ext[:sets][:L]
    L_ac = m.ext[:sets][:L_ac]
    W_ac = m.ext[:sets][:W_ac]
    B_ac = m.ext[:sets][:B_ac]
    B_arcs = m.ext[:sets][:B_arcs]
    bus_ij = m.ext[:sets][:bus_ij]
    bus_ji = m.ext[:sets][:bus_ji]
    bus_ij_ji = m.ext[:sets][:bus_ij_ji]


    # DC network
    CV = m.ext[:sets][:CV]
    ND = m.ext[:sets][:ND]
    BD = m.ext[:sets][:BD]
    ND_arcs = m.ext[:sets][:ND_arcs]
    CV_arcs = m.ext[:sets][:CV_arcs]   
    BD_dc_fr = m.ext[:sets][:BD_dc_fr]
    BD_dc_to = m.ext[:sets][:BD_dc_to]
    BD_dc = m.ext[:sets][:BD_dc]

    busdc_ij = m.ext[:sets][:busdc_ij]
    busdc_ji = m.ext[:sets][:busdc_ji]
    T= m.ext[:sets][:t]

    S = m.ext[:sets][:S]
    S_ac = m.ext[:sets][:S_ac]
    E= m.ext[:sets][:E]
    E_ac= m.ext[:sets][:E_ac]

    # Extract parameters
    # AC network
    gen_bus = m.ext[:parameters][:gen_bus]
    load_bus = m.ext[:parameters][:load_bus]
    
    #shunt_bus = m.ext[:parameters][:shunt_bus]
    vmmin = m.ext[:parameters][:vmmin]
    vmmax = m.ext[:parameters][:vmmax]
    vamin = m.ext[:parameters][:vamin]
    vamax = m.ext[:parameters][:vamax]
    rb =  m.ext[:parameters][:rb]
    xb =  m.ext[:parameters][:xb] 
    gb =  m.ext[:parameters][:gb]
    bb =  m.ext[:parameters][:bb] 
    #gs =  m.ext[:parameters][:gs]
    #bs =  m.ext[:parameters][:bs] 
    gfr = m.ext[:parameters][:gb_sh_fr]
    bfr = m.ext[:parameters][:bb_sh_fr]
    gto = m.ext[:parameters][:gb_sh_to]
    bto = m.ext[:parameters][:bb_sh_to]
    pmaxbranch = m.ext[:parameters][:pmaxbranch]
    angmin = m.ext[:parameters][:angmin]
    angmax = m.ext[:parameters][:angmax]
    b_shift = m.ext[:parameters][:b_shift]
    b_tap = m.ext[:parameters][:b_tap]
    pd = m.ext[:parameters][:pd]
    qd = m.ext[:parameters][:qd]
    pmax = m.ext[:parameters][:pmax]
    pmin = m.ext[:parameters][:pmin]
    #qmax = m.ext[:parameters][:qmax]
    #qmin = m.ext[:parameters][:qmin]
    gen_cost = m.ext[:parameters][:gen_cost]
    ij_ji_ang_max = m.ext[:parameters][:ij_ji_ang_max]
    ij_ji_ang_min = m.ext[:parameters][:ij_ji_ang_min]
    demand = m.ext[:parameters][:demand]
    wind= m.ext[:parameters][:wind]
    upramp = m.ext[:parameters][:upramp]
    downramp = m.ext[:parameters][:downramp]
    MUT = m.ext[:parameters][:MUT]
    MDT = m.ext[:parameters][:MDT]

    # DC network
    # DC bus
    busdc_vm_max = m.ext[:parameters][:busdc][:vm_max]
    busdc_vm_min = m.ext[:parameters][:busdc][:vm_min]
    busdc_vm_set = m.ext[:parameters][:busdc][:vm_set]
    busdc_c = m.ext[:parameters][:busdc][:c]
    busdc_p = m.ext[:parameters][:busdc][:p]

    # Converters
    conv_busdc = m.ext[:parameters][:convdc][:busdc]
    conv_bus = m.ext[:parameters][:convdc][:bus]
    conv_status = m.ext[:parameters][:convdc][:status]
    conv_loss_a = m.ext[:parameters][:convdc][:loss_a]
    conv_loss_b = m.ext[:parameters][:convdc][:loss_b]
    conv_loss_c_inv = m.ext[:parameters][:convdc][:loss_c_inv]
    conv_loss_c_rec = m.ext[:parameters][:convdc][:loss_c_rec]
    conv_p_ac_max = m.ext[:parameters][:convdc][:p_ac_max]
    conv_p_ac_min = m.ext[:parameters][:convdc][:p_ac_min]
    #conv_q_ac_max = m.ext[:parameters][:convdc][:q_ac_max]
    #conv_q_ac_min = m.ext[:parameters][:convdc][:q_ac_min]
    conv_p_dc_max = m.ext[:parameters][:convdc][:p_dc_max]
    conv_p_dc_min = m.ext[:parameters][:convdc][:p_dc_min]
    conv_droop = m.ext[:parameters][:convdc][:droop]
    conv_i_max = m.ext[:parameters][:convdc][:i_max]
    conv_vm_min = m.ext[:parameters][:convdc][:vm_min]
    conv_vm_max = m.ext[:parameters][:convdc][:vm_max]
    conv_vm_dc_set = m.ext[:parameters][:convdc][:vm_dc_set]
    conv_p_g = m.ext[:parameters][:convdc][:p_g]
    #conv_q_g = m.ext[:parameters][:convdc][:q_g]
    conv_b_f = m.ext[:parameters][:convdc][:b_f]
    conv_tf_r = m.ext[:parameters][:convdc][:r_tf]
    conv_tf_g = m.ext[:parameters][:convdc][:g_tf]
    conv_tf_x = m.ext[:parameters][:convdc][:x_tf]
    conv_tf_b = m.ext[:parameters][:convdc][:b_tf]
    conv_pr_r = m.ext[:parameters][:convdc][:r_pr]
    conv_pr_g = m.ext[:parameters][:convdc][:g_pr]
    conv_pr_x = m.ext[:parameters][:convdc][:x_pr]
    conv_pr_b = m.ext[:parameters][:convdc][:b_pr]
    conv_tf_tap = m.ext[:parameters][:convdc][:tap_tf]
    conv_is_tf = m.ext[:parameters][:convdc][:is_tf]
    conv_is_pr = m.ext[:parameters][:convdc][:is_pr]
    conv_is_filter = m.ext[:parameters][:convdc][:is_filter]

    # DC branches
    brdc_rate_a = m.ext[:parameters][:branchdc][:rate_a]
    brdc_rate_b = m.ext[:parameters][:branchdc][:rate_b]
    brdc_rate_c = m.ext[:parameters][:branchdc][:rate_c]
    brdc_status = m.ext[:parameters][:branchdc][:status]
    brdc_r = m.ext[:parameters][:branchdc][:r]
    brdc_g = m.ext[:parameters][:branchdc][:g]
    brdc_l = m.ext[:parameters][:branchdc][:l]
    brdc_dcpoles = m.ext[:parameters][:branchdc][:dcpoles]

    # Electrolyzer
    Epmax = m.ext[:parameters][:Epmax]
    Epmin = m.ext[:parameters][:Epmin]
    Eefficiency = m.ext[:parameters][:Eefficiency]
    Eloadfactor = m.ext[:parameters][:Eloadfactor] 
    Eflowmax = m.ext[:parameters][:Eflowmax]
    Estoragemax = m.ext[:parameters][:Estoragemax]
    Estoragemin = m.ext[:parameters][:Estoragemin]
    Estorageinitial = m.ext[:parameters][:Estorageinitial]
    Estorageend = m.ext[:parameters][:Estorageend]
    Edeployment = m.ext[:parameters][:Edeployment]
    Ereservecost = m.ext[:parameters][:Ereservecost]
    Eeficiencycarga = m.ext[:parameters][:Eeficiencycarga]
    Eeficiencydischarge = m.ext[:parameters][:Eeficiencydischarge]
    Estartupcost = m.ext[:parameters][:Estartupcost]
    Ecompressorpower = m.ext[:parameters][:Ecompressorpower]
    electrolyzer_bus = m.ext[:parameters][:electrolyzer_bus]
    Hydrogencost = m.ext[:parameters][:Hydrogencost]
    baseKG = m.ext[:parameters][:baseKG]

    #storage
    Spmax = m.ext[:parameters][:Spmax]
    Sstoragemax = m.ext[:parameters][:Sstoragemax]
    Sdod = m.ext[:parameters][:Sdod]
    Sefficiencycarga = m.ext[:parameters][:Sefficiencycarga]
    Sefficiencydischarge = m.ext[:parameters][:Sefficiencydischarge]
    Senergyinitial = m.ext[:parameters][:Senergyinitial]
    Senergyend = m.ext[:parameters][:Senergyend]
    Sdeployment = m.ext[:parameters][:Sdeployment]
    Sreservecost = m.ext[:parameters][:Sreservecost]
    storage_bus = m.ext[:parameters][:storage_bus]              



    ##### Create variables 
    # AC components
    # Bus variables



    va = m.ext[:variables][:va] = @variable(m, [i=N,t=T], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va") # voltage angle

    # Generator variables
    pg = m.ext[:variables][:pg] = @variable(m, [g=G,t=T], base_name = "pg") # active power generatio

    # Branch variables
    pb = m.ext[:variables][:pb] = @variable(m, [(b,i,j) in B_ac,t=T],lower_bound = -pmaxbranch[b], upper_bound = pmaxbranch[b], base_name = "pb") # from side active power flow (i->j)


    # Status variable generators
    zg = m.ext[:variables][:zg] = @variable(m, [g=G,t=T],binary= true , base_name = "zg") # status variable generator
    betag = m.ext[:variables][:betag] = @variable(m, [g=G,t=T], binary= true, base_name = "betag") #
    gammag = m.ext[:variables][:gammag] = @variable(m, [g=G,t=T], binary= true, base_name = "gammag") #
    

    # DC components
    
    # Branches
    brdc_p = m.ext[:variables][:brdc_p] = @variable(m, [(d,e,f)=BD_dc, t=T], lower_bound=-brdc_rate_a[d], upper_bound=brdc_rate_a[d], base_name="brdc_p")
  
    # Converters
    conv_p_ac = m.ext[:variables][:conv_p_ac] = @variable(m, [cv=CV,t=T], lower_bound=-conv_p_ac_max[cv], upper_bound=conv_p_ac_max[cv], base_name="conv_p_ac") # converter active power
    conv_p_dc = m.ext[:variables][:conv_p_dc] = @variable(m, [cv=CV,t=T], lower_bound=-conv_p_dc_max[cv], upper_bound=conv_p_dc_max[cv], base_name="conv_p_dc") # converter active power

    # Electrolyzer variables
    pe= m.ext[:variables][:pe] = @variable(m, [e=E,t=T], lower_bound=0, upper_bound=Epmax[e], base_name="pe") # Electrolyzer power consumption
    pe_compressor= m.ext[:variables][:pe_compressor] = @variable(m, [e=E,t=T], lower_bound=0, base_name="pe_compressor") # Electrolyzer power consumption for compressor
    pes= m.ext[:variables][:pes] = @variable(m, [e=E,t=T], lower_bound=0, base_name="pes") # Electrolyzer reserve power
    hfe= m.ext[:variables][:hfe] = @variable(m, [e=E,t=T], lower_bound=0, upper_bound=Eflowmax[e], base_name="hfe") # Electrolyzer hydrogen flow
    hfginyect= m.ext[:variables][:hfginyect] = @variable(m, [e=E,t=T], lower_bound=0, upper_bound=Eflowmax[e], base_name="hfginyect") # Electrolyzer hydrogen flow injected to the hydrogen grid
    hfgconsum= m.ext[:variables][:hfgconsum] = @variable(m, [e=E,t=T], lower_bound=0, upper_bound=Eflowmax[e], base_name="hfgconsum") # Electrolyzer hydrogen flow consumed from the hydrogen grid
    hss= m.ext[:variables][:hss] = @variable(m, [e=E,t=T], lower_bound=Estoragemin[e], upper_bound=Estoragemax[e], base_name="hss") # Electrolyzer storage level
    ze= m.ext[:variables][:ze] = @variable(m, [e=E,t=T], binary=true, base_name="ze") # Electrolyzer status
    zesu= m.ext[:variables][:zesu] =@variable(m, [e=E,t=T], binary=true, base_name="zesu") #Electrolyzer start up 
    zestb= m.ext[:variables][:zestb] = @variable(m, [e=E,t=T], binary=true, base_name="zestb") #Electrolyzer stanby indicator

    # Storage variables
    psc = m.ext[:variables][:psc] = @variable(m, [s=S,t=T],lower_bound=0, base_name="psc") #Charging power of the batteries
    psd = m.ext[:variables][:psd] = @variable(m, [s=S,t=T],lower_bound=0, base_name="psd") #Discharging power of the batteries
    es = m.ext[:variables][:es] = @variable(m, [s=S,t=T], lower_bound= Sstoragemax[s]*(1-Sdod[s]), upper_bound=Sstoragemax[s] , base_name="es") #Energy bounds of the batteries
    zs = m.ext[:variables][:zs] = @variable(m, [s=S,t=T], binary=true, base_name="zb") #standby indicator of the batteries

    ##### Objective
    max_gen_ncost = m.ext[:parameters][:gen_max_ncost]
    if max_gen_ncost == 1
        m.ext[:objective] = @objective(m, Min,
                sum(gen_cost[g][1]
                        for g in G)+sum(Hydrogencost[e]*baseKG*(-hfginyect[e,t]+hfgconsum[e,t]) for e in E, t in T)
        )
    elseif max_gen_ncost == 2
       m.ext[:objective] = @objective(m, Min,
            sum(gen_cost[g][1] * pg[g, t] + gen_cost[g][2] 
                    for g in G, t in T)+sum(Hydrogencost[e]*baseKG*(-hfginyect[e,t]+hfgconsum[e,t]) for e in E, t in T)
        )
    elseif max_gen_ncost == 3
        m.ext[:objective] = @NLobjective(m, Min,
                sum(gen_cost[g][1]*pg[g,t]^2 + gen_cost[g][2]*pg[g,t] + gen_cost[g][3]
                        for g in G, t in T)+sum(Hydrogencost[e]*baseKG*(-hfginyect[e,t]+hfgconsum[e,t]) for e in E, t in T)
        )
    elseif max_gen_ncost == 4
        m.ext[:objective] = @NLobjective(m, Min,
                sum(gen_cost[g][1]*pg[g,t]^3 + gen_cost[g][2]*pg[g,t]^2 + gen_cost[g][3]*pg[g,t] + gen_cost[g][4]
                        for g in G, t in T)+sum(Hydrogencost[e]*baseKG*(-hfginyect[e,t]+hfgconsum[e,t]) for e in E, t in T)
        )
    end

    ####################################################################################################
    ####################    AC NETWORK AND AC/DC CONSTRAINTS
    ####################################################################################################
    

    # Bus angle difference limits

        # Voltage angle on reference bus = 0
    m.ext[:constraints][:varef] = @constraint(m, [n_sl=N_sl,t=T], va[n_sl,t] == 0)

    m.ext[:constraints][:power_balance] = @constraint(m, [n=N,t=T],
        sum(pg[g,t] for g in G if gen_bus[g] == n) +wind[n][t]+sum(psd[s,t] for s in S if storage_bus[s] == n)-sum(psc[s,t] for s in S if storage_bus[s] == n)  - sum(pb[(br,i,j),t] for (br,i,j) in B_arcs[n])- sum(conv_p_ac[cv,t] for cv in CV if conv_bus[cv] == n) -sum(pe[e,t] for e in E if electrolyzer_bus[e] == n ) -sum(pe_compressor[e,t] for e in E if electrolyzer_bus[e] == n )-demand[n][t]== 0
        )


    # Power flow constraints in from and to direction

    #It is assumed that the reactance is 0.13 p.u. based on the test system data
    m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr,t=T],
        pb[(b, i, j),t] ==  1/0.13*(va[i,t] - va[j,t])
        )

         #It is assumed that the reactance is 0.13 p.u. based on the test system data
    m.ext[:constraints][:directional_flow_ac] = @constraint(m, [(b,i,j) = B_ac_fr,t=T],
        pb[(b, i, j),t] ==  -pb[(b, j, i),t]
        )

    # # Power flow constraints in from and to direction
    # m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr,t=T], pb[(b, i, j),t] ==  - 1/0.13*  (va[i,t] - va[j,t])) # active power i to j
    # m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to,t=T], pb[(b, j, i),t] ==  - 1/0.13* (va[j,t] - va[i,t])) # active power j to i

    # # Branch angle limits
    # m.ext[:constraints][:thetaij] = @constraint(m, [(b,i,j) = B_ac_fr,t=T], va[i,t] - va[j,t] <= angmax[b])
    # m.ext[:constraints][:thetaji] = @constraint(m, [(b,i,j) = B_ac_fr,t=T], va[i,t] - va[j,t] >= angmin[b])
    # m.ext[:constraints][:thetaij] = @constraint(m, [(b,j,i) = B_ac_to,t=T], va[j,t] - va[i,t] <= angmax[b])
    # m.ext[:constraints][:thetaji] = @constraint(m, [(b,j,i) = B_ac_to,t=T], va[j,t] - va[i,t] >= angmin[b])
  

     # active power i to j (I am not sure if this constraints has to be)
    #m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to,t=T],      
    # pb[(b, i, j),t] ==  1/0.13*(va[i,t] - va[j,t])
    #)


    m.ext[:constraints][:max_gen_power] = @constraint(m, [g=G,t=T],
        pg[g,t] <= pmax[g]*zg[g,t]
        )

    m.ext[:constraints][:min_gen_power] = @constraint(m, [g=G,t=T],
        pg[g,t] >= pmin[g]*zg[g,t]
        )


    # Nodal power balance DC
    m.ext[:constraints][:nodal_p_dc_balance] = @constraint(m, [nd=ND,t=T],
        - sum(conv_p_dc[cv,t] for cv in CV if conv_busdc[cv] == nd)
        - sum(brdc_p[(d,f,e),t] for (d,f,e) in ND_arcs[nd])
        == 0)

    # Converter AC-side and DC-side power balance
        m.ext[:constraints][:conv_p_loss] = @constraint(m, [cv=CV,t=T],
        conv_p_ac[cv,t] + conv_p_dc[cv,t]==0
    )


                    #It is assumed that the reactance is 0.13 p.u. based on the test system data
    m.ext[:constraints][:direccional_flow_dc] = @constraint(m, [(d,f,e) = BD_dc_fr,t=T],
        brdc_p[(d,f,e),t] ==  -brdc_p[(d, e, f),t]
        )

    
       

    #Unit commitment constraints

    #Unit commitment constraints

    Tlabels = collect(T)         # ensure indexable by position
    NT = length(Tlabels)
    #Ramp-up limit for generators
        m.ext[:constraints][:up_and_down_2_a] = @constraint(m, [g in G, k in 2:NT],
        pg[g, Tlabels[k]] - pg[g, Tlabels[k-1]] <= upramp[g]+(pmin[g]- upramp[g])*betag[g, Tlabels[k]]
    )

 
    #Ramp-down limit for generators
          m.ext[:constraints][:up_and_down_2_b] = @constraint(m, [g in G, k in 2:NT],
        pg[g, Tlabels[k-1]] - pg[g, Tlabels[k]] <= downramp[g]+pmin[g]*gammag[g, Tlabels[k]]
    )

 
  #Initial status generators (unit commitment)
        m.ext[:constraints][:up_and_down_3] = @constraint(m, [g in G],
        1 - zg[g, Tlabels[1]] + betag[g, Tlabels[1]] - gammag[g, Tlabels[1]] == 0
    )
    

    #Status transition generators (unit commitment)
    m.ext[:constraints][:up_and_down_4]=@constraint(m, [g in G, k in 2:NT], zg[g, Tlabels[k-1]] - zg[g, Tlabels[k]] + betag[g, Tlabels[k]] - gammag[g, Tlabels[k]] == 0)

    
    #Status generators (unit commitment)
        m.ext[:constraints][:up_and_down_4] = @constraint(m, [g in G, k in 2:NT],
        zg[g, Tlabels[k-1]] - zg[g, Tlabels[k]] + betag[g, Tlabels[k]] - gammag[g, Tlabels[k]] == 0
    )

    #Startup and shutdown exclusivity
        m.ext[:constraints][:up_and_down_4_1] = @constraint(m, [g=G, t=T],
        betag[g,t] + gammag[g,t] <= 1
    )

    #Minimum up time constraint
        m.ext[:constraints][:up_and_down_6] = Dict()
        for g in G
            if MUT[g] > 0  # Only proceed if MUT[g] is positive
                for k in MUT[g]:NT
                    m.ext[:constraints][:up_and_down_6][g, Tlabels[k]] = @constraint(m,
                    zg[g, Tlabels[k]] >= sum(betag[g, Tlabels[k-τ]] for τ in 0:MUT[g]-1)
                    )
                end
            end
        end



    #Minimum down time constraint
    m.ext[:constraints][:up_and_down_7] = Dict()
    for g in G
        if MDT[g] > 0  # Only proceed if MDT[g] is positive
            for k in MDT[g]:NT
                m.ext[:constraints][:up_and_down_7][g, Tlabels[k]] = @constraint(m,
                1 - zg[g, Tlabels[k]] >= sum(gammag[g, Tlabels[k-τ]] for τ in 0:MDT[g]-1)
                )
            end
        end
    end

#Storage constraints

    m.ext[:constraints][:upper_bound_storage_charging] = @constraint(m, [s = S, t = T],
    psc[s,t] <= Spmax[s] * (1 - zs[s,t])
    )

    m.ext[:constraints][:upper_bound_storage_discharging] = @constraint(m, [s = S, t = T],
    psd[s,t] <= Spmax[s] * zs[s,t]
    )


        # End energy value of the batteries
        
        m.ext[:constraints][:end_energy_value] = @constraint(m, [s in S],
            Senergyend[s] - es[s, Tlabels[NT]] ==
                Sefficiencycarga[s] * psc[s, Tlabels[NT]] -
                psd[s, Tlabels[NT]] / Sefficiencydischarge[s]
        )

        # Charging–discharging energy balance of the batteries
        m.ext[:constraints][:energy_balance] = @constraint(m, [s in S, k in 1:NT-1],
        es[s, Tlabels[k+1]] - es[s, Tlabels[k]] ==
            Sefficiencycarga[s] * psc[s, Tlabels[k]] -
            psd[s, Tlabels[k]] / Sefficiencydischarge[s]
    )

# Electrolyzer constraints

    #Hydrogen production constraint
        m.ext[:constraints][:hydrogen_production] = @constraint(m, [e in E, t in T],
        hfe[e,t] == pe[e,t]/Eefficiency[e]-0.05*Epmax[e]*zestb[e,t]/Eefficiency[e]
        )

    # Bounds on electrolyzer power considering standby
    m.ext[:constraints][:e_lower_bound] = @constraint(m, [e = E, t = T],
    pe[e,t] >= Epmin[e] * ze[e,t] + 0.05 * Epmax[e] * zestb[e,t]
    )    

    m.ext[:constraints][:e_upper_bound] = @constraint(m, [e = E, t = T],
    pe[e,t] <= Epmax[e] * ze[e,t] + 0.05 * Epmax[e] * zestb[e,t]
    )

    # Status constraint of electrolyzers (first time period)
    m.ext[:constraints][:e_status_ini] = @constraint(m, [e in E],
        zesu[e, Tlabels[1]] >= (ze[e, Tlabels[1]] - 1) + (zestb[e, Tlabels[1]] - 0)
    )

    # Startup constraint of electrolyzers (subsequent periods)
    m.ext[:constraints][:e_status_subsequent] = @constraint(m, [e in E, k in 2:NT],
    zesu[e, Tlabels[k]] >=
        (ze[e, Tlabels[k]]     - ze[e, Tlabels[k-1]]) +
        (zestb[e, Tlabels[k]]  - zestb[e, Tlabels[k-1]])
    )

    # Standby constraint of electrolyzers
    m.ext[:constraints][:e_standby] = @constraint(m, [e = E, t = T],
    zestb[e,t] + ze[e,t] <= 1
    )

        
    # Initial value of the hydrogen storage (relation between period 1 and 2)
        m.ext[:constraints][:initial_hydrogen_value] = @constraint(m, [e in E],
            hss[e, Tlabels[1]] == Estorageinitial[e]
        )
    # End value of the hydrogen storage
        m.ext[:constraints][:end_hydrogen_value] = @constraint(m, [e in E],
        Estorageend[e] ==
            hss[e, Tlabels[NT]] +
            hfe[e, Tlabels[NT]] -
            hfginyect[e, Tlabels[NT]] / Eeficiencycarga[e] +
            hfgconsum[e, Tlabels[NT]] * Eeficiencydischarge[e] -
            Eloadfactor[e] * Epmax[e] / Eefficiency[e]
        )

    # Charging–discharging of the hydrogen storage
        m.ext[:constraints][:charging_discharging_hydrogen] = @constraint(m, [e in E, k in 1:NT-1],
            hss[e, Tlabels[k+1]] ==
                hss[e, Tlabels[k]] +
                hfe[e, Tlabels[k]] -
                hfginyect[e, Tlabels[k]] / Eeficiencycarga[e] +
                hfgconsum[e, Tlabels[k]] * Eeficiencydischarge[e] -
                Eloadfactor[e] * Epmax[e] / Eefficiency[e]
        )
    #Compresor power constraint
        m.ext[:constraints][:compressor_power] = @constraint(m, [e in E, t in T],
        pe_compressor[e,t] == Ecompressorpower[e]*(hfe[e,t]+hfgconsum[e,t]+hfginyect[e,t])
        )

        


    return m 
end
