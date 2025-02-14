# 24 GW of BESS by 2030 according to https://www.energy-storage.news/uk-battery-energy-storage-market-to-grow-to-24gw-by-2030-says-rystad-energy/
## Secure economic dispatch of great britain power system with 24 GW of BESS and 5 GW of electrolyzers
#Author: Juan Camilo Castano
#Date: 2024-02-12
#@time begin
## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate()

## Step 1: input data
using CSV
using DataFrames
using YAML
using Plots


data = YAML.load_file(joinpath(@__DIR__, "data.yaml"))

data_generators = CSV.read("Generators_data.csv", DataFrame)
data_electrolyzers = CSV.read("Electrolyzers_data.csv", DataFrame)
data_bess = CSV.read("BESS_data.csv", DataFrame)
ts = CSV.read("Load_data.csv", DataFrame)
data_system_parameters = CSV.read("System_parameters.csv", DataFrame)

## Step 2: create model & pass data to model
using JuMP
using Ipopt
using Gurobi


m = Model(optimizer_with_attributes(Gurobi.Optimizer))


function define_sets!(m::Model, data::Dict, ts::DataFrame)
    #Step 2a: Create sets
    m.ext[:sets] = Dict()
    #Time steps
    J=m.ext[:sets][:J] = 1:24 # time periods

    #Dispatchable generators per type
    IDtype = m.ext[:sets][:IDtype] = [id for id in keys(data["dispatchableGenerators"])]
    #Dispatchable generators per unit
    ID = Array{Union{Nothing,String}}(nothing,0)
    for idtype in IDtype, i in 1:data["dispatchableGenerators"][idtype]["numberOfUnits"]
       ID = m.ext[:sets][:ID] = push!(ID,string(idtype,"_$(i)"))
    end

    #Electrolyzer per type
    IDtype_E = m.ext[:sets][:IDtype] = [id for id in keys(data["Electrolyzer"])]
    #Electrolyzers per unit
    ID_E = Array{Union{Nothing,String}}(nothing,0)
    for idtype_E in IDtype_E, i in 1:data["Electrolyzer"][idtype_E]["numberOfUnits"]
       ID_E = m.ext[:sets][:ID_E] = push!(ID_E,string(idtype_E,"_$(i)"))
    end

    #batteries energy storage systems
      ID_BESS = Array{Union{Nothing,String}}(nothing,0)
      for i in 1:data["BESS"]["NumberOfUnits"]
         ID_BESS = m.ext[:sets][:ID_BESS] = push!(ID_BESS,string("BESS","_$(i)"))
      end

      # Variable generators
      IV = m.ext[:sets][:IV] = [iv for iv in keys(data["variableGenerators"])]

      # All variable and dispatchable generators, per unit
      I = m.ext[:sets][:I] = union(m.ext[:sets][:IV],m.ext[:sets][:ID])
   return m
end

# Step 2b: add time series

function process_time_series_data!(m::Model, ts::DataFrame)
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:D] = ts.Load_typical_winter[1:24]
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
   #extract sets
   ID = m.ext[:sets][:ID]
   ID_E = m.ext[:sets][:ID_E]
   ID_BESS = m.ext[:sets][:ID_BESS]
   IV = m.ext[:sets][:IV]

   #create parameters dictionary
   m.ext[:parameters] = Dict()

   #System parameters
   d=data
   m.ext[:parameters][:rocofmax]=data["rocofmax"]
   m.ext[:parameters][:hydrogenCost]=data["hydrogenCost"]
   m.ext[:parameters][:FO]=data["FO"]
   m.ext[:parameters][:deltaf]=data["deltaf"]


   #parameters of dispatchable generators per unit
   d = data["dispatchableGenerators"]

   #Parameter Fuel cost
   m.ext[:parameters][:FCOST]=Dict()
   m.ext[:parameters][:GmaxD]=Dict()
   m.ext[:parameters][:GminD]=Dict()
   m.ext[:parameters][:Dtg]=Dict()
   m.ext[:parameters][:res_cost_g]=Dict()
   m.ext[:parameters][:inertia_Constant]=Dict()
   for i in ID
      if length(i)==6 
         m.ext[:parameters][:FCOST][i] =  d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
         m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"] 
         m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
         m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
         m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
         m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-2)]["innertiaCostant"]
      else 
         if length(i)==7 
            m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-3)]["marginalcost"] #	β Marginal fuel cost
            m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-3)]["maxPowerOutput"]
            m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-3)]["minStableOperatingPoint"]
            m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-3)]["deploymentTime"]
            m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-3)]["reserveCosts"]
            m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-3)]["innertiaCostant"]
         else
            if length(i)==9
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
               m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
               m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-2)]["innertiaCostant"]
            else
               if length(i)==8
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-4)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-4)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-4)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-4)]["deploymentTime"] 
               m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-4)]["reserveCosts"]
               m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-4)]["innertiaCostant"]          
               else
               end         
            end
         end
   end
   end

   converted_dict_FCOST = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:FCOST])
   m.ext[:parameters][:FCOST] = converted_dict_FCOST
   converted_dict_GmaxD = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:GmaxD])
   m.ext[:parameters][:GmaxD] = converted_dict_GmaxD
   converted_dict_GminD = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:GminD])
   m.ext[:parameters][:GminD] = converted_dict_GminD
   converted_dict_Dtg = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:Dtg])
   m.ext[:parameters][:Dtg] = converted_dict_Dtg
   converted_dict_res_cost_g = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:res_cost_g]) 
   m.ext[:parameters][:res_cost_g] = converted_dict_res_cost_g
   converted_dict_inertia_Constant = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:inertia_Constant])
   m.ext[:parameters][:inertia_Constant] = converted_dict_inertia_Constant

   #Parameter variable generators

   #Parameter electrolyzers
   d = data["Electrolyzer"]
   m.ext[:parameters][:PEmax] = Dict(i => d[SubString(i,1:length(i)-2)]["installed_capacity"] for i in ID_E)
   m.ext[:parameters][:PEmin] = Dict(i => d[SubString(i,1:length(i)-2)]["minimun_consumption"] for i in ID_E)
   m.ext[:parameters][:Eeff] = Dict(i => d[SubString(i,1:length(i)-2)]["efficiency"] for i in ID_E)
   m.ext[:parameters][:heffc]=Dict(i => d[SubString(i,1:length(i)-2)]["effc"] for i in ID_E)
   m.ext[:parameters][:heffd]=Dict(i => d[SubString(i,1:length(i)-2)]["effd"] for i in ID_E)
   m.ext[:parameters][:Eload_factor] = Dict(i => d[SubString(i,1:length(i)-2)]["load_factor_electrolyzer"] for i in ID_E) # load factor of electrolyzers
   m.ext[:parameters][:Max_h_f] = Dict(i => d[SubString(i,1:length(i)-2)]["max_hydrogen_flow"] for i in ID_E) #maximum hydrogen flow in kg/h
   m.ext[:parameters][:Max_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["max_hydrogen_storage"] for i in ID_E) #maximum hydrogen storage in kg
   m.ext[:parameters][:Min_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["min_hydrogen_storage"] for i in ID_E) #minimum hydrogen storage in kg
   m.ext[:parameters][:Ini_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["initial_hydrogen_storage"] for i in ID_E) #initial hydrogen storage in kg
   m.ext[:parameters][:End_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["final_hydrogen_storage"] for i in ID_E) #final hydrogen storage in kg
   m.ext[:parameters][:Dte] = Dict(i => d[SubString(i,1:length(i)-2)]["deploymentTime"] for i in ID_E) #final hydrogen storage in kg
   m.ext[:parameters][:res_cost_e]= Dict(i => d[SubString(i,1:length(i)-2)]["reserveCosts"] for i in ID_E) #final hydrogen storage in kg
   



   

   #Parameter BESS
   d = data
   m.ext[:parameters][:PBmax]=Dict()
   m.ext[:parameters][:EBmax]=Dict()
   m.ext[:parameters][:DOD_max]=Dict()
   m.ext[:parameters][:Beffc]=Dict()
   m.ext[:parameters][:Beffd]=Dict()
   m.ext[:parameters][:Ini_e_b]=Dict()
   m.ext[:parameters][:End_e_b]=Dict()
   m.ext[:parameters][:Dtb]=Dict()
   m.ext[:parameters][:res_cost_b]=Dict()
   
   for i in ID_BESS
      if length(i)==6 
         m.ext[:parameters][:PBmax][i] =  d[SubString(i,1:length(i)-2)]["Pmax"]
         m.ext[:parameters][:EBmax][i] =  d[SubString(i,1:length(i)-2)]["Emax"]
         m.ext[:parameters][:DOD_max][i] =  d[SubString(i,1:length(i)-2)]["DODmax"]
         m.ext[:parameters][:Beffc][i] =  d[SubString(i,1:length(i)-2)]["effc"]
         m.ext[:parameters][:Beffd][i] =  d[SubString(i,1:length(i)-2)]["effd"]
         m.ext[:parameters][:Ini_e_b][i] =  d[SubString(i,1:length(i)-2)]["Einit"]
         m.ext[:parameters][:End_e_b][i] =  d[SubString(i,1:length(i)-2)]["Efinal"]
         m.ext[:parameters][:Dtb][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
         m.ext[:parameters][:res_cost_b][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
      else 
         if length(i)==7 
         m.ext[:parameters][:PBmax][i] =  d[SubString(i,1:length(i)-3)]["Pmax"]
         m.ext[:parameters][:EBmax][i] =  d[SubString(i,1:length(i)-3)]["Emax"]
         m.ext[:parameters][:DOD_max][i] =  d[SubString(i,1:length(i)-3)]["DODmax"]
         m.ext[:parameters][:Beffc][i] =  d[SubString(i,1:length(i)-3)]["effc"]
         m.ext[:parameters][:Beffd][i] =  d[SubString(i,1:length(i)-3)]["effd"]
         m.ext[:parameters][:Ini_e_b][i] =  d[SubString(i,1:length(i)-3)]["Einit"]
         m.ext[:parameters][:End_e_b][i] =  d[SubString(i,1:length(i)-3)]["Efinal"]
         m.ext[:parameters][:Dtb][i] =  d[SubString(i,1:length(i)-3)]["deploymentTime"]
         m.ext[:parameters][:res_cost_b][i] =  d[SubString(i,1:length(i)-3)]["reserveCosts"]
         else
            
         end
   end
   end


   converted_dict_PBmax = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:PBmax])
   m.ext[:parameters][:PBmax] = converted_dict_PBmax

   converted_dict_EBmax = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:EBmax])
   m.ext[:parameters][:EBmax] = converted_dict_EBmax

   converted_dict_Ini_e_b = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:Ini_e_b])
   m.ext[:parameters][:Ini_e_b] = converted_dict_Ini_e_b

   converted_dict_End_e_b = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:End_e_b])
   m.ext[:parameters][:End_e_b] = converted_dict_End_e_b

   #OJO: no se convirtio el diccionario de eficiencia, deployment time, DOD and reserve costs en Dict{String, Int64}
   #return model
   return m
end

#call functions
define_sets!(m, data, ts)
process_time_series_data!(m, ts)
process_parameters!(m, data)


function build_model_interval_1!(m::Model)
   # Create m.ext entries "variables", "expressions" and "constraints"
   m.ext[:variables] = Dict()
   m.ext[:expressions] = Dict()
   m.ext[:constraints] = Dict()

   # Extract sets
   I = m.ext[:sets][:I]
   J = m.ext[:sets][:J]
   ID= m.ext[:sets][:ID]
   IV = m.ext[:sets][:IV]
   ID_E = m.ext[:sets][:ID_E]
   ID_BESS = m.ext[:sets][:ID_BESS]

   # Extract time series data and convert them in PU values
   D= m.ext[:timeseries][:D]
   Pbase=maximum(D) 
   D= m.ext[:timeseries][:D]/Pbase

   #Extract paremeters of the system   m.ext[:parameters][:rocofmax]=data["rocofmax"]
   rocofmax = m.ext[:parameters][:rocofmax]
   hydrogenCost = m.ext[:parameters][:hydrogenCost]
   FO = m.ext[:parameters][:FO]
   deltaf = m.ext[:parameters][:deltaf]

   # Extract parameters Generators
   CostFuel=m.ext[:parameters][:FCOST]
   res_cost_g=m.ext[:parameters][:res_cost_g]

   GmaxD=m.ext[:parameters][:GmaxD]
   GmaxD=Dict(key => value / Pbase for (key, value) in GmaxD)

   GminD=m.ext[:parameters][:GminD]
   GminD =Dict(key => value / Pbase for (key, value) in GminD)

   #Dtg=m.ext[:parameters][:Dtg]
   Dtg=15

   inertia_Constant=m.ext[:parameters][:inertia_Constant]



   #Define Mbase
   Max_h_s = m.ext[:parameters][:Max_h_s]
   Mbase=maximum( Max_h_s)[2]
   Max_h_s =Dict(key => value /Mbase  for (key, value) in  Max_h_s)


   #Extra parameters Electrolyzer

   PEmax = m.ext[:parameters][:PEmax]
   PEmax=Dict(key => value / Pbase for (key, value) in PEmax)

   PEmin = m.ext[:parameters][:PEmin]
   PEmin = Dict(key=> value / Pbase for(key,value) in PEmin)

   Eeff = m.ext[:parameters][:Eeff]
   Eeff =Dict(key => value*Mbase/Pbase for (key, value) in Eeff)

   heffc = m.ext[:parameters][:heffc]
   heffc =Dict(key => value /Mbase  for (key, value) in heffc)

   heffd = m.ext[:parameters][:heffd]
   heffd =Dict(key => value /Mbase  for (key, value) in heffd)

   Eload_factor = m.ext[:parameters][:Eload_factor]

   Max_h_f = m.ext[:parameters][:Max_h_f]
   Max_h_f = Dict(key => value /Mbase  for (key, value) in Max_h_f)

   Min_h_s = m.ext[:parameters][:Min_h_s]
   Min_h_s =Dict(key => value /Mbase  for (key, value) in  Min_h_s)


   Ini_h_s = m.ext[:parameters][:Ini_h_s]
   Ini_h_s =Dict(key => value /Mbase  for (key, value) in Ini_h_s)


   End_h_s = m.ext[:parameters][:End_h_s]
   End_h_s =Dict(key => value /Mbase  for (key, value) in End_h_s)


   #Dte = m.ext[:parameters][:Dte]
   Dte=0.2
   res_cost_e= m.ext[:parameters][:res_cost_e]

   
   
   #Extra parameteres BESS

   PBmax = m.ext[:parameters][:PBmax]
   PBmax=Dict(key => value /Pbase  for (key, value) in PBmax)

   EBmax = m.ext[:parameters][:EBmax]
   EBmax =Dict(key => value /Pbase  for (key, value) in EBmax)

   DOD_max = m.ext[:parameters][:DOD_max]
   Beffc = m.ext[:parameters][:Beffc]
   Beffd = m.ext[:parameters][:Beffd]
   Ini_e_b = m.ext[:parameters][:Ini_e_b]
   Ini_e_b =Dict(key => value /Pbase  for (key, value) in Ini_e_b)

   End_e_b = m.ext[:parameters][:End_e_b]
   End_e_b =Dict(key => value /Pbase  for (key, value) in End_e_b)

   #Dtb = m.ext[:parameters][:Dtb]
   Dtb=0.2
   res_cost_b = m.ext[:parameters][:res_cost_b]

   # create variables 
   g = m.ext[:variables][:g] = @variable(m, [i=ID,j=J],lower_bound=GminD[i], base_name="generation") #Power generation generators
   rg = m.ext[:variables][:rg] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="rg") #Reserve provided by generators
   rb = m.ext[:variables][:rb] = @variable(m, [i=ID_BESS,j=J],lower_bound=0,base_name="rb") #Reserve provided by batteries
   re = m.ext[:variables][:re] = @variable(m, [i=ID_E,j=J],lower_bound=0, base_name="re") #REserve provided by electrolyzers
   pl = m.ext[:variables][:pl] = @variable(m, [j=J],base_name="pl") #loss of generation
   pbc = m.ext[:variables][:pbc] = @variable(m, [i=ID_BESS,j=J],lower_bound=0,upper_bound=PBmax[i], base_name="pbc") #Charging power of the batteries
   pbd = m.ext[:variables][:pbd] = @variable(m, [i=ID_BESS,j=J],lower_bound=0,upper_bound=PBmax[i], base_name="pbd") #Discharging power of the batteries
   eb = m.ext[:variables][:eb] = @variable(m, [i=ID_BESS,j=J], lower_bound= EBmax[i]*(1-DOD_max[i]), upper_bound=EBmax[i] , base_name="eb") #Energy bounds of the batteries
   hfe= m.ext[:variables][:hfe] = @variable(m, [i=ID_E,j=J],lower_bound=0, upper_bound= Max_h_f[i], base_name="hfe") #Hydrogen flow limit of the hydrogen produced by electrolyzers
   hfg = m.ext[:variables][:hfg] = @variable(m, [i=ID_E,j=J],lower_bound=-Max_h_f[i],upper_bound= Max_h_f[i],base_name="hfg") #Hydrogen flow limit of the hydrogen flowing trhow the hydrogen pipeline
   hfsc = m.ext[:variables][:hfsc] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound= Max_h_f[i],base_name="hfsc") #Hydrogen flow storage charging
   hfsd = m.ext[:variables][:hfsd] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound= Max_h_f[i],base_name="hfsd") #Hydrogen flow storage discharging
   pe = m.ext[:variables][:pe] = @variable(m,  [i=ID_E,j=J],lower_bound= PEmin[i], upper_bound=PEmax[i], base_name="pe") #Power consumption electrolyzer
   hss = m.ext[:variables][:hss] = @variable(m, [i=ID_E,j=J],lower_bound=  Min_h_s[i], upper_bound= Max_h_s[i], base_name="hss") #hydrogen storage limit



   #create affine expressions

   g_costs=m.ext[:expressions][:g_costs] = @expression(m, [i=ID,j=J],g[i,j]*CostFuel[i]*Pbase
   )
   rg_costs=m.ext[:expressions][:rg_costs] = @expression(m, [i=ID,j=J],rg[i,j]*res_cost_g[i]*Pbase
      )
   rb_costs=m.ext[:expressions][:rb_costs] = @expression(m, [i=ID_BESS,j=J],rb[i,j]*res_cost_b[i]*Pbase
      )
   re_costs=m.ext[:expressions][:re_costs] = @expression(m, [i=ID_E,j=J],re[i,j]*res_cost_e[i]*Pbase
      )
   h_costs=m.ext[:expressions][:h_costs] = @expression(m, [i=ID_E,j=J],hfg[i,j]*hydrogenCost*Mbase
      )

   #Create objective function

   obj= m.ext[:objective] = @objective(m,Min, sum(g_costs)+sum(rg_costs)+sum(rb_costs)+sum(re_costs)+sum(h_costs)
   )

   #Constraints (market clearing constraint)
   con1=m.ext[:constraints][:con1] = @constraint(m, [j=J],
   sum(g[i,j] for i in ID) -sum(pbc[i,j] for i in ID_BESS)+sum(pbd[i,j] for i in ID_BESS) == D[j]+sum(pe[i,j] for i in ID_E)
   )

   #Constraint upper bound generators power
   con2=m.ext[:constraints][:con2] = @constraint(m, [i=ID,j=J],
   g[i,j]+rg[i,j].<=GmaxD[i]
   )

   #Constraint lost of generation

   con3=m.ext[:constraints][:con3] = @constraint(m, [i=ID,j=J],
   pl[j].>=g[i,j]
   )

   #Constraint upper bound reserve provided by BESS
   con4=m.ext[:constraints][:con4] = @constraint(m, [i=ID_BESS,j=J],
   rb[i,j].<=PBmax[i]+pbc[i,j]-pbd[i,j]
   )

   #Constraint upper bound reserve provided by Electrolyzer
   con5=m.ext[:constraints][:con5] = @constraint(m, [i=ID_E,j=J],re[i,j].<=pe[i,j]-PEmin[i])

   #Constraint end energy value of the batteries
   con6=m.ext[:constraints][:con6] = @constraint(m, [i=ID_BESS,j=J[end]],End_e_b[i]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])


   #Constraint initial value energy of the batteries
   con7=m.ext[:constraints][:con7] = @constraint(m, [i=ID_BESS,j=J[1]],eb[i,j+1]-Ini_e_b[i]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])

   #Constraint charging-discharging batteries
   con8=m.ext[:constraints][:con8] = @constraint(m, [i=ID_BESS,j=J[2:end-1]],eb[i,j+1]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])

   ##Constraint Electrolyzer
   #Constraint  hydrogen production electrolyzer
   con9=m.ext[:constraints][:con9] = @constraint(m, [i=ID_E,j=J],hfe[i,j]==pe[i,j]/Eeff[i])

   #Mass hydrogen equation
   #The term Eload_factorl[i]*PEmax[i]/Eeff[i] correspond to the hydrogen demand
   con10=m.ext[:constraints][:con10] = @constraint(m, [i=ID_E,j=J],hfe[i,j]+hfsc[i,j]-hfsd[i,j]==hfg[i,j]+Eload_factor[i]*PEmax[i]/(Eeff[i]*Mbase/Pbase))

   #Hydrogen storage constraints


   #Constraint end hydrogen value of the hydrogen storage
   con11=m.ext[:constraints][:con11] = @constraint(m, [i=ID_E,j=J[end]],End_h_s[i]==hss[i,j]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])


   #Constraint initial value of the hydrogen storage
   con12=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[1]],hss[i,j+1]==Ini_h_s[i]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])

   #Constraint charging-discharging of the hydrogen storage
   con13=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[2:end-1]],hss[i,j+1]==hss[i,j]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])

   #constraint ROCOF
   con17=m.ext[:constraints][:con17] = @constraint(m, [i=ID,j=J], pl[j]*FO/(2*(sum(values(inertia_Constant))-inertia_Constant[i])).<=rocofmax)


   #Constraint maximum frequency variation
   con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(values(inertia_Constant))-inertia_Constant[i])/FO*(sum(rb[:, j])/Dtb+sum(re[:, j])/Dte+(sum(rg[:, j])-rg[i, j])/Dtg)).>=pl[j]^2/(4*deltaf))

   #constraints nadir occurrence time

   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=0)

   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<= 0.00000000000001+sum(re[:, j])+sum(rb[:,j])*(Dte/Dtb)+(sum(rg[:, j])-rg[i, j])*Dte/Dtg)

   #=
   #Restricciones intervalo II
   #Constraint maximum frequency variation
   con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(values(inertia_Constant))-inertia_Constant[i])/FO-sum(re[:,j])*Dte/(4*deltaf))*(sum(rb[:, j])/Dtb+(sum(rg[:, j])-rg[i, j])/Dtg).>=(pl[j]-sum(re[:, j]))^2/(4*deltaf))

   #constraints nadir occurrence time

   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(re[:, j])  + sum(rb[:, j])*Dte/Dtb + (sum(rg[:, j]) - rg[i, j]) * Dte / Dtg)

   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000000000000001 + sum(re[:, j]) + sum(rb[:, j]) + (sum(rg[:, j]) - rg[i, j]) * Dtb / Dtg)
   

 
   #Restricciones intervalo III

   #Constraint maximum frequency variation
   con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(values(inertia_Constant))-inertia_Constant[i])/FO-(sum(rb[:, j])*Dtb/(4*deltaf))-(sum(re[:, j])*Dte/(4*deltaf)))*((sum(rg[:, j])-rg[i, j])/Dtg).>=(pl[j]-sum(rb[:, j])-sum(re[:, j]))^2/(4*deltaf))

   #constraints nadir occurrence time

   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j])*Dtb/Dtg)

   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000001 + sum(rb[:,j])+sum(re[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j]))
   =#
end


function build_model_interval_2!(m::Model)
   #Extract paremeters of the system   m.ext[:parameters][:rocofmax]=data["rocofmax"]
   build_model_interval_1!(m)
   FO = m.ext[:parameters][:FO]
   deltaf = m.ext[:parameters][:deltaf]
   inertia_Constant=m.ext[:parameters][:inertia_Constant]
   Dtb=0.2
   Dte=0.2
   Dtg=15

   # Extract sets
   J = m.ext[:sets][:J]
   ID= m.ext[:sets][:ID]
   I= m.ext[:sets][:I]

   #Extract variables
   rb = m.ext[:variables][:rb]
   re = m.ext[:variables][:re]
   rg = m.ext[:variables][:rg]
   pl = m.ext[:variables][:pl]

   
   for j in J
      for i in ID
         delete(m,m.ext[:constraints][:con14][i,j])
         delete(m,m.ext[:constraints][:con15][i,j])
         delete(m,m.ext[:constraints][:con16][i,j])
      end
   end

   #Constraint maximum frequency variation
   con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(values(inertia_Constant))-inertia_Constant[i])/FO-sum(re[:,j])*Dte/(4*deltaf))*(sum(rb[:, j])/Dtb+(sum(rg[:, j])-rg[i, j])/Dtg).>=(pl[j]-sum(re[:, j]))^2/(4*deltaf))

   #constraints nadir occurrence time

   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(re[:, j])  + sum(rb[:, j])*Dte/Dtb + (sum(rg[:, j]) - rg[i, j]) * Dte / Dtg)

   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000000000000001 + sum(re[:, j]) + sum(rb[:, j]) + (sum(rg[:, j]) - rg[i, j]) * Dtb / Dtg)

   return (m)
end


function build_model_interval_3!(m::Model)
      build_model_interval_1!(m)
      #Extract paremeters of the system   m.ext[:parameters][:rocofmax]=data["rocofmax"]
      FO = m.ext[:parameters][:FO]
      deltaf = m.ext[:parameters][:deltaf]
      inertia_Constant=m.ext[:parameters][:inertia_Constant]
      Dtb=0.2
      Dte=0.2
      Dtg=15

      # Extract sets
      J = m.ext[:sets][:J]
      ID= m.ext[:sets][:ID]
      I= m.ext[:sets][:I]

      #Extract variables
      rb = m.ext[:variables][:rb]
      re = m.ext[:variables][:re]
      rg = m.ext[:variables][:rg]
      pl = m.ext[:variables][:pl]

      
      for j in J
         for i in ID
            delete(m,m.ext[:constraints][:con14][i,j])
            delete(m,m.ext[:constraints][:con15][i,j])
            delete(m,m.ext[:constraints][:con16][i,j])
         end
      end

   #Constraint maximum frequency variation
   con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(values(inertia_Constant))-inertia_Constant[i])/FO-(sum(rb[:, j])*Dtb/(4*deltaf))-(sum(re[:, j])*Dte/(4*deltaf)))*((sum(rg[:, j])-rg[i, j])/Dtg).>=(pl[j]-sum(rb[:, j])-sum(re[:, j]))^2/(4*deltaf))

   #constraints nadir occurrence time

   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j])*Dtb/Dtg)

   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000001 + sum(rb[:,j])+sum(re[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j]))
   
   return (m)
end

build_model_interval_1!(m)

open("output_1.txt", "w") do file
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
       optimize!(m)
       opt_time = solve_time(m)
   end
end





build_model_interval_2!(m)

open("output_2.txt", "w") do file
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
       optimize!(m)
       opt_time = solve_time(m)
   end
end





build_model_interval_3!(m)

open("output_3.txt", "w") do file
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
       optimize!(m)
       opt_time = solve_time(m)
   end
end





Max_h_s = m.ext[:parameters][:Max_h_s]
Mbase=maximum( Max_h_s)[2]
D=m.ext[:timeseries][:D]
Pbase=maximum(D)

g = value.(m.ext[:variables][:g])*Pbase
rg= value.(m.ext[:variables][:rg])*Pbase
re= value.(m.ext[:variables][:re])*Pbase
rb= value.(m.ext[:variables][:rb])*Pbase
pl= value.(m.ext[:variables][:pl])*Pbase
pbc= value.(m.ext[:variables][:pbc])*Pbase
pbd= value.(m.ext[:variables][:pbd])*Pbase
eb= value.(m.ext[:variables][:eb])*Pbase
pe= value.(m.ext[:variables][:pe])*Pbase
hfe= value.(m.ext[:variables][:hfe])*Mbase
hfg= value.(m.ext[:variables][:hfg])*Mbase
hfsc= value.(m.ext[:variables][:hfsc])*Mbase
hfsd= value.(m.ext[:variables][:hfsd])*Mbase
hss= value.(m.ext[:variables][:hss])*Mbase


# Extract sets
I = m.ext[:sets][:I]
J = m.ext[:sets][:J]
ID= m.ext[:sets][:ID]
IV = m.ext[:sets][:IV]
ID_E = m.ext[:sets][:ID_E]
ID_BESS = m.ext[:sets][:ID_BESS]
Ivec = [i for  i in I]
Bvec = [i for  i in ID_BESS]

gvec = [g[i,j] for  i in ID, j in J]
rgvec = [rg[i,j] for  i in ID, j in J]
revec = [re[i,j] for  i in ID_E, j in J]
rbvec = [rb[i,j] for  i in ID_BESS, j in J]
plvec = [pl[j] for j in J]
pevec = [pe[i,j] for  i in ID_E, j in J]
pbcvec = [pbc[i,j] for  i in ID_BESS, j in J]
pbdvec = [pbd[i,j] for  i in ID_BESS, j in J]
ebvec = [eb[i,j] for  i in ID_BESS, j in J]
hfevec = [hfe[i,j] for  i in ID_E, j in J]
hfgvec = [hfg[i,j] for  i in ID_E, j in J]
hfscvec = [hfsc[i,j] for  i in ID_E, j in J]
hfsdvec = [hfsd[i,j] for  i in ID_E, j in J]
hssvec = [hss[i,j] for  i in ID_E, j in J]

using StatsPlots
p2 = groupedbar(transpose(gvec[:,:]),bar_position = :stack,label=permutedims(Ivec),legend=:outertopright);
plot(p2, layout = (1,2), size=(2400, 1900))
p3 = groupedbar(transpose(rgvec[:,:]),bar_position = :stack,label=permutedims(Bvec),legend=:outertopright);
plot(p3, layout = (1,2), size=(2400, 1900))
# Plot the vectors
plot(hfevec[1,:], label = "hfevec", lw = 2)
plot!(hfgvec[1,:], label = "hfgvec", lw = 2)
plot!(hfscvec[1,:], label = "hfscvec", lw = 2)
plot!(hfsdvec[1,:], label = "hfsdvec", lw = 2)
plot!(hssvec[1,:], label = "hssvec", lw = 2)
# Display the plot
plot!()


plot(hfscvec[1,:], label = "hfscvec", lw = 2)
plot!(hfsdvec[1,:], label = "hfsdvec", lw = 2)
plot!(hssvec[1,:], label = "hssvec", lw = 2)
# Display the plot
plot!()

plot(hfscvec[1,:], label = "hfscvec", lw = 2)
plot!(hfsdvec[1,:], label = "hfsdvec", lw = 2)
plot!(hssvec[1,:], label = "hssvec", lw = 2)
# Display the plot
plot!()

plot(pbcvec[1,:], label = "pbcvec", lw = 2)
plot!(pbdvec[1,:], label = "pbdvec", lw = 2)
plot!(ebvec[1,:], label = "ebvec", lw = 2)
plot!()

#end