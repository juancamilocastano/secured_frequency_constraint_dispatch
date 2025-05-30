# 24 GW of BESS by 2030 according to https://www.energy-storage.news/uk-battery-energy-storage-market-to-grow-to-24gw-by-2030-says-rystad-energy/
## Secure economic dispatch of great britain power system with 24 GW of BESS and 5 GW of electrolyzers
#Author: Juan Camilo Castano
#Date: 2024-02-12
#@time begin
## Step 0: Activate environment - ensure consistency accross computers

#References of grid topologies: Economic evaluation of a power-to-hydrogen system providing frequency regulation reserves: a case study of Denmark
#Hydrogen Electrolyzer Load Modelling for Steady-State Power System Studies
#From green hydrogen to electricity: A review on recent advances, challenges, and opportunities on power-to-hydrogen-to-power systems (they call it supply chains)
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
capacity_factor_renewable=CSV.read("Capacity.csv", DataFrame)
data_system_parameters = CSV.read("System_parameters.csv", DataFrame)

output_file = joinpath(@__DIR__, "data.txt")

# Open the file in write mode and save the data
open(output_file, "w") do file
    write(file, string(data))
end

## Step 2: create model & pass data to model
using JuMP
using Ipopt
using Gurobi


m = Model(optimizer_with_attributes(Gurobi.Optimizer))
#https://www.sciencedirect.com/science/article/pii/S0378779624005649 according to this paper, the tolerant gap is between 0.1% and 0.01%
set_optimizer_attribute(m, "MIPGap", 0.001) 

function define_sets!(m::Model, data::Dict, ts::DataFrame)
    #Step 2a: Create sets
    m.ext[:sets] = Dict()
    #Time steps
    J=m.ext[:sets][:J] = 1:nrow(ts) # time periods

    #Dispatchable generators per type
    IDtype = m.ext[:sets][:IDtype] = [id for id in keys(data["dispatchableGenerators"])]
    #Dispatchable generators per unit
    ID = Array{Union{Nothing,String}}(nothing,0)
    for idtype in IDtype, i in 1:data["dispatchableGenerators"][idtype]["numberOfUnits"]
       ID = m.ext[:sets][:ID] = push!(ID,string(idtype,"_$(i)"))
    end

    ID_Nuclear = Array{Union{Nothing,String}}(nothing,0)
    for i in 1:data["dispatchableGenerators"]["Nuclear"]["numberOfUnits"]
      ID_Nuclear = m.ext[:sets][:ID_Nuclear] = push!(ID_Nuclear,string("Nuclear_$(i)"))
   end

   ID_Pump = Array{Union{Nothing,String}}(nothing,0)
   for i in 1:data["Pump"]["numberOfUnits"]
     ID_Pump = m.ext[:sets][:ID_Pump] = push!(ID_Pump,string("Pump_$(i)"))
  end

   #Defines the set of generators that can provide frequency reserve
   ID_GR=m.ext[:sets][:ID_GR] =setdiff(ID,ID_Nuclear)




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
   number_hours=nrow(ts)
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:D] = ts.Load_typical_winter[1:number_hours]
    m.ext[:timeseries][:SC] =  capacity_factor_renewable.Solar_Capacity_Factor[1:number_hours]
    m.ext[:timeseries][:WC] =  capacity_factor_renewable.Wind_Capacity_Factor[1:number_hours]
    m.ext[:timeseries][:HVC] = capacity_factor_renewable.Hydrogen_variable[1:number_hours]
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
   #extract sets
   ID = m.ext[:sets][:ID]
   ID_E = m.ext[:sets][:ID_E]
   ID_BESS = m.ext[:sets][:ID_BESS]
   IV = m.ext[:sets][:IV]
   ID_Nuclear = m.ext[:sets][:ID_Nuclear]
   ID_GR = m.ext[:sets][:ID_GR]
   ID_Pump = m.ext[:sets][:ID_Pump]

   #create parameters dictionary
   m.ext[:parameters] = Dict()

   #System parameters
   d=data
   m.ext[:parameters][:rocofmax]=data["rocofmax"]
   m.ext[:parameters][:hydrogenCost]=data["hydrogenCost"]
   m.ext[:parameters][:FO]=data["FO"]
   m.ext[:parameters][:deltaf]=data["deltaf"]

   #Parameter renewable
   d=data["variableGenerators"]
   m.ext[:parameters][:Installed_W]=d["Wind"]["installedCapacity"]
   m.ext[:parameters][:Installed_S]=d["Solar"]["installedCapacity"]


   #parameters of dispatchable generators per unit
   d = data["dispatchableGenerators"]

   #Parameter Fuel cost
   m.ext[:parameters][:FCOST]=Dict()
   m.ext[:parameters][:GmaxD]=Dict()
   m.ext[:parameters][:GminD]=Dict()
   m.ext[:parameters][:Dtg]=Dict()
   m.ext[:parameters][:res_cost_g]=Dict()
   m.ext[:parameters][:inertia_Constant]=Dict()
   m.ext[:parameters][:startupCost]=Dict()
   m.ext[:parameters][:MUT]=Dict()
   m.ext[:parameters][:MDT]=Dict()
   m.ext[:parameters][:maxFrdelivarable]=Dict()
   m.ext[:parameters][:notloadCost]=Dict()
   m.ext[:parameters][:upramprate]=Dict()
   m.ext[:parameters][:downramprate]=Dict()



   for i in ID
      if length(i)==6 
         m.ext[:parameters][:FCOST][i] =  d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
         m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"] 
         m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
         m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
         m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
         m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-2)]["innertiaCostant"]
        m.ext[:parameters][:startupCost][i]= d[SubString(i,1:length(i)-2)]["startupCost"]
        m.ext[:parameters][:MUT][i]= d[SubString(i,1:length(i)-2)]["minUpTime"]
        m.ext[:parameters][:MDT][i]= d[SubString(i,1:length(i)-2)]["minDownTime"]
        m.ext[:parameters][:maxFrdelivarable][i]= d[SubString(i,1:length(i)-2)]["maxFrdelivarable"]
        m.ext[:parameters][:notloadCost][i]= d[SubString(i,1:length(i)-2)]["notloadCost"]
        m.ext[:parameters][:upramprate][i]= d[SubString(i,1:length(i)-2)]["upramprate"]
         m.ext[:parameters][:downramprate][i]= d[SubString(i,1:length(i)-2)]["downramprate"]
  
            

  
        

      else 
         if length(i)==7 
            m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-3)]["marginalcost"] #	β Marginal fuel cost
            m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-3)]["maxPowerOutput"]
            m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-3)]["minStableOperatingPoint"]
            m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-3)]["deploymentTime"]
            m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-3)]["reserveCosts"]
            m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-3)]["innertiaCostant"]
            m.ext[:parameters][:startupCost][i]= d[SubString(i,1:length(i)-3)]["startupCost"]
            m.ext[:parameters][:MUT][i]= d[SubString(i,1:length(i)-3)]["minUpTime"]
            m.ext[:parameters][:MDT][i]= d[SubString(i,1:length(i)-3)]["minDownTime"]
            m.ext[:parameters][:maxFrdelivarable][i]= d[SubString(i,1:length(i)-3)]["maxFrdelivarable"]
            m.ext[:parameters][:notloadCost][i]= d[SubString(i,1:length(i)-3)]["notloadCost"]
            m.ext[:parameters][:upramprate][i]= d[SubString(i,1:length(i)-3)]["upramprate"]
            m.ext[:parameters][:downramprate][i]= d[SubString(i,1:length(i)-3)]["downramprate"]

         else
            if length(i)==9
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
               m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
               m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-2)]["innertiaCostant"]
               m.ext[:parameters][:startupCost][i]= d[SubString(i,1:length(i)-2)]["startupCost"]
               m.ext[:parameters][:MUT][i]= d[SubString(i,1:length(i)-2)]["minUpTime"]
               m.ext[:parameters][:MDT][i]= d[SubString(i,1:length(i)-2)]["minDownTime"]
               m.ext[:parameters][:maxFrdelivarable][i]= d[SubString(i,1:length(i)-2)]["maxFrdelivarable"]
               m.ext[:parameters][:notloadCost][i]= d[SubString(i,1:length(i)-2)]["notloadCost"]
               m.ext[:parameters][:upramprate][i]= d[SubString(i,1:length(i)-2)]["upramprate"]
               m.ext[:parameters][:downramprate][i]= d[SubString(i,1:length(i)-2)]["downramprate"]
            else
               if length(i)==8
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-4)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-4)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-4)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-4)]["deploymentTime"] 
               m.ext[:parameters][:res_cost_g][i] =  d[SubString(i,1:length(i)-4)]["reserveCosts"]
               m.ext[:parameters][:inertia_Constant][i]= d[SubString(i,1:length(i)-4)]["innertiaCostant"]
               m.ext[:parameters][:startupCost][i]= d[SubString(i,1:length(i)-4)]["startupCost"]
               m.ext[:parameters][:MUT][i]= d[SubString(i,1:length(i)-4)]["minUpTime"]
               m.ext[:parameters][:MDT][i]= d[SubString(i,1:length(i)-4)]["minDownTime"]
               m.ext[:parameters][:maxFrdelivarable][i]= d[SubString(i,1:length(i)-4)]["maxFrdelivarable"]
               m.ext[:parameters][:notloadCost][i]= d[SubString(i,1:length(i)-4)]["notloadCost"]  
               m.ext[:parameters][:upramprate][i]= d[SubString(i,1:length(i)-4)]["upramprate"]
               m.ext[:parameters][:downramprate][i]= d[SubString(i,1:length(i)-4)]["downramprate"]	        
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
   #converted_dict_res_cost_g = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:res_cost_g]) 
   #m.ext[:parameters][:res_cost_g] = converted_dict_res_cost_g
   converted_dict_inertia_Constant = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:inertia_Constant])
   m.ext[:parameters][:inertia_Constant] = converted_dict_inertia_Constant
   converted_dict_startupCost = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:startupCost])
   m.ext[:parameters][:startupCost] = converted_dict_startupCost
   converted_dict_minUpTime = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:MUT])
   m.ext[:parameters][:MUT] = converted_dict_minUpTime
   converted_dict_minDownTime = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:MDT])
   m.ext[:parameters][:MDT] = converted_dict_minDownTime
   converted_dict_maxFrdelivarable = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:maxFrdelivarable])
   m.ext[:parameters][:maxFrdelivarable] = converted_dict_maxFrdelivarable
   converted_dict_notloadCost = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:notloadCost])
   m.ext[:parameters][:notloadCost] = converted_dict_notloadCost
   converte_dict_upramprate = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:upramprate])
   m.ext[:parameters][:upramprate] = converte_dict_upramprate
   converted_dict_downramprate = Dict{String, Int64}(k => v for (k, v) in m.ext[:parameters][:downramprate])
   m.ext[:parameters][:downramprate] = converted_dict_downramprate


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
   m.ext[:parameters][:start_up_cost_e] = Dict(i => d[SubString(i,1:length(i)-2)]["startupCost"] for i in ID_E) # Startu up cost of electrolyzers in Euros
   m.ext[:parameters][:compresor_power] = Dict(i => d[SubString(i,1:length(i)-2)]["power_compressor"] for i in ID_E) # Power consumption of compressor un MWh/kg 
   

   #Parameters Pump

   d=data
   m.ext[:parameters][:PFCOST]=Dict()
   m.ext[:parameters][:PGmaxD]=Dict()
   m.ext[:parameters][:PGminD]=Dict()
   m.ext[:parameters][:PDtg]=Dict()
   m.ext[:parameters][:Pres_cost_g]=Dict()
   m.ext[:parameters][:Pinertia_Constant]=Dict()
   m.ext[:parameters][:PstartupCost]=Dict()
   m.ext[:parameters][:PMUT]=Dict()
   m.ext[:parameters][:PMDT]=Dict()
   m.ext[:parameters][:PmaxFrdelivarable]=Dict()
   m.ext[:parameters][:PnotloadCost]=Dict()
   m.ext[:parameters][:Pupramprate]=Dict()
   m.ext[:parameters][:Pdownramprate]=Dict()
   m.ext[:parameters][:PEBmax]=Dict()
   m.ext[:parameters][:PDOD_max]=Dict()
   m.ext[:parameters][:PBeffc]=Dict()
   m.ext[:parameters][:PBeffd]=Dict()
   m.ext[:parameters][:PIni_e_b]=Dict()
   m.ext[:parameters][:PEnd_e_b]=Dict()
   m.ext[:parameters][:PDtb]=Dict()
   m.ext[:parameters][:minStableOperatingPointC]=Dict()
   m.ext[:parameters][:minStableOperatingPointD]=Dict()

   for i in ID_Pump
      m.ext[:parameters][:PFCOST][i] =  d[SubString(i,1:length(i)-2)]["marginalcost"] 
      m.ext[:parameters][:PGmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"] 
      m.ext[:parameters][:PGminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPointD"]
      m.ext[:parameters][:PDtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
      m.ext[:parameters][:Pres_cost_g][i] =  d[SubString(i,1:length(i)-2)]["reserveCosts"]
      m.ext[:parameters][:Pinertia_Constant][i]= d[SubString(i,1:length(i)-2)]["innertiaCostant"]
      m.ext[:parameters][:PstartupCost][i]= d[SubString(i,1:length(i)-2)]["startupCost"]
      m.ext[:parameters][:PMUT][i]= d[SubString(i,1:length(i)-2)]["minUpTime"]
      m.ext[:parameters][:PMDT][i]= d[SubString(i,1:length(i)-2)]["minDownTime"]
      m.ext[:parameters][:PmaxFrdelivarable][i]= d[SubString(i,1:length(i)-2)]["maxFrdelivarable"]
      m.ext[:parameters][:PnotloadCost][i]= d[SubString(i,1:length(i)-2)]["notloadCost"]
      m.ext[:parameters][:Pupramprate][i]= d[SubString(i,1:length(i)-2)]["upramprate"]
      m.ext[:parameters][:Pdownramprate][i]= d[SubString(i,1:length(i)-2)]["downramprate"]
      m.ext[:parameters][:PEBmax][i]= d[SubString(i,1:length(i)-2)]["Emax"]
      m.ext[:parameters][:PDOD_max][i]= d[SubString(i,1:length(i)-2)]["DODmax"]
      m.ext[:parameters][:PBeffc][i]= d[SubString(i,1:length(i)-2)]["effc"]
      m.ext[:parameters][:PBeffd][i]= d[SubString(i,1:length(i)-2)]["effd"]
      m.ext[:parameters][:PIni_e_b][i]= d[SubString(i,1:length(i)-2)]["Einit"]
      m.ext[:parameters][:PEnd_e_b][i]= d[SubString(i,1:length(i)-2)]["Efinal"]
      m.ext[:parameters][:PDtb][i]= d[SubString(i,1:length(i)-2)]["deploymentTime"]
      m.ext[:parameters][:minStableOperatingPointC][i]= d[SubString(i,1:length(i)-2)]["minStableOperatingPointC"]
      m.ext[:parameters][:minStableOperatingPointD][i]= d[SubString(i,1:length(i)-2)]["minStableOperatingPointD"]
   end



   #Parameter BESS
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
ID_Nuclear = m.ext[:sets][:ID_Nuclear]
ID_GR= m.ext[:sets][:ID_GR]
ID_Pump= m.ext[:sets][:ID_Pump]

# Extract time series data and convert them in PU values
D= m.ext[:timeseries][:D]
Pbase=maximum(D) 
D= m.ext[:timeseries][:D]/Pbase
#Extrac capacity factor renewable
SC= m.ext[:timeseries][:SC]
WC= m.ext[:timeseries][:WC] #data taken from https://researchdata.reading.ac.uk/191/
HVC= m.ext[:timeseries][:HVC] 
Installed_S = m.ext[:parameters][:Installed_S]/Pbase
Installed_W = m.ext[:parameters][:Installed_W]/Pbase


#Extract paremeters of the system   m.ext[:parameters][:rocofmax]=data["rocofmax"]
FO_base = m.ext[:parameters][:FO]
rocofmax = m.ext[:parameters][:rocofmax]/FO_base #rocofmax in PU
HVC=ones(24,1) #Fixed hydrogen costs
hydrogenCost=  m.ext[:parameters][:hydrogenCost]*HVC #Fixed hydrogen costs
#hydrogenCost = m.ext[:parameters][:hydrogenCost]*1.2*HVC #Variable hydrogen costs
deltaf = m.ext[:parameters][:deltaf]/FO_base
FO = m.ext[:parameters][:FO]/FO_base

# Extract parameters Generators
CostFuel=m.ext[:parameters][:FCOST]
res_cost_g=m.ext[:parameters][:res_cost_g]

GmaxD=m.ext[:parameters][:GmaxD]
GmaxD=Dict(key => value / Pbase for (key, value) in GmaxD)

GminD=m.ext[:parameters][:GminD]
GminD =Dict(key => value / Pbase for (key, value) in GminD)

startupCost=m.ext[:parameters][:startupCost]

MUT=m.ext[:parameters][:MUT]

MDT=m.ext[:parameters][:MDT]

maxFrdelivarable=m.ext[:parameters][:maxFrdelivarable]

notloadCost=m.ext[:parameters][:notloadCost]

downramprate=m.ext[:parameters][:downramprate]
downramprate=Dict(key => value / Pbase for (key, value) in downramprate)

upramprate=m.ext[:parameters][:upramprate]
upramprate=Dict(key => value / Pbase for (key, value) in upramprate)
#Dtg=m.ext[:parameters][:Dtg]
Dtg=15



#Define Mbase
Max_h_s = m.ext[:parameters][:Max_h_s]
Mbase=maximum( Max_h_s)[2]
Max_h_s =Dict(key => value /Mbase  for (key, value) in  Max_h_s)


#Extra parameters Electrolyzer

PEmax = m.ext[:parameters][:PEmax]
PEmax=Dict(key => value / Pbase for (key, value) in PEmax)

maxFrdelivarable=m.ext[:parameters][:maxFrdelivarable]
maxFrdelivarable=Dict(key => value / Pbase for (key, value) in maxFrdelivarable)

PEmin = m.ext[:parameters][:PEmin]
PEmin = Dict(key=> value / Pbase for(key,value) in PEmin)

Eeff = m.ext[:parameters][:Eeff]
Eeff =Dict(key => value*Mbase/Pbase for (key, value) in Eeff)

heffc = m.ext[:parameters][:heffc]
heffc =Dict(key => value   for (key, value) in heffc)

heffd = m.ext[:parameters][:heffd]
heffd =Dict(key => value  for (key, value) in heffd)

Eload_factor = m.ext[:parameters][:Eload_factor]
Eload_factor =Dict(key => value*0  for (key, value) in Eload_factor)

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
start_up_cost_e = m.ext[:parameters][:start_up_cost_e]
compresor_power = m.ext[:parameters][:compresor_power]
compresor_power =Dict(key => value *Mbase/Pbase  for (key, value) in compresor_power)



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

#Extract parameters Pump
CostFuel_Pump=m.ext[:parameters][:PFCOST]
res_cost_Pump=m.ext[:parameters][:Pres_cost_g]

PGmaxD=m.ext[:parameters][:PGmaxD]
PGmaxD=Dict(key => value / Pbase for (key, value) in PGmaxD)

PGminD=m.ext[:parameters][:PGminD]
PGminD =Dict(key => value / Pbase for (key, value) in PGminD)

PstartupCost=m.ext[:parameters][:PstartupCost]

PMUT=m.ext[:parameters][:PMUT]
PMDT=m.ext[:parameters][:PMDT]

PmaxFrdelivarable=m.ext[:parameters][:PmaxFrdelivarable]
PmaxFrdelivarable=Dict(key => value / Pbase for (key, value) in PmaxFrdelivarable)
PnotloadCost=m.ext[:parameters][:PnotloadCost]

Pdownramprate=m.ext[:parameters][:Pdownramprate]
Pdownramprate=Dict(key => value / Pbase for (key, value) in Pdownramprate)

Pupramprate=m.ext[:parameters][:Pupramprate]
Pupramprate=Dict(key => value / Pbase for (key, value) in Pupramprate)

PEBmax=m.ext[:parameters][:PEBmax]
PEBmax=Dict(key => value /Pbase  for (key, value) in PEBmax)

PDOD_max=m.ext[:parameters][:PDOD_max]
PBeffc=m.ext[:parameters][:PBeffc]
PBeffd=m.ext[:parameters][:PBeffd]

PIni_e_b=m.ext[:parameters][:PIni_e_b]
PIni_e_b=Dict(key => value /Pbase  for (key, value) in PIni_e_b)

PEnd_e_b=m.ext[:parameters][:PEnd_e_b]
PEnd_e_b=Dict(key => value /Pbase  for (key, value) in PEnd_e_b)

Pres_cost_g=m.ext[:parameters][:Pres_cost_g]

minStableOperatingPointC=m.ext[:parameters][:minStableOperatingPointC]
minStableOperatingPointC=Dict(key => value /Pbase  for (key, value) in minStableOperatingPointC)

minStableOperatingPointD=m.ext[:parameters][:minStableOperatingPointD]
minStableOperatingPointD=Dict(key => value /Pbase  for (key, value) in minStableOperatingPointD)

inertia_Constant=m.ext[:parameters][:inertia_Constant]
Inertia_Vector= Dict(k => inertia_Constant[k] * GmaxD[k] for k in keys(inertia_Constant) ∩ keys(GmaxD))
Pinertia_Constant=m.ext[:parameters][:Pinertia_Constant]
PInertia_Vector= Dict(k => Pinertia_Constant[k] * PGmaxD[k] for k in keys(Pinertia_Constant) ∩ keys(PGmaxD))

#Create folder to save the results

# create variables 

zuc = m.ext[:variables][:zuc] = @variable(m, [i=ID,j=J], binary=true, base_name="commitment")
v = m.ext[:variables][:v] = @variable(m, [i=ID,j=J], binary=true, base_name="start_up")
w = m.ext[:variables][:w] = @variable(m, [i=ID,j=J], binary=true, base_name="shoot_down")
g = m.ext[:variables][:g] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="generation") #Power generation generators
x = m.ext[:variables][:x] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="x") #Auxiliary variable rotate second order cone
y = m.ext[:variables][:y] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="y") #Auxiliary variable rotate second order cone
z = m.ext[:variables][:z] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="z") #Auxiliary variable rotate second order cone
xp = m.ext[:variables][:xp] = @variable(m, [i=ID_Pump,j=J],lower_bound=0, base_name="x") #Auxiliary variable rotate second order cone
yp = m.ext[:variables][:yp] = @variable(m, [i=ID_Pump,j=J],lower_bound=0, base_name="y") #Auxiliary variable rotate second order cone
zp = m.ext[:variables][:zp] = @variable(m, [i=ID_Pump,j=J],lower_bound=0, base_name="z") #Auxiliary variable rotate second order cone
rg = m.ext[:variables][:rg] = @variable(m, [i=ID,j=J],lower_bound=0,upper_bound=maxFrdelivarable[i], base_name="rg") #Reserve provided by generators
rb = m.ext[:variables][:rb] = @variable(m, [i=ID_BESS,j=J],lower_bound=0, base_name="rb") #Reserve provided by batteries
re = m.ext[:variables][:re] = @variable(m, [i=ID_E,j=J],lower_bound=0, base_name="re") #REserve provided by electrolyzers
pl = m.ext[:variables][:pl] = @variable(m, [j=J],lower_bound=0,base_name="pl") #loss of generation
RCU = m.ext[:variables][:RCU] = @variable(m, [j=J],lower_bound=0, base_name="RCU") #Renewable curtailment
pbc = m.ext[:variables][:pbc] = @variable(m, [i=ID_BESS,j=J],lower_bound=0, base_name="pbc") #Charging power of the batteries
pbd = m.ext[:variables][:pbd] = @variable(m, [i=ID_BESS,j=J],lower_bound=0, base_name="pbd") #Discharging power of the batteries
eb = m.ext[:variables][:eb] = @variable(m, [i=ID_BESS,j=J], lower_bound= EBmax[i]*(1-DOD_max[i]), upper_bound=EBmax[i] , base_name="eb") #Energy bounds of the batteries
zb = m.ext[:variables][:zb] = @variable(m, [i=ID_BESS,j=J], binary=true, base_name="on_off_b")

hfe= m.ext[:variables][:hfe] = @variable(m, [i=ID_E,j=J],lower_bound=0, upper_bound= Max_h_f[i], base_name="hfe") #Hydrogen flow limit of the hydrogen produced by electrolyzers
hfgdinyec = m.ext[:variables][:hfgdinyec] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound=Max_h_f[i],base_name="hfgdinyec") #Hydrogen flow limit of the hydrogen flowing trhow the hydrogen pipeline
hfgdcon= m.ext[:variables][:hfgdcon] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound=0,base_name="hfgdcon") #Hydrogen flow limit of the hydrogen flowing trhow the hydrogen pipeline
#zhf = m.ext[:variables][:zhf] = @variable(m, [i=ID_E,j=J], binary=true, base_name="on_off_b")
ze = m.ext[:variables][:ze] = @variable(m, [i=ID_E,j=J], binary=true, base_name="on_off_E")
zesu = m.ext[:variables][:zesu] = @variable(m, [i=ID_E,j=J], binary=true, base_name="on_off_E_startup")
zestb= m.ext[:variables][:zestb] = @variable(m, [i=ID_E,j=J], binary=true, base_name="on_off_E_stand_by")


pe = m.ext[:variables][:pe] = @variable(m,  [i=ID_E,j=J],lower_bound= 0, upper_bound=PEmax[i], base_name="pe") #Power consumption electrolyzer
pe_c= m.ext[:variables][:pe_c] = @variable(m,  [i=ID_E,j=J], base_name="pe_c") #Power consumption of compressor electrolyzer 
hss = m.ext[:variables][:hss] = @variable(m, [i=ID_E,j=J],lower_bound=Min_h_s[i], upper_bound= Max_h_s[i], base_name="hss") #hydrogen storage limit


#Pump variables

Ppc= m.ext[:variables][:Ppc] = @variable(m, [i=ID_Pump,j=J],lower_bound=0,upper_bound=PGmaxD[i],base_name="Ppc") #Charging variable

Ppd= m.ext[:variables][:Ppd] = @variable(m, [i=ID_Pump,j=J],lower_bound=0,upper_bound=PGmaxD[i],base_name="Ppd") #Discharging variable

Pener= m.ext[:variables][:Pener] = @variable(m, [i=ID_Pump,j=J],lower_bound=PEBmax[i]*(1-PDOD_max[i]),upper_bound=PEBmax[i],base_name="Pener") #Energy variable

zpcommit= m.ext[:variables][:zpcommit] = @variable(m, [i=ID_Pump,j=J], binary=true, base_name="zp") #Charging variable

rgp= m.ext[:variables][:rgp] = @variable(m, [i=ID_Pump,j=J],lower_bound=0,upper_bound=PmaxFrdelivarable[i], base_name="rgp") #Reserve provided by generators


#create affine expressions

g_costs=m.ext[:expressions][:g_costs] = @expression(m, [i=ID,j=J],g[i,j]*CostFuel[i]*Pbase
)
rg_costs=m.ext[:expressions][:rg_costs] = @expression(m, [i=ID,j=J],rg[i,j]*res_cost_g[i]*Pbase
   )

rgp_costs=m.ext[:expressions][:prg_costs] = @expression(m, [i=ID_Pump,j=J],rgp[i,j]*Pres_cost_g[i]*Pbase
   )
rb_costs=m.ext[:expressions][:rb_costs] = @expression(m, [i=ID_BESS,j=J],rb[i,j]*res_cost_b[i]*Pbase
   )
re_costs=m.ext[:expressions][:re_costs] = @expression(m, [i=ID_E,j=J],re[i,j]*res_cost_e[i]*Pbase
   )
h_costs=m.ext[:expressions][:h_costs] = @expression(m, [i=ID_E,j=J],(hfgdcon[i,j]-hfgdinyec[i,j])*hydrogenCost[j]*Mbase
   )
scu_UC=m.ext[:expressions][:scu_UC] = @expression(m, [i=ID,j=J], startupCost[i]*v[i,j])

scu_UC_e=m.ext[:expressions][:scu_UC_e] = @expression(m, [i=ID_E,j=J], start_up_cost_e[i]*zesu[i,j])


#Create folder to save the results
RE_costs=res_cost_e["E_500_1"]
RG_costs=res_cost_g["CCGT_77"]
Installed_W_F=Installed_W*Pbase
Installed_S_F=Installed_S*Pbase
folder_name_plot="UC_CRE_$(RE_costs)_CRG_$(RG_costs)_IW_$(Installed_W_F)_IS_$(Installed_S_F)_S3"
mkdir(folder_name_plot)


#It creates the expression of the inertia
# Combine index sets (assuming ID and ID_Pump are disjoint)

Combined_ID = vcat(collect(ID), collect(ID_Pump)) 
Inertia_Expression=m.ext[:expressions][:Inertia_Expression] = @expression(m, [i in Combined_ID, j in J],
    i in ID ? Inertia_Vector[i] * zuc[i,j] : PInertia_Vector[i]
)


# Inertia expression without considering pump
#Inertia_Expression=m.ext[:expressions][:Inertia_Expression] = @expression(m, [i=ID,j=J],Inertia_Vector[i]*zuc[i,j])

#Create the objective function
obj= m.ext[:objective] = @objective(m,Min, (sum(g_costs)+sum(rg_costs)+sum(rgp_costs)+sum(rb_costs)+sum(re_costs)+sum(scu_UC)+sum(h_costs)+sum(scu_UC_e)) #Objective function
)

con1=m.ext[:constraints][:con1] = @constraint(m, [j=J],
WC[j]*Installed_W + SC[j]*Installed_S-RCU[j] + sum(g[i,j] for i in ID) - sum(pbc[i,j] for i in ID_BESS) + sum(pbd[i,j] for i in ID_BESS)-sum(Ppc[i,j] for i in ID_Pump) +sum(Ppd[i,j] for i in ID_Pump) == D[j] + sum(pe[i,j] for i in ID_E) + sum(pe_c[i,j] for i in ID_E)
)

#constraint renewable courtailment
con1_1=m.ext[:constraints][:con1_1] = @constraint(m, [j=J],
RCU[j] <= WC[j]*Installed_W + SC[j]*Installed_S
)

#Constraignt upper bound reserve generators
con1_2=m.ext[:constraints][:con1_2] = @constraint(m, [i=ID,j=J],
rg[i,j].<=maxFrdelivarable[i]*zuc[i,j]
)
#Constraignt upper bound generators power
con2_1=m.ext[:constraints][:con2_1] = @constraint(m, [i=ID,j=J],
g[i,j]+rg[i,j].<=GmaxD[i]*zuc[i,j]
)
#Constraint upper bound generators power
con2_2=m.ext[:constraints][:con2_2] = @constraint(m, [i=ID,j=J],
GminD[i]*zuc[i,j].<=g[i,j]+rg[i,j]
)


#Constraint upper bound reserve pump
con1_2_Pump=m.ext[:constraints][:con1_2_Pump] = @constraint(m, [i=ID_Pump,j=J],
rgp[i,j].<=PmaxFrdelivarable[i]*zpcommit[i,j]
)
#Constraint upper bound pump discharging mode
con2_1_Pump=m.ext[:constraints][:con2_1_Pump] = @constraint(m, [i=ID_Pump,j=J],
Ppd[i,j]+rgp[i,j].<=PGmaxD[i]*zpcommit[i,j]
)


#Constraint lower bound pump discharging mode
con2_1_1_Pump=m.ext[:constraints][:con2_1_1_Pump] = @constraint(m, [i=ID_Pump,j=J],
minStableOperatingPointD[i]*zpcommit[i,j].<=Ppd[i,j]+rgp[i,j]
)


#Constraint upper bound pump charging mode
con2_2_Pump=m.ext[:constraints][:con2_2_Pump] = @constraint(m, [i=ID_Pump,j=J],
Ppc[i,j].<=PGmaxD[i]*(1-zpcommit[i,j])
)

#constraint lower bound pump charging mode
con2_2_1_Pump=m.ext[:constraints][:con2_2_1_Pump] = @constraint(m, [i=ID_Pump,j=J],
Ppc[i,j].>=minStableOperatingPointC[i]*(1-zpcommit[i,j])
)



#Constraint charging and discharging mode of the batteries
con2_1_1=m.ext[:constraints][:con2_1_1] = @constraint(m, [i=ID_BESS,j=J],
pbc[i,j].<=PBmax[i]*(1-zb[i,j])
)
con2_1_2=m.ext[:constraints][:con2_1_2] = @constraint(m, [i=ID_BESS,j=J],
pbd[i,j].<=PBmax[i]*zb[i,j]
)


con2_2_a = m.ext[:constraints][:con2_2_a] = @constraint(m, [i=ID,j=J[2:end]],
g[i,j] +rg[i,j]- g[i,j-1] <= upramprate[i] ) # constraint taken from Mathematical Programming for Power Systems Operation

con2_2_b = m.ext[:constraints][:con2_2_b] = @constraint(m, [i=ID,j=J[2:end]],
g[i,j-1] - g[i,j] <= downramprate[i] )


#Initial status generators (unit commitment)
con2_3=m.ext[:constraints][:con2_3] = @constraint(m, [i=ID,j=J[1]],1-zuc[i,j]+v[i,j]-w[i,j]==0)

#Status generators (unit commitment)
con2_4=m.ext[:constraints][:con2_4] = @constraint(m, [i=ID,j=J[2:end]],
zuc[i,j-1]-zuc[i,j]+v[i,j]-w[i,j]==0
)

con_2_4_1=m.ext[:constraints][:con_2_4_1]=@constraint(m, [i=ID,j=J],+v[i,j]+w[i,j]<=1
)


#On status constraint for nuclear
con2_5=m.ext[:constraints][:con2_5] = @constraint(m, [i=ID_Nuclear,j=J],
zuc[i,j]==1
)






#Minimum Up Time constraint
con2_6 = m.ext[:constraints][:con2_6] = Dict()
for i in ID
   if MDT[i] > 0  # Only proceed if MUT[i] is positive, to avoid issues with MDT=0
       for j in J[MUT[i]:end]
           con2_6[i,j] =  @constraint(m,
           zuc[i,j] >= sum(v[i,j-jj] for jj in 0:MUT[i]-1)
           )
       end
   end
end


#Minimum Down Time constraint
con2_7 = m.ext[:constraints][:con2_7] = Dict()
for i in ID
   if MDT[i] > 0  # Only proceed if MDT[i] is positive
       for j in J[MDT[i]:end]
           con2_7[i,j] = @constraint(m,
               1 - zuc[i,j] >= sum(w[i,j-jj] for jj in 0:MDT[i]-1)
           )
       end
   end
end





#Constraint lost of generation

con3=m.ext[:constraints][:con3] = @constraint(m, [i=ID,j=J],
pl[j].>=g[i,j]
)
con3_Pump=m.ext[:constraints][:con3] = @constraint(m, [i=ID_Pump,j=J],
pl[j].>=Ppd[i,j]
)

#= Constraints that include the loss of wind and solar generation
con3_1=m.ext[:constraints][:con3_1] = @constraint(m, [i=ID,j=J],
pl[j].>=SC[j]*Installed_S
)

con3_2=m.ext[:constraints][:con3_1] = @constraint(m, [i=ID,j=J],
pl[j].>=WC[j]*Installed_W*0
)
=#

#Constraint upper bound reserve provided by BESS
con4=m.ext[:constraints][:con4] = @constraint(m, [i=ID_BESS,j=J],
rb[i,j].<=PBmax[i]+pbc[i,j]-pbd[i,j]
)

#Constraint upper bound reserve provided by Electrolyzer
con5=m.ext[:constraints][:con5] = @constraint(m, [i=ID_E,j=J],re[i,j].<=pe[i,j]-PEmin[i]*ze[i,j] 
)



#Constraint end energy value of the batteries
con6=m.ext[:constraints][:con6] = @constraint(m, [i=ID_BESS,j=J[end]],End_e_b[i]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])
#con6=m.ext[:constraints][:con6] = @constraint(m, [i=ID_BESS,j=J[end]],End_e_b[i]==eb[i,j])


#Constraint initial value energy of the batteries
#con7=m.ext[:constraints][:con7] = @constraint(m, [i=ID_BESS,j=J[1]],eb[i,j+1]-Ini_e_b[i]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])
con7=m.ext[:constraints][:con7] = @constraint(m, [i=ID_BESS,j=J[1]],eb[i,j]==Ini_e_b[i])

#Constraint charging-discharging batteries
con8=m.ext[:constraints][:con8] = @constraint(m, [i=ID_BESS,j=J[1:end-1]],eb[i,j+1]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])

#Constraint charging-discharging pump
con6_Pump=m.ext[:constraints][:con6_Pump] = @constraint(m, [i=ID_Pump,j=J[end]],PEnd_e_b[i]-Pener[i,j]==Ppc[i,j]*PBeffc[i]-Ppd[i,j]/PBeffd[i])
con7_Pump=m.ext[:constraints][:con7_Pump] = @constraint(m, [i=ID_Pump,j=J[1]],Pener[i,j]==PIni_e_b[i])
con8_Pump=m.ext[:constraints][:con8_Pump] = @constraint(m, [i=ID_Pump,j=J[1:end-1]],Pener[i,j+1]-Pener[i,j]==Ppc[i,j]*PBeffc[i]-Ppd[i,j]/PBeffd[i])
#


##Constraint Electrolyzer
#Constraint  hydrogen production electrolyzer considering standby
con9=m.ext[:constraints][:con9] = @constraint(m, [i=ID_E,j=J],hfe[i,j]==pe[i,j]/Eeff[i]-0.05*PEmax[i]/Eeff[i]*zestb[i,j]
)

#Constraint  hydrogen production electrolyzer without considering standby
#con9=m.ext[:constraints][:con9] = @constraint(m, [i=ID_E,j=J],hfe[i,j]==pe[i,j]/Eeff[i])


#Hydrogen storage constraints




#Constraint end hydrogen value of the hydrogen storage
con11=m.ext[:constraints][:con11] = @constraint(m, [i=ID_E,j=J[end]],End_h_s[i]==hss[i,j]+hfe[i,j]-hfgdinyec[i,j]/heffc[i]+hfgdcon[i,j]*heffd[i]-Eload_factor[i]*PEmax[i]/(Eeff[i]))
#con11=m.ext[:constraints][:con11] = @constraint(m, [i=ID_E,j=J[end]],End_h_s[i]==hss[i,j])


#Constraint initial value of the hydrogen storage
con12=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[1]],hss[i,j+1]==Ini_h_s[i]+hfe[i,j]-hfgdinyec[i,j]/heffc[i]+hfgdcon[i,j]*heffd[i]-Eload_factor[i]*PEmax[i]/(Eeff[i]))
#con12=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[1]],hss[i,j]==Ini_h_s[i])

#Constraint charging-discharging of the hydrogen storage
con13=m.ext[:constraints][:con13] = @constraint(m, [i=ID_E,j=J[1:end-1]],hss[i,j+1]==hss[i,j]+hfe[i,j]-hfgdinyec[i,j]/heffc[i]+hfgdcon[i,j]*heffd[i]-Eload_factor[i]*PEmax[i]/(Eeff[i]))

#Bounds electrolyzer considering standby
con13_1=m.ext[:constraints][:con13_1] = @constraint(m, [i=ID_E,j=J],pe[i,j].<=PEmax[i]*ze[i,j]+0.05*PEmax[i]*zestb[i,j])
con13_2=m.ext[:constraints][:con13_2] = @constraint(m, [i=ID_E,j=J],pe[i,j].>=PEmin[i]*ze[i,j]+0.05*PEmax[i]*zestb[i,j])

#Bounds electrolyzer without considering standby
#con13_1=m.ext[:constraints][:con13_1] = @constraint(m, [i=ID_E,j=J],pe[i,j].<=PEmax[i]*ze[i,j])
#con13_2=m.ext[:constraints][:con13_2] = @constraint(m, [i=ID_E,j=J],pe[i,j].>=PEmin[i]*ze[i,j])

#Status constraint of electrolyzers
con13_2_1=m.ext[:constraints][:con13_2_1] = @constraint(m, [i=ID_E,j=J[1]],zesu[i,j].>=ze[i,j]-1+(zestb[i,j]-0)) 

#Startup constraint of electrolyzers
con13_2_2=m.ext[:constraints][:con13_2_2] = @constraint(m, [i=ID_E,j=J[2:end]],zesu[i,j].>=(ze[i,j]-ze[i,j-1])+(zestb[i,j]-zestb[i,j-1])) 

#Startup constraint of electrolyzers
con13_2_3=m.ext[:constraints][:con13_2_3] = @constraint(m, [i=ID_E,j=J],zestb[i,j]+ze[i,j].<=1
) #Standby constraint of electrolyzers

#con13_2_4=m.ext[:constraints][:con13_2_4] = @constraint(m, [i=ID_E,j=J], ze[i,j]==0)
#Power consumption in standby=5% of nominal power according to chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2306.10962


#Status constraints electrolyzer


#Power consumption of the electrolyzer
con13_8=m.ext[:constraints][:con13_8] = @constraint(m, [i=ID_E,j=J],pe_c[i,j]==compresor_power[i]*(hfe[i,j]+hfgdcon[i,j]+hfgdinyec[i,j]))


#Constraint maximum frequency variation
con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO*(sum(rb[:, j])/Dtb+sum(re[:, j])/Dte+(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j])/Dtg)).>=pl[j]^2/(4*deltaf))
con14_pump= m.ext[:constraints][:con14_pump] = @constraint(m, [i=ID_Pump,j=J],((sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO*(sum(rb[:, j])/Dtb+sum(re[:, j])/Dte+(sum(rg[:, j]))/Dtg)).>=pl[j]^2/(4*deltaf))

#constraints nadir occurrence time

con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=0)

con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<= 0.00001+sum(re[:, j])+sum(rb[:,j])*(Dte/Dtb)+(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j])*Dte/Dtg)
con16_pump=m.ext[:constraints][:con16_pump] = @constraint(m, [i=ID_Pump,j=J],pl[j].<= 0.00001+sum(re[:, j])+sum(rb[:,j])*(Dte/Dtb)+(sum(rg[:, j]))*Dte/Dtg)

#constraint ROCOF

con17=m.ext[:constraints][:con17] = @constraint(m, [i=ID,j=J], pl[j]*FO/2 .<=rocofmax*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j]))
con17_pump=m.ext[:constraints][:con17_pump] = @constraint(m, [i=ID_Pump,j=J], pl[j]*FO/2 .<=rocofmax*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j]))

#QSS frequency constraint
con18=m.ext[:constraints][:con18] = @constraint(m, [i=ID,j=J],pl[j].<=0.00001+sum(re[:, j])+sum(rb[:, j])+(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j]))
con18_pump=m.ext[:constraints][:con18_pump] = @constraint(m, [i=ID_Pump,j=J],pl[j].<=0.00001+sum(re[:, j])+sum(rb[:, j])+(sum(rg[:, j])))




Model_1_time=@elapsed begin
# Specify the full path for the output file
output_file_path = joinpath(folder_name_plot, "output_1.txt")
# Open the file with the full path and write the output
open(output_file_path, "w") do file
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
      #set_optimizer_attribute(m, "Threads", 4)
       optimize!(m)
       opt_time = solve_time(m)
   end
end
end



#=
#Constraints nadir interval II
for j in J
    for i in ID
       delete(m,m.ext[:constraints][:con14][i,j])
       delete(m,m.ext[:constraints][:con15][i,j])
       delete(m,m.ext[:constraints][:con16][i,j])
    end
end


   #Constraints frequency nadir as rotate second order cone
   #con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO-sum(re[:,j])*Dte/(4*deltaf))*(sum(rb[:, j])/Dtb+(sum(rg[:, j])-rg[i, j])/Dtg).>=(pl[j]-sum(re[:, j]))^2/(4*deltaf))
   
   con14_1=m.ext[:constraints][:con14_1] = @constraint(m,[i=ID,j=J],y[i,j] ==2*deltaf*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO-sum(re[:, j])*Dte/2)
   con14_2=m.ext[:constraints][:con14_2] = @constraint(m,[i=ID,j=J],z[i,j] ==sum(re[:, j])/Dte-(sum(rg[:, j])-rg[i, j])/Dtg)
   con14_3=m.ext[:constraints][:con14_3] = @constraint(m,[i=ID,j=J],[y[i,j]; z[i,j]; pl[j]-sum(re[:, j])] in RotatedSecondOrderCone())
   #constraints nadir occurrence time
   con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(re[:, j])  + sum(rb[:, j])*Dte/Dtb + (sum(rg[:, j]) - rg[i, j]) * Dte / Dtg)
   con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000000000000001 + sum(re[:, j]) + sum(rb[:, j]) + (sum(rg[:, j]) - rg[i, j]) * Dtb / Dtg)

   output_file_path = joinpath(folder_name_plot, "output_2.txt")
# Open the file with the full path and write the output
open(output_file_path, "w") do file
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
       optimize!(m)
       opt_time = solve_time(m)
   end
end
 =#
 
Model_3_time=@elapsed begin
 #Constraints nadir interval III
for j in J
    for i in ID
      #=
      #When running with 3 nadir intervals
       delete(m,m.ext[:constraints][:con14_1][i,j])
       delete(m,m.ext[:constraints][:con14_2][i,j])
       delete(m,m.ext[:constraints][:con14_3][i,j]) 
       delete(m,m.ext[:constraints][:con15][i,j])
       delete(m,m.ext[:constraints][:con16][i,j])
      =#
      #when running with 2 nadir intervals
       delete(m,m.ext[:constraints][:con14][i,j])
       delete(m,m.ext[:constraints][:con15][i,j])
       delete(m,m.ext[:constraints][:con16][i,j])
      
    end

    for i in ID_Pump
      delete(m,m.ext[:constraints][:con14_pump][i,j])
      delete(m,m.ext[:constraints][:con16_pump][i,j])
    end
 end



#Constraints frequency nadir as rotate second order cone
 con14_1=m.ext[:constraints][:con14_1] = @constraint(m,[i=ID,j=J],y[i,j] ==2*deltaf*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO-sum(re[:, j])*Dte/2-sum(rb[:, j])*Dtb/2)
 con14_2=m.ext[:constraints][:con14_2] = @constraint(m,[i=ID,j=J],z[i,j] ==(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j])/Dtg)
 con14_3=m.ext[:constraints][:con14_3] = @constraint(m,[i=ID,j=J],[y[i,j]; z[i,j]; pl[j]-sum(re[:, j])-sum(rb[:,j])] in RotatedSecondOrderCone())



 con14_1_Pump=m.ext[:constraints][:con14_1_Pump] = @constraint(m,[i=ID_Pump,j=1],yp[i,j] ==2*deltaf*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j]*zpcommit[i,j])/FO-sum(re[:, j])*Dte/2-sum(rb[:, j])*Dtb/2)
 con14_2_Pump=m.ext[:constraints][:con14_2_Pump] = @constraint(m,[i=ID_Pump,j=1],zp[i,j] ==(sum(rg[:, j])+sum(rgp[:,j])-rgp[i, j])/Dtg)
 con14_3_Pump=m.ext[:constraints][:con14_3_Pump] = @constraint(m,[i=ID_Pump,j=1],[yp[i,j]; zp[i,j]; pl[j]-sum(re[:, j])-sum(rb[:,j])] in RotatedSecondOrderCone())
 #constraints nadir occurrence time
 con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j])*Dtb/Dtg)
 con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000001 + sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])+sum(rgp[:,j])-rg[i, j]))

   con15_pump= m.ext[:constraints][:con15_pump] = @constraint(m, [i=ID_Pump,j=J],pl[j].>=sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])+sum(rgp[:,j])-rgp[i,j])*Dtb/Dtg)
   con16_pump=m.ext[:constraints][:con16_pump] = @constraint(m, [i=ID_Pump,j=J],pl[j].<=0.000001 + sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])+sum(rgp[:,j])-rgp[i,j]))



 output_file_path = joinpath(folder_name_plot, "output_3.txt")
# Open the file with the full path and write the output
open(output_file_path, "w") do file 
   # Redirect stdout to the file within the block
   redirect_stdout(file) do
      #set_optimizer_attribute(m, "Threads", 4)
        set_optimizer(m, Gurobi.Optimizer)
        optimize!(m)
        opt_time = solve_time(m)
   end
end

end






Post_Processing_time = @elapsed begin


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
hfgdinyec= value.(m.ext[:variables][:hfgdinyec])*Mbase
hfgdcon= value.(m.ext[:variables][:hfgdcon])*Mbase
hss= value.(m.ext[:variables][:hss])*Mbase
zucvalues=  value.(m.ext[:variables][:zuc])
pe_c= value.(m.ext[:variables][:pe_c])*Pbase
wvalues=  value.(m.ext[:variables][:w])
vvalues=  value.(m.ext[:variables][:v])
RCU= value.(m.ext[:variables][:RCU])*Pbase
ze= value.(m.ext[:variables][:ze])
zesu= value.(m.ext[:variables][:zesu])
zestb= value.(m.ext[:variables][:zestb])
zb= value.(m.ext[:variables][:zb])
Ppd= value.(m.ext[:variables][:Ppd])*Pbase
Ppc= value.(m.ext[:variables][:Ppc])*Pbase
Pener= value.(m.ext[:variables][:Pener])*Pbase
zpcommit= value.(m.ext[:variables][:zpcommit])
rgp= value.(m.ext[:variables][:rgp])*Pbase



HD = Dict(key => Mbase*Eload_factor[key] * PEmax[key] / Eeff[key] for key in keys(Eload_factor))

Ivec = [i for  i in I]
Inuc = [i for  i in ID_Nuclear]
Bvec = [i for  i in ID_BESS]
Evec = [i for  i in ID_E]
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
hfgdinyecvec = [hfgdinyec[i,j] for  i in ID_E, j in J]
hfgdconvec = [hfgdcon[i,j] for  i in ID_E, j in J]
hssvec = [hss[i,j] for  i in ID_E, j in J]
HDvec = [HD[i] for  i in ID_E, j in J]
plvec= [pl[j] for j in J]
RCUvec = [RCU[j] for j in J]
zucvector = [zucvalues[i,j] for i in ID, j in J]
wvector = [wvalues[i,j] for i in ID, j in J]
vvector = [vvalues[i,j] for i in ID, j in J]
zevector= [ze[i,j] for i in ID_E, j in J]
zesuvector=[zesu[i,j] for i in ID_E, j in J]
zestbvector=[zestb[i,j] for i in ID_E, j in J]
zpcommitvector= [zpcommit[i,j] for i in ID_Pump, j in J]
Ppdvector= [Ppd[i,j] for i in ID_Pump, j in J]
Ppcvector= [Ppc[i,j] for i in ID_Pump, j in J]
Penervector= [Pener[i,j] for i in ID_Pump, j in J]
rgpvector= [rgp[i,j] for i in ID_Pump, j in J]
ps=SC*Installed_S*Pbase
pw=WC*Installed_W*Pbase
pe_cvec= [pe_c[i,j] for  i in ID_E, j in J]




using StatsPlots
save_path = joinpath(folder_name_plot, "Wind_and_curtailment.png")
P_curtailment = plot(RCUvec, label = "RCU", lw = 2, color = :red,title = "Wind generation and renewable curtailment")
plot!(pw, label = "Wind power", lw = 2, color = :blue)
savefig(P_curtailment, save_path)
display(P_curtailment)


p_generators = groupedbar(transpose(gvec[:,:]),
    bar_position = :stack,
    label = permutedims(ID),
    legend = :outertopright,
    xlabel = "Time [h]",
    ylabel = "Power [MW]",
    title = "Generated Power of Units")
plot(p_generators, layout = (1,2), size=(500, 500))

p_RB = groupedbar(transpose(rbvec[:,:]),
    bar_position = :stack,
    label = permutedims(Bvec),
    legend = :outertopright,
    xlabel = "Time [h]",
    ylabel = "Procured reserve [MW]",
    title = "Hourly procured reserve BESS");
plot(p_RB, layout = (1,2), size=(500, 500))

p_RG = groupedbar(transpose(rgvec[:,:]),
    bar_position = :stack,
    label = permutedims(Ivec),
    legend = :outertopright,
    xlabel = "Time [h]",
    ylabel = "Procured reserve [MW]",
    title = "Hourly procured reserve generators");
plot(p_RG, layout = (1,2), size=(500, 500))

p_re = groupedbar(transpose(revec[:,:]),
    bar_position = :stack,
    label = permutedims(Evec),
    legend = :outertopright,
    xlabel = "Time [h]",
    ylabel = "Procured reserve [MW]",
    title = "Hourly procured reserve Electrolyzers");
plot(p_re, layout = (1,2), size=(500, 500))



p_h_storage = plot(hfgdinyecvec[1,:], 
     legend = :outerright,
     linewidth = 2,
     color = :blue,
     title = "Stored hydrogen and hydrogen storage flows", 
     xlabel = "Time [H]", 
     ylabel = "Mass flow [Kg/h], Mass [Kg]", 
     label = "hfgdinyecvec",
     ylims = (0, maximum([maximum(hfgdinyecvec[1,:]), maximum(hssvec[1,:]), maximum(hfgdconvec[1,:])]) + 10))  # Added ylims

plot!(hssvec[1,:],
      linewidth = 2,
      color = :red,
      label = "hssvec")

plot!(hfgdconvec[1,:], 
      linewidth = 2,
      color = :orange,
      label = "hfgdconvec")



display(p_h_storage)

save_path = joinpath(folder_name_plot, "Hydrogen_storage_flows.png")
savefig(p_h_storage, save_path)

# Plot the vectors
p_HF = plot(hfevec[1,:], label = "hfevec", lw = 2)
plot!(hfgdinyecvec[1,:], label = "hfgdinyecvec", lw = 2)
plot!(hfgdconvec[1,:], label = "hfgdconvec", lw = 2)
plot!(hssvec[1,:], label = "hssvec", lw = 2)
plot!(HDvec[1, :], 
color = :red, 
seriestype = :scatter, 
markershape = :circle,
label = "HD")

plot!(xlabel = "Time [h]",
      ylabel = "Hydrogen mass [kg]",
      legend = :outertopright,
      title = "Hourly Hydrogen mass flow")

display(p_HF)

for i in 1:20
p_HF_11 = plot(hfevec[i,:], label = "hfevec", lw = 2)
plot!(hfgdinyecvec[i,:], label = "hfgdinyecvec", lw = 2)
plot!(hfgdconvec[i,:], label = "hfgdconvec", lw = 2)
plot!(hssvec[i,:], label = "hssvec", lw = 2)
plot!(HDvec[i, :], 
color = :red, 
seriestype = :scatter, 
markershape = :circle,
label = "HD")

plot!(xlabel = "Time [h]",
      ylabel = "Hydrogen mass [kg] ",
      legend = :outertopright,
      title = "Hourly Hydrogen mass flow electrolyzer $i")

display(p_HF_11)

save_path = joinpath(folder_name_plot, "Hydrogen_flows_plot_$i.png")
savefig(p_HF_11, save_path)
end




save_path = joinpath(folder_name_plot, "Produced_power_plot.png")
savefig(p_generators, save_path)

save_path = joinpath(folder_name_plot, "BESS_reserve_plot.png")
savefig(p_RB, save_path)

save_path = joinpath(folder_name_plot, "Electrolyzer_Freserve_plot.png")
savefig(p_re, save_path)

save_path = joinpath(folder_name_plot, "Generator_r_plot.png")
savefig(p_RG, save_path)

save_path = joinpath(folder_name_plot, "Hydrogen_flows_plot.png")
savefig(p_HF, save_path)


sum_reserve_g = vec(sum(rgvec, dims=1))
sum_reserve_b = vec(sum(rbvec, dims=1))
sum_reserve_e = vec(sum(revec, dims=1))
total_reserve = sum_reserve_g + sum_reserve_b + sum_reserve_e

#It creates a folder to store the results

number_hours=nrow(ts)



pr=bar(1:number_hours, [sum_reserve_g sum_reserve_b sum_reserve_e],
    bar_width = 0.8,           # Ancho de las barras
    label = ["Reserve Generators" "Reserve BESS" "Reserve Electrolyzer"],  # Etiquetas para la leyenda
    fillcolor = [:red :green :blue],
    title = "Hourly Procured reserve",
    xlabel = "Time [h]",       # Etiqueta del eje X
    ylabel = "Procured Reserve [MW]",
    alpha = 0.5,               # Transparencia de las barras
    layout = (1, 1),
   legend = :outerright
    )           # Disposición del gráfico


plot!(pr, 1:number_hours, plvec,
    label = "Loss of power",
    linewidth = 2,
    color = :purple,
    linestyle = :solid)

plot!(pr, 1:number_hours, total_reserve,
label = "Total reserve",
linewidth = 2,
color = :red,
linestyle = :dash)

display(pr)
save_path = joinpath(folder_name_plot, "reserves_plot.png")
savefig(pr, save_path)




p_bat = plot(pbcvec[1,:], 
     legend = :outerright,
     linewidth = 2,
     color = :blue,           # Changed from purple to blue for clarity
     title = "Power and Energy of the BESS", 
     xlabel = "Time [H]", 
     ylabel = "Power [MW], Energy [MWh]", 
     label = "Charging power")

plot!(ebvec[1,:],
      linewidth = 2,
      color = :red,         # Changed from red to orange for better contrast
      label = "Battery energy" # Fixed typo in "energyr"
)

# Add pbdvec[1,:] to the same plot
plot!(pbdvec[1,:], 
      linewidth = 2,
      color = :orange,      # Changed from green to darkgreen for distinction
      label = "Discharging power")

# Display the combined plot
display(p_bat)

save_path = joinpath(folder_name_plot, "BESS_power_energy_plot.png")
savefig(p_bat, save_path)



# plot pump values

p_pump = plot(Ppcvector[1,:], 
     legend = :outerright,
     linewidth = 2,
     color = :blue,           # Changed from purple to blue for clarity
     title = "Power and Energy of the Pump", 
     xlabel = "Time [H]", 
     ylabel = "Power [MW], Energy [MWh]", 
     label = "Charging power")

     plot!(Penervector[1,:],
      linewidth = 2,
      color = :red,         # Changed from red to orange for better contrast
      label = "Pump energy") # Fixed typo in "energy"

      plot!(Ppdvector[1,:],
      linewidth = 2,
      color = :orange,      # Changed from green to darkgreen for distinction
      label = "Discharging power")
   display(p_pump)
save_path = joinpath(folder_name_plot, "Pump_power_energy_plot.png")
savefig(p_pump, save_path)


Sum_Inertia_Vector=Dict()

for j in J
   Sum_Inertia_Vector[j]=sum(Inertia_Vector[i]*value.(zuc[i,j]) for i in ID)+ PInertia_Vector["Pump_1"]*value.(zpcommit["Pump_1",j])
end



x_Inertia = collect(keys(Sum_Inertia_Vector))
y_Inertia = collect(values(Sum_Inertia_Vector))

# Create scatter plot
p_inertia=scatter(x_Inertia, y_Inertia, 
    title = "Total Inertia vs Time",
    xlabel = "Time[h]",
    ylabel = "Total sum of system's inertia [s]",
    legend = false,
    markersize = 5,
    ylims = (0, maximum(y_Inertia)+3)
)
save_path = joinpath(folder_name_plot, "inertia_plot.png")
savefig(p_inertia, save_path)

Sum_Inertia_Vector_energy=Dict()

for j in J
   Sum_Inertia_Vector_energy[j]=sum(Inertia_Vector[i]*value.(zuc[i,j])*Pbase/1000 for i in ID)+ PInertia_Vector["Pump_1"]*value.(zpcommit["Pump_1",j])*Pbase/1000
end



x_Inertia_energy = collect(keys(Sum_Inertia_Vector_energy))
y_Inertia_energy = collect(values(Sum_Inertia_Vector_energy))

# Create scatter plot
p_inertia_energy=scatter(x_Inertia_energy, y_Inertia_energy, 
    title = "Total Inertia vs Time",
    xlabel = "Time[h]",
    ylabel = "Total sum of system's inertia [GWs]",
    legend = false,
    markersize = 5,
    ylims = (0, maximum(y_Inertia_energy)+100)
)
save_path = joinpath(folder_name_plot, "inertia_plot_energy.png")
savefig(p_inertia_energy, save_path)





P_D_W_S_N=plot(D*Pbase, 
      linewidth = 2,
      color = :red,
     title = "Demand, wind and solar generation", 
     xlabel = "Time [H]", 
     ylabel = "Power [MW]", 
     label = "Demand")

plot!(ps,
linewidth = 2,
color = :green,
label = "Solar"
)

# Add pbdvec[1,:] to the same plot with a different cross style
plot!(pw, 
linewidth = 2,
color = :blue,
label = "Wind"
)

plot!(D*Pbase-ps-pw, 
linewidth = 2,
color = :orange,
label = "Net demand"
)
# Save the combined plot

display(P_D_W_S_N)

save_path = joinpath(folder_name_plot, "Demand_renewables.png")
savefig(P_D_W_S_N, save_path)

Total_Demand_Electro_storage=(D'*Pbase+sum(pevec, dims=1)+sum(pbcvec, dims=1)+sum(Ppcvector, dims=1)+sum(pe_cvec, dims=1))'

P_D_W_S_N_all=plot(Total_Demand_Electro_storage, 
      linewidth = 2,
      color = :red,
     title = "Demand_all, wind and solar generation", 
     xlabel = "Time [H]", 
     ylabel = "Power [MW]", 
     label = "Demand")

plot!(ps,
linewidth = 2,
color = :green,
label = "Solar"
)

# Add pbdvec[1,:] to the same plot with a different cross style
plot!(pw, 
linewidth = 2,
color = :blue,
label = "Wind"
)

plot!(Total_Demand_Electro_storage-ps-pw, 
linewidth = 2,
color = :orange,
label = "Net demand all"
)
# Save the combined plot

display(P_D_W_S_N_all)

save_path = joinpath(folder_name_plot, "Total_Demand_renewables.png")
savefig(P_D_W_S_N_all, save_path)

end

OU_v=sum(zucvector', dims=2)
online_units=bar(OU_v, 
    xlabel="Hours",                     # Eje X
    ylabel="Number of online units",    # Eje Y
    title="Online units at each hour",  # Título
    legend=false,                       # Sin leyenda
    bar_width=0.8,                     # Ancho de las barras
    size=(500, 300))                    # Tamaño del gráfico
display(online_units)
# Guardar el gráfico (opcional)
save_path = joinpath(folder_name_plot, "Online_units.png")
savefig(online_units, save_path)


#Making file with the computing time
output_file_path = joinpath(folder_name_plot, "Computing_time.txt")
# Ensure the folder exists (optional, but recommended)
mkpath(folder_name_plot)
# Use output_file_path in open

open(output_file_path, "w") do file
    write(file, "Computing time model 1= $(string(Model_1_time)), Computing time model 3= $(string(Model_3_time)), Post processing time= $(string(Post_Processing_time))")  # Convert to string and write
end

#=
open(output_file_path, "w") do file
   write(file, "Computing time model 1= $(string(Model_1_time)), Post processing time= $(string(Post_Processing_time))")  # Convert to string and write
end
=#
using JSON
save_path = joinpath(folder_name_plot, "variables.json")
#Saving variables values in a json file
results = Dict("g" => g,
"rg" => rg,"rgp" => rgp, "re" => re, "rb" => rb, "pl" => pl, "pbc" => pbc, "pbd" => pbd,
               "eb" => eb, "Ppc"=>Ppc, "Ppd"=>Ppd,"Pener"=>Pener, "pe" => pe, "hfe" => hfe, "pe_c"=>pe_c, "hfgdinyec" => hfgdinyec,
               "hfgdcon" => hfgdcon, "hss" => hss, "RCU" => RCU,
               "zucvalues" => zucvalues, "wvalues" => wvalues, "vvalues" => vvalues,
               "ze" => ze, "zesu" => zesuvector,
               "zestb" => zestbvector, "zb" => zb, "costs_g" => sum(value.(g_costs))
               ,"zpcommit"=>zpcommit,"costs_rg" => sum(value.(rg_costs)), "costs_rb" => sum(value.(rb_costs)),
               "costs_re" => sum(value.(re_costs)), "costs_scu_UC" => sum(value.(scu_UC)),
               "costs_scu_UC_e" => sum(value.(scu_UC_e)), "costs_h" => sum(value.(h_costs)),
               "costs_total_reserve" => sum(value.(rb_costs)) + sum(value.(re_costs)) + sum(value.(rg_costs))+sum(value.(rgp_costs)),
               "hydrogenCost" => hydrogenCost,
               "Inertia_average"=>y_Inertia,
               "Inertia_average_energy"=>y_Inertia_energy,
               "Inertia_x_values"=>x_Inertia,
               "Total_Demand_Electro_storage"=>Total_Demand_Electro_storage,
               "Wind power" => pw,
               "Solar power" => ps,
         
               )

open(save_path, "w") do io
   JSON.print(io, results, 2)  # Pretty-print with indentation
end


status = termination_status(m)