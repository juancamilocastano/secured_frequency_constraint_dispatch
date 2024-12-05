# 24 GW of BESS by 2030 according to https://www.energy-storage.news/uk-battery-energy-storage-market-to-grow-to-24gw-by-2030-says-rystad-energy/
## Secure economic dispatch of great britain power system with 24 GW of BESS and 5 GW of electrolyzers
#Author: Juan Camilo Castano
#Date: 2024-02-12

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate()

## Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data.yaml"))

data_generators = CSV.read("Generators_data.csv", DataFrame)
data_electrolyzers = CSV.read("Electrolyzers_data.csv", DataFrame)
data_bess = CSV.read("BESS_data.csv", DataFrame)
ts = CSV.read("Load_data.csv", DataFrame)
data_system_parameters = CSV.read("System_parameters.csv", DataFrame)

## Step 2: create model & pass data to model
using JuMP
using Ipopt
m = Model(optimizer_with_attributes(Ipopt.Optimizer))


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

    #Battery energy storage systems
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

   #parameters of dispatchable generators per unit
   d = data["dispatchableGenerators"]

   #Parameter Fuel cost
   m.ext[:parameters][:FCOST]=Dict()
   m.ext[:parameters][:GmaxD]=Dict()
   m.ext[:parameters][:GminD]=Dict()
   m.ext[:parameters][:Dtg]=Dict()
   for i in ID
      if length(i)==6 
         m.ext[:parameters][:FCOST][i] =  d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
         m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"] 
         m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
         m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
      else 
         if length(i)==7 
            m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-3)]["marginalcost"] #	β Marginal fuel cost
            m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-3)]["maxPowerOutput"]
            m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-3)]["minStableOperatingPoint"]
            m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-3)]["deploymentTime"]
         else
            if length(i)==9
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-2)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-2)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
            else
               if length(i)==8
               m.ext[:parameters][:FCOST][i] = d[SubString(i,1:length(i)-4)]["marginalcost"] #	β Marginal fuel cost
               m.ext[:parameters][:GmaxD][i] =  d[SubString(i,1:length(i)-4)]["maxPowerOutput"]
               m.ext[:parameters][:GminD][i] =  d[SubString(i,1:length(i)-4)]["minStableOperatingPoint"]
               m.ext[:parameters][:Dtg][i] =  d[SubString(i,1:length(i)-4)]["deploymentTime"]           
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

   #Parameter variable generators

   #Parameter electrolyzers
   d = data["Electrolyzer"]
   m.ext[:parameters][:PEmax] = Dict(i => d[SubString(i,1:length(i)-2)]["installed_capacity"] for i in ID_E)
   m.ext[:parameters][:PEmin] = Dict(i => d[SubString(i,1:length(i)-2)]["minimun_consumption"] for i in ID_E)
   m.ext[:parameters][:Eeff] = Dict(i => d[SubString(i,1:length(i)-2)]["efficiency"] for i in ID_E)
   m.ext[:parameters][:Eload_factor] = Dict(i => d[SubString(i,1:length(i)-2)]["load_factor_electrolyzer"] for i in ID_E) # load factor of electrolyzers
   m.ext[:parameters][:Max_h_f] = Dict(i => d[SubString(i,1:length(i)-2)]["max_hydrogen_flow"] for i in ID_E) #maximum hydrogen flow in kg/h
   m.ext[:parameters][:Max_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["max_hydrogen_storage"] for i in ID_E) #maximum hydrogen storage in kg
   m.ext[:parameters][:Min_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["min_hydrogen_storage"] for i in ID_E) #minimum hydrogen storage in kg
   m.ext[:parameters][:Ini_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["initial_hydrogen_storage"] for i in ID_E) #initial hydrogen storage in kg
   m.ext[:parameters][:End_h_s] = Dict(i => d[SubString(i,1:length(i)-2)]["final_hydrogen_storage"] for i in ID_E) #final hydrogen storage in kg
   m.ext[:parameters][:Dte] = Dict(i => d[SubString(i,1:length(i)-2)]["deploymentTime"] for i in ID_E) #final hydrogen storage in kg



   d = data

   #Parameter Fuel cost

   m.ext[:parameters][:PBmax]=Dict()
   m.ext[:parameters][:EBmax]=Dict()
   m.ext[:parameters][:Beff]=Dict()
   m.ext[:parameters][:Ini_e_b]=Dict()
   m.ext[:parameters][:End_e_b]=Dict()
   m.ext[:parameters][:Dtb]=Dict()

   for i in ID_BESS
      if length(i)==6 
         m.ext[:parameters][:PBmax][i] =  d[SubString(i,1:length(i)-2)]["Pmax"]
         m.ext[:parameters][:EBmax][i] =  d[SubString(i,1:length(i)-2)]["Emax"]
         m.ext[:parameters][:Beff][i] =  d[SubString(i,1:length(i)-2)]["eff"]
         m.ext[:parameters][:Ini_e_b][i] =  d[SubString(i,1:length(i)-2)]["Einit"]
         m.ext[:parameters][:End_e_b][i] =  d[SubString(i,1:length(i)-2)]["Efinal"]
         m.ext[:parameters][:Dtb][i] =  d[SubString(i,1:length(i)-2)]["deploymentTime"]
      else 
         if length(i)==7 
         m.ext[:parameters][:PBmax][i] =  d[SubString(i,1:length(i)-3)]["Pmax"]
         m.ext[:parameters][:EBmax][i] =  d[SubString(i,1:length(i)-3)]["Emax"]
         m.ext[:parameters][:Beff][i] =  d[SubString(i,1:length(i)-3)]["eff"]
         m.ext[:parameters][:Ini_e_b][i] =  d[SubString(i,1:length(i)-3)]["Einit"]
         m.ext[:parameters][:End_e_b][i] =  d[SubString(i,1:length(i)-3)]["Efinal"]
         m.ext[:parameters][:Dtb][i] =  d[SubString(i,1:length(i)-3)]["deploymentTime"]
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

   #OJO: no se convirtio el diccionario de eficiencia y de deployment time en Dict{String, Int64}
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

# Extract time series data
D= m.ext[:timeseries][:D]
