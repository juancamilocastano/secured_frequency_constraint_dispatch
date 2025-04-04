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

    ID_Nuclear = Array{Union{Nothing,String}}(nothing,0)
    for i in 1:data["dispatchableGenerators"]["Nuclear"]["numberOfUnits"]
      ID_Nuclear = m.ext[:sets][:ID_Nuclear] = push!(ID_Nuclear,string("Nuclear_$(i)"))
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
    m.ext[:timeseries][:SC] =  capacity_factor_renewable.Solar_Capacity_Factor[1:24]
    m.ext[:timeseries][:WC] =  capacity_factor_renewable.Wind_Capacity_Factor[1:24]
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

# Extract time series data and convert them in PU values
D= m.ext[:timeseries][:D]
Pbase=maximum(D) 
D= m.ext[:timeseries][:D]/Pbase
#Extrac capacity factor renewable
SC= m.ext[:timeseries][:SC]
WC= m.ext[:timeseries][:WC]
Installed_S = m.ext[:parameters][:Installed_S]/Pbase
Installed_W = m.ext[:parameters][:Installed_W]/Pbase


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

startupCost=m.ext[:parameters][:startupCost]


MUT=m.ext[:parameters][:MUT]


MDT=m.ext[:parameters][:MDT]


maxFrdelivarable=m.ext[:parameters][:maxFrdelivarable]


notloadCost=m.ext[:parameters][:notloadCost]



#Dtg=m.ext[:parameters][:Dtg]
Dtg=15

inertia_Constant=m.ext[:parameters][:inertia_Constant]
Inertia_Vector= Dict(k => inertia_Constant[k] * GmaxD[k] for k in keys(inertia_Constant) ∩ keys(GmaxD))

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
heffc =Dict(key => value   for (key, value) in heffc)

heffd = m.ext[:parameters][:heffd]
heffd =Dict(key => value  for (key, value) in heffd)

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



#Create folder to save the results

# create variables 

zuc = m.ext[:variables][:zuc] = @variable(m, [i=ID,j=J], binary=true, base_name="commitment")
v = m.ext[:variables][:v] = @variable(m, [i=ID,j=J], binary=true, base_name="start_up")
w = m.ext[:variables][:w] = @variable(m, [i=ID,j=J], binary=true, base_name="shoot_down")
g = m.ext[:variables][:g] = @variable(m, [i=ID,j=J],lower_bound=GminD[i], base_name="generation") #Power generation generators
x = m.ext[:variables][:x] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="x") #Auxiliary variable rotate second order cone
y = m.ext[:variables][:y] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="y") #Auxiliary variable rotate second order cone
z = m.ext[:variables][:z] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="z") #Auxiliary variable rotate second order cone
rg = m.ext[:variables][:rg] = @variable(m, [i=ID,j=J],lower_bound=0, base_name="rg") #Reserve provided by generators
rb = m.ext[:variables][:rb] = @variable(m, [i=ID_BESS,j=J],lower_bound=0, base_name="rb") #Reserve provided by batteries
re = m.ext[:variables][:re] = @variable(m, [i=ID_E,j=J],lower_bound=0, base_name="re") #REserve provided by electrolyzers
pl = m.ext[:variables][:pl] = @variable(m, [j=J],base_name="pl") #loss of generation
pbc = m.ext[:variables][:pbc] = @variable(m, [i=ID_BESS,j=J],lower_bound=0,upper_bound=PBmax[i], base_name="pbc") #Charging power of the batteries
pbd = m.ext[:variables][:pbd] = @variable(m, [i=ID_BESS,j=J],lower_bound=0,upper_bound=PBmax[i], base_name="pbd") #Discharging power of the batteries
eb = m.ext[:variables][:eb] = @variable(m, [i=ID_BESS,j=J], lower_bound= EBmax[i]*(1-DOD_max[i]), upper_bound=EBmax[i] , base_name="eb") #Energy bounds of the batteries
hfe= m.ext[:variables][:hfe] = @variable(m, [i=ID_E,j=J],lower_bound=0, upper_bound= Max_h_f[i], base_name="hfe") #Hydrogen flow limit of the hydrogen produced by electrolyzers
hfg = m.ext[:variables][:hfg] = @variable(m, [i=ID_E,j=J],lower_bound=-Max_h_f[i],upper_bound= Max_h_f[i],base_name="hfg") #Hydrogen flow limit of the hydrogen flowing trhow the hydrogen pipeline
#hfg = m.ext[:variables][:hfg] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound= 0,base_name="hfg") #Hydrogen flow limit of the hydrogen flowing trhow the hydrogen pipeline
hfsc = m.ext[:variables][:hfsc] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound= Max_h_f[i],base_name="hfsc") #Hydrogen flow storage charging
hfsd = m.ext[:variables][:hfsd] = @variable(m, [i=ID_E,j=J],lower_bound=0,upper_bound= Max_h_f[i],base_name="hfsd") #Hydrogen flow storage discharging
pe = m.ext[:variables][:pe] = @variable(m,  [i=ID_E,j=J],lower_bound= PEmin[i], upper_bound=PEmax[i], base_name="pe") #Power consumption electrolyzer
pe_c= m.ext[:variables][:pe_c] = @variable(m,  [i=ID_E,j=J], base_name="pe_c") #Power consumption of compressor electrolyzer 
hss = m.ext[:variables][:hss] = @variable(m, [i=ID_E,j=J],lower_bound=Min_h_s[i], upper_bound= Max_h_s[i], base_name="hss") #hydrogen storage limit




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
scu_UC=m.ext[:expressions][:scu_UC] = @expression(m, [i=ID,j=J], startupCost[i]*v[i,j])


#Create folder to save the results
RE_costs=res_cost_e["E_500_1"]
RG_costs=res_cost_g["CCGT_77"]
Installed_W_F=Installed_W*Pbase
Installed_S_F=Installed_S*Pbase
folder_name_plot="Results_UC_CRE_$(RE_costs)_CRG_$(RG_costs)_IW_$(Installed_W_F)_IS_$(Installed_S_F)"
mkdir(folder_name_plot)


#It creates the expression of the inertia
Inertia_Expression=m.ext[:expressions][:Inertia_Expression] = @expression(m, [i=ID,j=J],Inertia_Vector[i]*zuc[i,j])

#Create objective function

obj= m.ext[:objective] = @objective(m,Min, sum(g_costs)+sum(rg_costs)+sum(rb_costs)+sum(re_costs)+sum(scu_UC)-sum(h_costs) #Objective function
)


con1=m.ext[:constraints][:con1] = @constraint(m, [j=J],
WC[j]*Installed_W+sum(g[i,j] for i in ID) -sum(pbc[i,j] for i in ID_BESS)+sum(pbd[i,j] for i in ID_BESS) == D[j]-SC[j]*Installed_S+sum(pe[i,j] for i in ID_E)+sum(pe_c[i,j] for i in ID_E)
)

#Constraint upper bound generators power
con2_1=m.ext[:constraints][:con2_1] = @constraint(m, [i=ID,j=J],
g[i,j]+rg[i,j].<=GmaxD[i]*zuc[i,j]
)
#Constraint upper bound generators power
con2_2=m.ext[:constraints][:con2_2] = @constraint(m, [i=ID,j=J],
GminD[i]*zuc[i,j].<=g[i,j]+rg[i,j]
)

#Initial status generators (unit commitment)
con2_3=m.ext[:constraints][:con2_3] = @constraint(m, [i=ID,j=J[1]],1-zuc[i,j]+v[i,j]-w[i,j]==0)

#Status generators (unit commitment)
con2_4=m.ext[:constraints][:con2_4] = @constraint(m, [i=ID,j=J[2:end]],
zuc[i,j-1]-zuc[i,j]+v[i,j]-w[i,j]==0
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
#=
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
con5=m.ext[:constraints][:con5] = @constraint(m, [i=ID_E,j=J],re[i,j].<=pe[i,j]-PEmin[i]
)

#Constraint end energy value of the batteries
con6=m.ext[:constraints][:con6] = @constraint(m, [i=ID_BESS,j=J[end]],End_e_b[i]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])
#con6=m.ext[:constraints][:con6] = @constraint(m, [i=ID_BESS,j=J[end]],End_e_b[i]==eb[i,j])


#Constraint initial value energy of the batteries
#con7=m.ext[:constraints][:con7] = @constraint(m, [i=ID_BESS,j=J[1]],eb[i,j+1]-Ini_e_b[i]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])
con7=m.ext[:constraints][:con7] = @constraint(m, [i=ID_BESS,j=J[1]],eb[i,j]==Ini_e_b[i])

#Constraint charging-discharging batteries
con8=m.ext[:constraints][:con8] = @constraint(m, [i=ID_BESS,j=J[1:end-1]],eb[i,j+1]-eb[i,j]==Beffc[i]*pbc[i,j]-pbd[i,j]/Beffd[i])

##Constraint Electrolyzer
#Constraint  hydrogen production electrolyzer
con9=m.ext[:constraints][:con9] = @constraint(m, [i=ID_E,j=J],hfe[i,j]==pe[i,j]/Eeff[i]
)

#Mass hydrogen equation
#The term Eload_factorl[i]*PEmax[i]/Eeff[i] correspond to the hydrogen demand
con10=m.ext[:constraints][:con10] = @constraint(m, [i=ID_E,j=J],hfe[i,j]+hfsd[i,j]==hfsc[i,j]+hfg[i,j]+Eload_factor[i]*PEmax[i]/(Eeff[i]))

#Hydrogen storage constraints


#Constraint end hydrogen value of the hydrogen storage
con11=m.ext[:constraints][:con11] = @constraint(m, [i=ID_E,j=J[end]],End_h_s[i]==hss[i,j]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])
#con11=m.ext[:constraints][:con11] = @constraint(m, [i=ID_E,j=J[end]],End_h_s[i]==hss[i,j])


#Constraint initial value of the hydrogen storage
con12=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[1]],hss[i,j+1]==Ini_h_s[i]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])
#con12=m.ext[:constraints][:con12] = @constraint(m, [i=ID_E,j=J[1]],hss[i,j]==Ini_h_s[i])

#Constraint charging-discharging of the hydrogen storage
con13=m.ext[:constraints][:con13] = @constraint(m, [i=ID_E,j=J[1:end-1]],hss[i,j+1]==hss[i,j]+hfsc[i,j]*heffc[i]-hfsd[i,j]/heffd[i])

#Constraint_commitment electrolyzer

con13_1=m.ext[:constraints][:con13_1] = @constraint(m, [i=ID_E,j=J],pe[i,j].<=PEmax[i])
con13_2=m.ext[:constraints][:con13_2] = @constraint(m, [i=ID_E,j=J],pe[i,j].>=PEmin[i])

#Power consumption in standby=5% of nominal power according to chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2306.10962




#Status constraints electrolyzer


#Power consumption of the electrolyzer
con13_8=m.ext[:constraints][:con13_8] = @constraint(m, [i=ID_E,j=J],pe_c[i,j]==compresor_power[i]*hfe[i,j])


#Constraint maximum frequency variation
con14= m.ext[:constraints][:con14] = @constraint(m, [i=ID,j=J],((sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO*(sum(rb[:, j])/Dtb+sum(re[:, j])/Dte+(sum(rg[:, j])-rg[i, j])/Dtg)).>=pl[j]^2/(4*deltaf))

#constraints nadir occurrence time

con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=0)

con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<= 0.00000000000001+sum(re[:, j])+sum(rb[:,j])*(Dte/Dtb)+(sum(rg[:, j])-rg[i, j])*Dte/Dtg)

#constraint ROCOF

con17=m.ext[:constraints][:con17] = @constraint(m, [i=ID,j=J], pl[j]*FO/2 .<=rocofmax*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j]))

#QSS frequency constraint
con18=m.ext[:constraints][:con18] = @constraint(m, [i=ID,j=J],pl[j].<=0.00000000000001+sum(re[:, j])+sum(rb[:, j])+(sum(rg[:, j])-rg[i, j]))




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
 end

#Constraints frequency nadir as rotate second order cone
 con14_1=m.ext[:constraints][:con14_1] = @constraint(m,[i=ID,j=J],y[i,j] ==2*deltaf*(sum(Inertia_Expression[:,j])-Inertia_Expression[i,j])/FO-sum(re[:, j])*Dte/2-sum(rb[:, j])*Dtb/2)
 con14_2=m.ext[:constraints][:con14_2] = @constraint(m,[i=ID,j=J],z[i,j] ==(sum(rg[:, j])-rg[i, j])/Dtg)
 con14_3=m.ext[:constraints][:con14_3] = @constraint(m,[i=ID,j=J],[y[i,j]; z[i,j]; pl[j]-sum(re[:, j])-sum(rb[:,j])] in RotatedSecondOrderCone())

 #constraints nadir occurrence time
 con15= m.ext[:constraints][:con15] = @constraint(m, [i=ID,j=J],pl[j].>=sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j])*Dtb/Dtg)
 con16=m.ext[:constraints][:con16] = @constraint(m, [i=ID,j=J],pl[j].<=0.000001 + sum(rb[:,j])+sum(re[:,j])+(sum(rg[:, j])-rg[i, j]))



 output_file_path = joinpath(folder_name_plot, "output_3.txt")
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
hfg= value.(m.ext[:variables][:hfg])*Mbase
hfsc= value.(m.ext[:variables][:hfsc])*Mbase
hfsd= value.(m.ext[:variables][:hfsd])*Mbase
hss= value.(m.ext[:variables][:hss])*Mbase
zucvalues=  value.(m.ext[:variables][:zuc])
pe_c= value.(m.ext[:variables][:pe_c])*Pbase




Ivec = [i for  i in I]
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
hfgvec = [hfg[i,j] for  i in ID_E, j in J]
hfscvec = [hfsc[i,j] for  i in ID_E, j in J]
hfsdvec = [hfsd[i,j] for  i in ID_E, j in J]
hssvec = [hss[i,j] for  i in ID_E, j in J]
plvec= [pl[j] for j in J]
ps=SC*Installed_S*Pbase
pw=WC*Installed_W*Pbase



using StatsPlots


p_generators = groupedbar(transpose(gvec[:,:]),
    bar_position = :stack,
    label = permutedims(Ivec),
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

p_re = groupedbar(transpose(revec[:,:]),
    bar_position = :stack,
    label = permutedims(Evec),
    legend = :outertopright,
    xlabel = "Time [h]",
    ylabel = "Procured reserve [MW]",
    title = "Hourly procured reserve Electrolyzers");
plot(p_re, layout = (1,2), size=(500, 500))

# Plot the vectors
p_HF = plot(hfevec[1,:], label = "hfevec", lw = 2)
plot!(hfgvec[1,:], label = "hfgvec", lw = 2)
plot!(hfscvec[1,:], label = "hfscvec", lw = 2)
plot!(hfsdvec[1,:], label = "hfsdvec", lw = 2)
plot!(hssvec[1,:], label = "hssvec", lw = 2)

plot!(xlabel = "Time [h]",
      ylabel = "Hydrogen mass [kg]",
      legend = :outertopright,
      title = "Hourly Hydrogen mass flow")



save_path = joinpath(folder_name_plot, "Produced_power_plot.png")
savefig(p_generators, save_path)

save_path = joinpath(folder_name_plot, "BESS_reserve_plot.png")
savefig(p_RB, save_path)

save_path = joinpath(folder_name_plot, "Electrolyzer_Freserve_plot.png")
savefig(p_re, save_path)

save_path = joinpath(folder_name_plot, "Hydrogen_flows_plot.png")
savefig(p_HF, save_path)


sum_reserve_g = vec(sum(rgvec, dims=1))
sum_reserve_b = vec(sum(rbvec, dims=1))
sum_reserve_e = vec(sum(revec, dims=1))
total_reserve = sum_reserve_g + sum_reserve_b + sum_reserve_e

#It creates a folder to store the results




pr=bar(1:24, [sum_reserve_g sum_reserve_b sum_reserve_e],
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


plot!(pr, 1:24, plvec,
    label = "Loss of power",
    linewidth = 2,
    color = :purple,
    linestyle = :solid)

plot!(pr, 1:24, total_reserve,
label = "Total reserve",
linewidth = 2,
color = :red,
linestyle = :dash)

display(pr)
save_path = joinpath(folder_name_plot, "reserves_plot.png")
savefig(pr, save_path)




p_bat=plot(pbcvec[1,:], 
     seriestype = :scatter,  # Use scatter instead of a line
     marker = :cross,        # Cross-shaped marker
     markersize = 3,         # Slightly smaller for a clean look
     markercolor = :black,   # Black for academic clarity
     legend = :outerright,
     title = "Power and Energy of the BESS", 
     xlabel = "Time [H]", 
     ylabel = "Power [MW], Energy [MWh]", 
     label = "Charging power")

plot!(ebvec[1,:],
linewidth = 2,
color = :red,
label = "Battery energyr"
)

# Add pbdvec[1,:] to the same plot with a different cross style
plot!(pbdvec[1,:], 
      seriestype = :scatter, 
      marker = :xcross,       # Diagonal cross for distinction
      markersize = 3,         # Consistent size
      markercolor = :blue,    # Blue to differentiate
      label = "Discharging power")  # Different label

# Save the combined plot

display(p_bat)

save_path = joinpath(folder_name_plot, "BESS_power_energy_plot.png")
savefig(p_bat, save_path)

Sum_Inertia_Vector=Dict()

for j in J
   Sum_Inertia_Vector[j]=sum(Inertia_Vector[i]*value.(zuc[i,j]) for i in ID)  
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

end


#Making file with the computing time
output_file_path = joinpath(folder_name_plot, "Computing_time.txt")
# Ensure the folder exists (optional, but recommended)
mkpath(folder_name_plot)
# Use output_file_path in open
open(output_file_path, "w") do file
    write(file, "Computing time model 1= $(string(Model_1_time)), Computing time model 3= $(string(Model_3_time)), Post processing time= $(string(Post_Processing_time))")  # Convert to string and write
end



