
using Pkg
using CSV
using JSON

json_string = read("variables.json", String)

data = JSON.parse(json_string)

using DataFrames

Sum_Inertia_Vector = data["Sum_Inertia_Vector"]  # Dict{String, Any} with 24 entries
Sum_Inertia_Vector_energy = data["Sum_Inertia_Vector_energy"]  # Dict{String, Any} with 24 entries
rocof_post_opt = data["rocof_post_opt"]  # Dict{String, Any} with 24 entries
Inertia_average_energy= data["Inertia_average_energy"]  # Dict{String, Any} with 24 entries
Total_Demand_Electro_storage_MW = data["Total_Demand_Electro_storage_MW"]  # Net_Demand in GW
Wind_power = data["Wind power"]/1000  # Dict{String, Any} with 24 entries
Deltaf_nadir_re= data["Deltaf_nadir_re"]
Deltaf_nadir_re = Dict(k => v * -1 for (k, v) in Deltaf_nadir_re)
Deltaf_nadir_rg_2 = data["Deltaf_nadir_rg_2"]
Deltaf_nadir_rg_2 = Dict(k => v === nothing ? 0 : v for (k, v) in Deltaf_nadir_rg_2)
Demand=Total_Demand_Electro_storage_MW[1]/1000
Net_Demand=Demand-Wind_power

rb_per_hour = data["rb_per_hour"]  # Dict{String, Any} with 24 entries
rb_per_hour = [x[1] for x in rb_per_hour] 
re_per_hour = data["re_per_hour"]  # Dict{String, Any} with 24 entries
re_per_hour = [x[1] for x in re_per_hour]
rg_per_hour = data["rg_per_hour"]  # Dict{String, Any} with 24 entries
rg_per_hour = [x[1] for x in rg_per_hour]



pl= data["pl"]  


df = DataFrame(Hour = collect(keys(Sum_Inertia_Vector)), Inertia = collect(values(Sum_Inertia_Vector)), RoCoF= collect(values(rocof_post_opt)),Deltaf_nadir_re = collect(values(Deltaf_nadir_re)),Deltaf_nadir_rg_2=collect(values(Deltaf_nadir_rg_2)), Inertia_energy = collect(values(Sum_Inertia_Vector_energy)))

df.Hour = parse.(Int, df.Hour)

df_sorted = sort(df, :Hour)

df_sorted.Net_Demand = Net_Demand
df_sorted.pl = pl
df_sorted.rb_per_hour = rb_per_hour
df_sorted.re_per_hour = re_per_hour
df_sorted.rg_per_hour = rg_per_hour


CSV.write("Z.csv", df_sorted)