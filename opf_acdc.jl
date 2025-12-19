## Step 0: Activate environment
using Pkg
# Pkg.activate(@__DIR__)
# Pkg.instantiate()
# Pkg.update()
# Pkg.add("Ipopt")
# Pkg.add("PowerModels")
# Pkg.add("PowerModelsACDC")
# Pkg.add("JuMP")
using PowerModels, PowerModelsACDC, Gurobi, JuMP, CSV, DataFrames, CairoMakie
# Pkg.add("Plots") # if Plots package not added yet, for plotting results


gurobi = optimizer_with_attributes(Gurobi.Optimizer)

##### Step 1: Import the grid data and initialize the JuMP model
# Select the MATPOWER case file
path = pwd()
case_file = joinpath(path,"test_system.m")

# For convenience, use the parser of Powermodels to convert the MATPOWER format file to a Julia dictionary
data = PowerModels.parse_file(case_file)
ts = CSV.read("Load_data_hvdc.csv", DataFrame)
tsw= CSV.read("Wind_data_hvdc.csv", DataFrame)

# Initialize the JuMP model (an empty JuMP model) with defined solver
m = Model(gurobi)

##### Step 2: create the JuMP model & pass data to model
include(joinpath(path,"init_model.jl")) # Define functions define_sets! and process_parameters!
define_sets!(m, data,ts,tsw) # Pass the sets to the JuMP model
process_parameters!(m, data,ts,tsw) # Pass the parameters to the JuMP model

##### Step 3: Build the model
include(joinpath(path,"build_ac_opf_acdc_modify.jl")) # Define build_ac_opf_acdc! function
build_ac_opf_acdc_modify!(m) # Pass the model to the build_ac_opf_acdc! function

##### Step 4: Solve the model
optimize!(m) # Solve the model


G=m.ext[:sets][:G]
E=m.ext[:sets][:E]
T=m.ext[:sets][:t]
baseKG=m.ext[:parameters][:baseKG]
baseMVA=m.ext[:parameters][:baseMVA]
pg=JuMP.value.(m.ext[:variables][:pg])*baseMVA
pgvec= [pg[g,t] for g in G, t in T]
pev=JuMP.value.(m.ext[:variables][:pe])*baseMVA
pevec= [pev[e,t] for e in E, t in T]
hfe=JuMP.value.(m.ext[:variables][:hfe])*baseKG
hfevec= [hfe[e,t] for e in E, t in T]
hss=JuMP.value.(m.ext[:variables][:hss])*baseKG
hssvec= [hss[e,t] for e in E, t in T]
hfgconsum=JuMP.value.(m.ext[:variables][:hfgconsum])*baseKG
hfgconsumvec= [hfgconsum[e,t] for e in E, t in T]
hfginyect=JuMP.value.(m.ext[:variables][:hfginyect])*baseKG
hfginyectvec= [hfginyect[e,t] for e in E, t in T]



fig1 = Figure()
ax = fig1[1, 1] = Axis(fig1,
    title = "Storage and Hydrogen flows",
    xlabel = "Time (hours)",
    ylabel = "Hydrogen Flow (Kg/h) and Storage (Kg)"
)

lines!(ax, hfginyectvec[1, :], label = "Hydrogen injected")
lines!(ax, hssvec[1, :], label = "Hydrogen stored")
lines!(ax, hfevec[1, :], label = "Hydrogen Electrolyzer")
lines!(ax, hfgconsumvec[1, :], label = "Hydrogen consumed")

# Legend in separate panel
fig1[1, 2] = Legend(fig1, ax, "Hydrogen Flows", framevisible = false)

fig1

fig2 = Figure()
ax2 = fig2[1, 1] = Axis(fig2,
    title = "Storage and Hydrogen flows",
    xlabel = "Time (hours)",
    ylabel = "Hydrogen Flow (Kg/h) and Storage (Kg)"
)

lines!(ax2, hfginyectvec[2, :], label = "Hydrogen injected")
lines!(ax2, hssvec[2, :], label = "Hydrogen stored")
lines!(ax2, hfevec[2, :], label = "Hydrogen Electrolyzer")
lines!(ax2, hfgconsumvec[2, :], label = "Hydrogen consumed")

fig2[1, 2] = Legend(fig2, ax2, "Hydrogen Flows", framevisible = false)  
fig2


println(objective_value(m))





solution_summary(m)
#println(objective_value(m)) # Print the objective value of the model

