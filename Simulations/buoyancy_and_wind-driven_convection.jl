using ArgParse

using Oceananigans
using Oceananigans.Units

using Random
Random.seed!(11)
using Printf

using Statistics

###########-------- COMMAND LINE ARGUMENTS ----------------#############
@info "Parse command line arguments..."
# Returns a dictionary of command line arguments
function parse_command_line_arguments()
    settings = ArgParseSettings()
    @add_arg_table! settings begin
        "casename"
            help = "Name of simulation case"
            required = true
            arg_type = String

        "--nTf"
            help = "Number of inertial periods the simulation runs"
            arg_type = Float64

        "--outdir"
            help = "Path of directory to save outputs under"
           # default = "/glade/work/zhihuaz/Data/FrontalZone"   
            # default = "/glade/derecho/scratch/zhihuaz/FrontalZone/Output"
            default = "./Data/"
            arg_type = String
    end
    return parse_args(settings)
end

args = parse_command_line_arguments()
for (arg,val) in args
    @info "    $arg => $val"
end

casename = args["casename"]
outdir   = args["outdir"]

###########-------- SIMULATION PARAMETERS ----------------#############
@info "Load in simulation parameters"

include("simparams.jl")
group = ifelse(casename[3]=='1', "sim2D1","sim2D2")

@info "Loading $(group) with parameters:"

pm = getproperty(SimParams(), Symbol(group))

state_parameters = (; pm.N₀², pm.M², pm.f, pm.σ, pm.B₀)
for (param, val) in pairs(state_parameters)
    @info "     $param => $val"
end

###########-------- GRID SET UP ----------------#############
@info "Set up grid...."

# Normalized height ranging from 0 to 1
@inline h(k) = (k - 1) / pm.Nz

# Linear near-surface generator
@inline ζ₀(k) = 1 + (h(k) - 1) / pm.refinement

# Bottom-intensified stretching function
@inline Σ(k) = (1 - exp(-pm.stretching * h(k))) / (1 - exp(-pm.stretching))

# Generating function
@inline z_faces(k) = pm.Lz * (ζ₀(k) * Σ(k) - 1)


grid = RectilinearGrid(GPU(), size=(pm.Nx, pm.Nz), x=(0, pm.Nx), z=z_faces, topology=(Periodic, Flat, Bounded))

###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
state_parameters = (; pm.N₀², pm.M², pm.f, pm.σ, pm.B₀)

# TODO: surface buoyancy flux - add diurnal variability

@inline surface_buoyancy_flux_amplitude(z,t,p) = p.B₀

surface_buoyancy_flux = FluxBoundaryCondition(surface_buoyancy_flux_amplitude, parameters=state_parameters)

b_bcs = FieldBoundaryConditions(top = surface_buoyancy_flux, bottom = GradientBoundaryCondition(pm.N₀²))

free_slip_u = GradientBoundaryCondition(0.0)
free_slip_v = GradientBoundaryCondition(0.0)

u_bcs = FieldBoundaryConditions(top = free_slip_u, bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = free_slip_v, bottom = GradientBoundaryCondition(0))

###########-------- BACKGROUND FIELDS------------############

# background buoyancy 
@inline background_buoyancy(x,z,t,p) = p.M²*x;
B_field = BackgroundField(background_buoyancy, parameters = state_parameters)

# set the background
@inline Vt(x, z, t, p) = (z) * (p.M² / p.f) 
V_field = BackgroundField(Vt, parameters = state_parameters)


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."

# Quadratic presented in paper
# @inline (damping_func)(x, z, t, b, p) = (-100 ≤ z ≤ -80) ? -p.σ * ((-80 - z) / 20)^2 *(b - p.N₀²*z) : 0  


#### velocity sponge layer ####
@inline (mask)(x, z) =exp(-(z - (-100))^2 / (2 * (6)^2))

# TODO: couldn't pass the paramter in N₀²
@inline (target_b)(x,y,z) = z*9e-5  #z*pm.N₀²


uvw_sponge = Relaxation(rate=pm.σ, mask=mask, target = 0)
b_sponge = Relaxation(rate=pm.σ, mask=mask, target = target_b)#(, parameters=state_parameters))



# TODO: Add buoyancy sponge layer 
# causes slow time stepping down
#sponge_forcing = (; u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, b=b_sponge)
# double check that linear target/gaussian doesn't work 
# try adding sponge layer via sponging
# sponge_b(x,y,z,b,p) = rate(b - target) 
# than add as a forcing
# let dev know that 2D doesn't work with stuff
# TODO: how much the sponge layer effects the results.

# TODO: center around day 15 with a difference of one inenrtial 

# NOTE: You can see the overturning of the ML in the u and v plots


# seems to be okay without buoyancy
sponge_forcing = (; u=uvw_sponge, v=uvw_sponge, w=uvw_sponge)#, b=b_sponge)


###########-------- Model -----------------#############
@info "Setting up model...."

model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=pm.f),
                            buoyancy = BuoyancyTracer(),
                            tracers = (:b),
                            boundary_conditions = (u = u_bcs, v = v_bcs, b=b_bcs),
                            forcing = sponge_forcing,
                            advection = WENO(),
                            background_fields = (; v = V_field, b = B_field),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=pm.ν₀, κ=pm.κ₀))


###########-------- INITIAL CONDITIONS ---------------#############
@info "Defining initial conditions...."

# Add random flucatuations
u,v,w = model.velocities

vᵢ = rand(size(v)...) * 0.001
uᵢ = rand(size(u)...) * 0.001
wᵢ = rand(size(w)...) * 0.001


vᵢ .-= mean(vᵢ)
uᵢ .-= mean(uᵢ)
wᵢ .-= mean(wᵢ)

set!(model, v=vᵢ, u = uᵢ, w = wᵢ)

# set the initial buoyancy
b = model.tracers.b
@inline init_buoyancy(x,z) = pm.N₀²*z #+ pm.M²*x
set!(model, b = init_buoyancy)

###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."

stop_time = ifelse(args["nTf"] == nothing, pm.nTf, args["nTf"])days

Δx  = minimum_xspacing(grid, Center(), Center(), Center())
Δy  = minimum_yspacing(grid, Center(), Center(), Center())
Δz  = minimum_zspacing(grid, Center(), Center(), Center())
Δt₀ = pm.cfl * min(Δx, Δy, Δz) / max(pm.M² / pm.f, 0.02)seconds
# simulation = Simulation(model, Δt=Δt₀, stop_time=10minute, wall_time_limit=12hours)
simulation = Simulation(model, Δt=Δt₀, stop_time=stop_time, wall_time_limit=12hours)

wizard = TimeStepWizard(cfl=pm.cfl, diffusive_cfl=pm.cfl, min_change=0, max_change=5, max_Δt=pm.max_Δt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))



# Create a progress message.
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, progress % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        100 * (time(sim) / sim.stop_time),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))
                        
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))



###########-------- DIAGNOSTICS --------------#############
@info "Add diagnostics..."

mU = Field(Average(model.velocities.u, dims=(1, 2)))
mV = Field(Average(model.velocities.v, dims=1))
mW = Field(Average(model.velocities.w, dims=1))
mN² = Field(Average(∂z(model.tracers.b), dims=(1, 2)))
mb = Field(Average(model.tracers.b, dims=1))

N²F = Field(∂z(model.tracers.b))
uF = Field(model.velocities.u)
vF = Field(model.velocities.v)
wF = Field(model.velocities.w)
bF = Field(model.tracers.b)


fields = Dict("u" => uF, "v" => vF, "w" => wF, "b" => bF, "N2" => N²F)

fields_mean = Dict("u" => mU, "v" => mV, "w" => mW, "N2" => mN², "b" => mb)

simulation.output_writers[:averages] = NetCDFOutputWriter(model, fields_mean;
                                                       filename = casename * "_averages.nc",
                                                       dir = outdir,
                                                       schedule = TimeInterval(pm.out_interval_mean),
                                                       overwrite_existing = true)


simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, fields, filename=casename*".nc", schedule=TimeInterval(pm.out_interval_mean), overwrite_existing = true)

###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory used on a GPU run
@info "Run...."
run!(simulation)
@info "Simulation completed in |" * prettytime(simulation.run_wall_time)