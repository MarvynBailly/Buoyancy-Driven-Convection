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
            help = "Number of days the simulation runs"
            arg_type = Float64

        "--outdir"
            help = "Path of directory to save outputs under"
            #default = "/glade/work/mbailly/Data/FrontalZone"   
            default = "/glade/derecho/scratch/mbailly/Buoyancy-Driven-Convection/Output"
            #default = "./Data/"
            arg_type = String

        "--flux"
            help = "Type of flux"
            default = "constant"
            arg_type = String

        "--SL"
            help = "Type of Sponge Layer to use"
            default = "pw"
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
sl_type  = args["SL"]
flux_type = args["flux"]

###########-------- SIMULATION PARAMETERS ----------------#############
include("simparams.jl")
simparams = SimParams()

group = casename[1:6]
group_symbol = Symbol(group)
possible_symbols = fieldnames(typeof(simparams))


if !(group_symbol in possible_symbols)
    error("The group symbol ", group_symbol, " is NOT in possible_symbols.")
end

@info "Loading $(group) with parameters:"
pm = getproperty(SimParams(), Symbol(group_symbol))

state_parameters = (; pm.N₀², pm.M², pm.f, pm.σ, pm.B₀, pm.flux_depth)
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

# assuming t is in seconds 
@inline function piecewise_flux(t,p)
    time_minutes_mod = mod(t/(60 * 60),24)
    if(time_minutes_mod <= 6)
        value = p.flux_depth #cool 
    elseif(time_minutes_mod > 6 && time_minutes_mod < 18)
        value = π*(p.B₀ - p.flux_depth)*sin((π / 12)*(time_minutes_mod - 6)) + p.flux_depth
    elseif(time_minutes_mod >= 18)
        value = p.flux_depth
    else
	    error("Buoyancy Flux Function Error")
        value = p.flux_depth
	end
    return value
end

if(flux_type == "constant")
    @info "     loading constant flux"
    @inline flux_forcing(t,p) = -p.B₀
elseif(flux_type == "pw1")
    @info "     loading $sl_type"
    @inline flux_forcing(t,p) = -piecewise_flux(t,p)
else
    error("Undefined flux type")
end

@inline surface_buoyancy_flux_amplitude(z,t,p) = flux_forcing(t,p)

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
@inline background_velocity(x, z, t, p) = (z) * (p.M² / p.f) 
V_field = BackgroundField(background_velocity, parameters = state_parameters)


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."

@inline pw_mask(x, z) = (-100 ≤ z ≤ -80) ? ((-80 - z) / 20)^2 : 0
@inline gs_mask(x, z) = exp(-(z - (-100))^2 / (2 * (6)^2))
@inline no_mask(x, z) = 0

@inline target_uvw(x, y, z) = 0
@inline target_b(x, y, z, p) = z*p.N₀²
if(sl_type == "pw")
    @info "     loading piecewise mask"
    @inline mask(x,z) = pw_mask(x,z)
elseif(sl_type == "gs")
    @info "     loading Gaussian mask"
    @inline mask(x,z) = gs_mask(x,z)
elseif(sl_type == "no")
    @info "     no sponge layer"
    @inline mask(x,z) = no_mask(x,z)
else
    error("Undefined sponger layer type")
end

@inline sponge_u(x, y, z, u, p)  = -p.σ * mask(x, z) * (u - target_uvw(x, y, z))
@inline sponge_v(x, y, z, v, p)  = -p.σ * mask(x, z) * (v - target_uvw(x, y, z))
@inline sponge_w(x, y, z, w, p)  = -p.σ * mask(x, z) * (w - target_uvw(x, y, z))
@inline sponge_b(x, y, z, b, p)  = -p.σ * mask(x, z) * (b - target_b(x, y, z, p))
 
Fu = Forcing(sponge_u, field_dependencies = (:u), parameters = state_parameters)
Fv = Forcing(sponge_v, field_dependencies = (:v), parameters = state_parameters)
Fw = Forcing(sponge_w, field_dependencies = (:w), parameters = state_parameters)
Fb = Forcing(sponge_b, field_dependencies = (:b), parameters = state_parameters)

sponge_forcing = (; u=Fu, v=Fv, w=Fw, b= Fb)
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

u,v,w = model.velocities
b = model.tracers.b

# Averages
mU = Field(Average(u,      dims=(1, 2)))
mV = Field(Average(v,      dims=(1, 2)))
mW = Field(Average(w,      dims=(1, 2)))
mN² = Field(Average(∂z(b),    dims=(1, 2)))
mb = Field(Average(b,         dims=(1, 2)))

# Perts
u′ = Field(u - mU)
v′ = Field(v - mV)
w′ = Field(w - mW)
b′ = Field(b - mb)


vw′dz = Field(Average(∂z(v′ * w′), dims=(1,2)))
fivw′dz = Field(1/pm.f * vw′dz)
uw′dz = Field(Average(∂z(u′ * w′), dims=(1,2)))
fiuw′dz = Field(1/pm.f * uw′dz)

# Mean Buoyancy Budget
κ = diffusivity(model.closure, model.diffusivity_fields, Val(:b))
bdz2 = Field(Average(∂z(∂z(b)),    dims=(1, 2)))
κbdz2 = Field(κ * bdz2)

wbdz = Field(Average(∂z(w′ * b′),    dims=(1, 2)))

uM2 = Field(mU * pm.M²)

# Mean Potential Vorticity
vt = Field(Average(v + model.background_fields.velocities.v,dims=(1,2)))
mvt = Field(Average(v + model.background_fields.velocities.v,dims=(1, 2)))
mvtdz = Field(Average(∂z(v + model.background_fields.velocities.v),dims=(1, 2)))
M2mvtdz = Field(pm.M² * mvtdz)

ζ = Field(∂x(v) - ∂y(u))
mζ = Field(Average(∂x(v) - ∂y(u), dims=(1,2)))
ζ′ = Field(ζ - mζ)
ζb′ = Field(Average(∂z(ζ′*b′), dims=(1,2)))

fbdz = Field(pm.f * mN²) 

# Mean Potential Vorticity Flux
fM2u = Field(pm.f*pm.M²*mU)

fwbdz = Field(pm.f * wbdz)

# TODO: Replace this with Tomas's method
b_flux1(model) = piecewise_flux(model.clock.time,state_parameters)




α = pm.αᵣ + pm.αₛ

fields_mean = Dict("u" => mU, "v" => mV, "w" => mW, "vt" => vt, "N2" => mN², "b" => mb, "fivw′dz" => fivw′dz, "fiuw′dz" => fiuw′dz, "κbdz2" => κbdz2, "wbdz" => wbdz, "uM2" => uM2, "M2mvtdz" => M2mvtdz, "ζb′" => ζb′, "fbdz" => fbdz, "fM2u" => fM2u, "fwbdz" => fwbdz, "b_flux" => b_flux1)

dims = Dict("b_flux" => ())

global_attributes = Dict("ν₀" => pm.ν₀, "κ₀" => pm.κ₀,
                         "B₀" => pm.B₀, "N₀²" => pm.N₀², "M²" => pm.M², "f" => pm.f,
                         "α" => α, "β" => pm.β)


simulation.output_writers[:averages] = NetCDFOutputWriter(model, fields_mean;
                                                       filename = casename * "_averages.nc",
                                                       dir = outdir,
                                                       dimensions = dims,
                                                       global_attributes = global_attributes,
                                                       schedule = TimeInterval(pm.out_interval_mean),
                                                       overwrite_existing = true)




# N²F = Field(∂z(b))
# uF = Field(u)
# vF = Field(v)
# wF = Field(w)
# bF = Field(b)

# fields = Dict("u" => uF, "v" => vF, "w" => wF, "b" => bF, "N2" => N²F)


# simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, dir = outdir, fields, filename=casename*".nc", schedule=TimeInterval(pm.out_interval_mean), overwrite_existing = true)

###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory used on a GPU run
@info "Run...."
run!(simulation)
@info "Simulation completed in |" * prettytime(simulation.run_wall_time)
