using Parameters

@with_kw struct SimParams
    commons =(;
    ##### Table Parameters #####
    Nx = 1024,   # numbers of points in hor dir
    Nz = 128,         # number of points in vert

    Lx = 1000meters, #1000meters        # domain horizontal extents
    Lz = 100meters,        # depth

    f = 1e-4,
    N₀²=9e-5,

    ν₀ = 1.0e-4, # [m² s⁻¹] viscosity
    κ₀ = 1.0e-4, # [m² s⁻¹] diffusivity - double check


    ##### Sponge Layer ##### 
    refinement = 11,  # controls spaczeting near surface (higher means finer spaced)
    stretching = 1.6, # controls rate of strecthing at bottom
    σ = 0.005,


    ##### Simulation #####
    nTf = 20days,

    cfl = 0.75,
    max_Δt = 5minutes,

    out_interval_mean = 5minutes

    )

    sim2D1 = (; commons...,
        #2D₁
        M² = 0,
        B₀ = -4.24e-8,
        αᵣ= 0.08,
        αₛ = 0.22,
        β = 0,
    )

    sim2D2 = (; commons...,
        #2D₂
        M² = -4.24e-7,
        B₀ = -4.24e-8,
        αᵣ= 0.00,
        αₛ = 0.21,
        β = -0.09,
    )

    sim2D3 = (; commons...,
        #2D₃
        M² = -2.12e-7,
        B₀ = -4.24e-8,
        αᵣ= 0.02,
        αₛ = 0.21,
        β = -0.04,
    )

    sim2D4 = (; commons...,
        #2D₄
        M² = -8.24e-7,
        B₀ = -4.24e-8,
        αᵣ= 0.00,
        αₛ = 0.21,
        β = -0.20,
)

    sim2D5 = (; commons...,
        #2D₅
        M² = -4.24e-7,
        B₀ = -8.48e-8,
        αᵣ= 0.01,
        αₛ = 0.11,
        β = -0.03,
    )

    sim2D6 = (; commons...,
        #2D₆
        M² = -4.24e-7,
        B₀ = -2.12e-8,
        αᵣ= 0.00,
        αₛ = 0.42,
        β = -0.12,
    )
end
