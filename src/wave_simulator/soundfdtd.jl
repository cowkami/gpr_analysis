## Sound Simulator by FDTD method
module SoundFDTD
export SoundSpace2D, update_grid


struct SoundSpace2DConstants
    # space grid size
    nx::Int
    ny::Int
    
    # discretization size
    dx::Real  # [m]
    dy::Real  # [m]
    dt::Real  # [s]

    # parameters
    v::Real  # phase velocity [m/s]
    ρ::Real  # density [kg/m2]
    κ::Real  # elasticity [N/m2]
end

mutable struct SoundSpace2D
    # constants
    constants::SoundSpace2DConstants
    time::Real  # [s]

    # discritized space matrices
    Vx::Array{Real, 2}  # velocities x directions
    Vy::Array{Real, 2}  # velocities z directions
    P::Array{Real, 2}   # sound pressures

    function SoundSpace2D(
        nx::Int,
        ny::Int, 
        dx::Real,
        dy::Real,
        dt::Real,
        v::Real,
        ρ::Real,
        κ::Real=0.,
    )::SoundSpace2D
        # check stability
        C = v * dt * sqrt(1/dx^2 + 1/dy^2)
        if C > 1/sqrt(2)
            @warn "it may not be stable to calculate." 
        end

        if κ == 0.
            κ = v * ρ
        end
        constants = SoundSpace2DConstants(nx, ny, dx, dy, dt, v, ρ, κ)
        Vx = zeros(nx + 1, ny)
        Vy = zeros(nx, ny + 1)
        P = zeros(nx, ny)
        new(constants, 0., Vx, Vy, P)
    end
end


function _update_V(s::SoundSpace2D)::SoundSpace2D
    c = s.constants
    s.Vx[2:c.nx,  :    ] -= (c.dt / (c.ρ * c.dx)) * (s.P[2:c.nx,  :    ] - s.P[1:c.nx - 1,  :        ])
    s.Vy[ :    , 2:c.ny] -= (c.dt / (c.ρ * c.dy)) * (s.P[ :    , 2:c.ny] - s.P[ :        , 1:c.ny - 1])
    s
end

function _update_P(s::SoundSpace2D)::SoundSpace2D
    c = s.constants
    s.P -= (
        (c.κ * c.dt / c.dx) * (s.Vx[2:c.nx+1,  :    ] - s.Vx[1:c.nx,  :    ]) +
        (c.κ * c.dt / c.dy) * (s.Vy[ :    , 2:c.ny+1] - s.Vy[ :    , 1:c.ny])
    )

    # sound source
    A = 1.e-7
    f = 1.e+3
    if s.time < s.constants.dt * 10
        s.P[50, 50:100] .+= (
            ((c.κ * c.dt) / (c.dx * c.dy)) * A * sin(2*pi*s.time/25)
        )
    end
    s
end

function update_grid(s::SoundSpace2D)::SoundSpace2D
    _update_V(s)
    _update_P(s)
    s.time += s.constants.dt
    s
end


end  # module