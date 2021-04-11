using Colors, Plots
using Logging
using ProgressBars

Logging.LogLevel(0)

include("soundfdtd.jl")
using .SoundFDTD


grid = SoundSpace2D(
    256,
    500,
    1.e-3,
    1.e-3,
    20.e-6,
    340,
    1.2
)

function normalize(x::Array{Real, 2})::Array{Real, 2}
    mx = max(x...)
    mi = min(x...)
    if mx == mi
        x
    else
        (x .- mi) / (mx - mi)
    end
end

function clip(
        x::Array{Real, 2},
        lower_bound::Real,
        upper_bound::Real
    )::Array{Real, 2}
    x[x .< lower_bound] .= lower_bound
    x[x .> upper_bound] .= upper_bound
    x
end


function play()
    anim = Animation()
    for i = ProgressBar(1:100)
        update_grid(grid)

        mat = clip(grid.P, -1.e30, 1e30)
        # mat = Gray.(255 .* normalize(grid.P))
        mat = Gray.(255 .* normalize(mat))
        frame(anim, plot(mat))
    end
    gif(anim, "fdtd.gif", fps=120)
end

play()