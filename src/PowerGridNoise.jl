module PowerGridNoise
    using DifferentialEquations
    using Distributions

    ## TODO: add rng to make the trajectories reproducible :-)

    include("frequency_dynamics.jl")
    export pg_frequency_model

    include("load_profile.jl")
    export load_profile_model, ornstein_uhlenbeck, ornstein_uhlenbeck_ϵ
end
