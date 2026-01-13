using WGLMakie, Bonito, JLD2, Rotations, Observables
WGLMakie.activate!(resize_to=:parent)

"""
    rotate_points(points, roll, pitch, yaw)

Rotate a vector of Point3f by Euler angles (in radians).
Uses ZYX convention: yaw (Z) → pitch (Y) → roll (X).
"""
function rotate_points(points::Vector{Point3f}, roll::Real, pitch::Real, yaw::Real)
    R = RotZYX(yaw, pitch, roll)
    return [Point3f(R * Vec3f(p...)) for p in points]
end

"""
    rotate_points_deg(points, roll, pitch, yaw)

Rotate a vector of Point3f by Euler angles (in degrees).
Uses ZYX convention: yaw (Z) → pitch (Y) → roll (X).
"""
function rotate_points_deg(points::Vector{Point3f}, roll::Real, pitch::Real, yaw::Real)
    return rotate_points(points, deg2rad(roll), deg2rad(pitch), deg2rad(yaw))
end

"""
    rotate_configurations(configurations, roll, pitch, yaw)

Rotate all configurations by Euler angles (in degrees).
"""
function rotate_configurations(configurations::Vector{Vector{Point3f}}, roll::Real, pitch::Real, yaw::Real)
    return [rotate_points_deg(config, roll, pitch, yaw) for config in configurations]
end

function get_point3f_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [Point3f(e) for e in eachcol(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...))]
end

open("docking_sim.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session
        conf_grad = cgrad(:Dark2_3, 3, categorical=true)
        confc1 = conf_grad[1]
        confc2 = conf_grad[2]
        confc3 = conf_grad[3]
        confmap = [confc1, confc2, confc3]

        @load "jld2s/1a30_docking.jld2" input output

        template_centers = input["template_centers"]
        template_radii = input["template_radii"]
        vis_range = 200:argmin(output["Es"])
        fig = Figure(size=(800, 800))

        configurations = [get_point3f_realization(state, template_centers) for state in output["states"][vis_range]]
        configurations = rotate_configurations(configurations, 30, -120, 0)

        index_slider = Slider(1:length(configurations))
        points = map(index_slider) do idx
            configurations[idx]
        end

        conf_scene = LScene(fig[1, 1], show_axis = false)
        colors = vcat([fill(confmap[i], length(template_radii[i])) for i in eachindex(template_radii)]...)
        meshscatter!(conf_scene, points, markersize = vcat(template_radii...), color = colors, colormap = :rainbow)

        slider_style = DOM.style("""
            .slider-container {
                display: flex;
                align-items: center;
                padding: 15px 20px;
            }
            input[type="range"] {
                -webkit-appearance: none;
                appearance: none;
                width: 100%;
                height: 10px;
                background: #e0e0e0;
                border-radius: 5px;
                outline: none;
                cursor: pointer;
            }
            input[type="range"]::-webkit-slider-runnable-track {
                height: 10px;
                border-radius: 5px;
                background: #e0e0e0;
            }
            input[type="range"]::-webkit-slider-thumb {
                -webkit-appearance: none;
                appearance: none;
                width: 14px;
                height: 14px;
                background: #4063d8;
                border: none;
                border-radius: 50%;
                cursor: pointer;
                margin-top: -2px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.3);
            }
            input[type="range"]::-webkit-slider-thumb:hover {
                background: #3050c0;
            }
            input[type="range"]::-moz-range-thumb {
                width: 14px;
                height: 14px;
                background: #4063d8;
                border: none;
                border-radius: 50%;
                cursor: pointer;
                box-shadow: 0 1px 3px rgba(0,0,0,0.3);
            }
            input[type="range"]::-moz-range-thumb:hover {
                background: #3050c0;
            }
            input[type="range"]::-moz-range-progress {
                background: lightblue;
                height: 10px;
                border-radius: 5px;
            }
            input[type="range"]::-moz-range-track {
                background: #e0e0e0;
                height: 10px;
                border-radius: 5px;
            }
        """)
        slider = DOM.div(slider_style, DOM.div(index_slider, class="slider-container"))
        Bonito.record_states(session, DOM.div(slider, fig))
    end
    show(io, MIME"text/html"(), app)
    println(io, """
        </body>
    </html>
    """)
end
