using WGLMakie, Bonito, JLD2, Rotations, Observables
WGLMakie.activate!(resize_to=:parent)

function get_point3f_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x)
        return [Point3f(e) for e in eachcol(hvcat((n_mol), [R * template_centers .+ t for (R,t) in x]...))]
    else
        @assert false "This function is not implemented for single hard spheres."
    end
end

open("tmv_sim.html", "w") do io
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
        confmap = [confc1, confc3, confc2]

        @load "jld2s/6r7m_with_p.jld2" input output
        fig = Figure(size=(800, 800))

        template_centers = input["template_centers"]
        template_radii = input["template_radii"]
        n_atoms_per_mol = length(template_radii)
        n_mol = input["n_mol"]
        configurations = [get_point3f_realization(state, template_centers) for state in output["states"]]

        index_slider = Slider(1:length(output["states"]))
        points = map(index_slider) do idx
            configurations[idx]
        end

        conf_ax = LScene(fig[1, 1], show_axis = false)
        conf_colors = vcat([[j for _ in 1:n_atoms_per_mol] for j in 1:n_mol]...)
        meshscatter!(conf_ax, points, markersize = vcat([template_radii for _ in 1:n_mol]...), color = conf_colors, colormap = confmap)

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
