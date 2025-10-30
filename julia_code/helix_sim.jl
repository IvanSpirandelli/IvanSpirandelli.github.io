using WGLMakie, Bonito, FileIO, Observables, JLD2, MorphoMol, Distances
WGLMakie.activate!()

open("interactive_plot_buttons_fixed.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session

        ma_grad = cgrad(:Set1_7, 7, categorical=true)
        fsolc = ma_grad[5]

        @load "jld2s/helix_assembly.jld2" output n template_centers template_radii rs pf bounds
        vis_range =  1:25:argmin(output["Es"])
        fig = Figure(fontsize = 16)

        # ---- Data Preparation ----
        xs = 1:length(vis_range)
        Sums = output["Es"][vis_range]
        Sums = Sums .- Sums[1]
        n_mol = 8
        configurations = [[Point3f(e) for e in eachcol(reshape(x, (3, length(x)÷3)))] for x in output["states"][vis_range]]

        # ---- State Management and Controls ----

        # The reactive index that drives the entire visualization
        obs_index = Observable(1)

        # UI Buttons
        prev_button = Button("Previous")
        next_button = Button("Next")

        # ---- Corrected Event Listeners ----
        # The `on` function listens for changes in `button.value`.
        # The block receives the new value (`is_clicked`, which will be `true`).
        # We explicitly check `if is_clicked` to ensure we only act on the press.

        on(prev_button.value) do is_clicked
            if is_clicked
                # Decrement index, but not below 1
                if obs_index[] > 1
                    obs_index[] = obs_index[] - 1
                end
            end
        end

        on(next_button.value) do is_clicked
            if is_clicked
                # Increment index, but not beyond the number of frames
                if obs_index[] < length(vis_range)
                    obs_index[] = obs_index[] + 1
                end
            end
        end

        # ---- Declarative Observables using @lift ----
        # These are linked to `obs_index` and will update synchronously
        # whenever a button click successfully changes the index.

        points = @lift(configurations[$obs_index])

        connections_for_all = @lift begin
            c = configurations[$obs_index]
            segments = Point3f[]
            for i in 1:n_mol
                for j in (i+1):n_mol
                    if euclidean(c[i], c[j]) < 2.0 + 0.25 * rs
                        push!(segments, c[i])
                        push!(segments, c[j])
                    end
                end
            end
            segments
        end

        Sums_mark = @lift(Point2f($obs_index, Sums[$obs_index]))

        # ---- Plotting ----
        cgl = GridLayout(fig[1, 1:2])
        mgl = GridLayout(fig[2, 1:2])

        conf_ax1 = LScene(cgl[1, 1], show_axis = false)
        meshscatter!(conf_ax1, points, markersize = vcat([template_radii for _ in 1:n_mol]...), color = :blue)
        meshscatter!(conf_ax1, points, markersize = vcat([template_radii .+ rs for _ in 1:n_mol]...), color = RGBAf(0.0, 0.5, 0.3, 0.15))

        conf_ax2 = LScene(cgl[1, 2], show_axis = false)
        meshscatter!(conf_ax2, points, markersize = vcat([template_radii .* 0.25 for _ in 1:n_mol]...), color = :blue)
        linesegments!(conf_ax2, connections_for_all, color = :black, linewidth = 2)

        Sums_ax = Axis(mgl[1, 1], title = "Morphometric approach to solvation free energy")
        lines!(Sums_ax, xs, Sums, color = fsolc)
        scatter!(Sums_ax, Sums_mark, color=:black, marker = Rect, markersize = [(1, 2500)])
        scatter!(Sums_ax, Sums_mark, color=:black, marker = Rect, markersize = [(2500, 1)])

        # ---- Finalizing the App Layout ----
        index_display = @lift(DOM.div("Iteration: $($obs_index) / $(length(vis_range))", style="margin: 0px 20px; font-weight: bold;"))
        controls = DOM.div(prev_button, index_display, next_button, style="display: flex; justify-content: center; align-items: center; margin-bottom: 10px;")

        Bonito.record_states(session, DOM.div(controls, fig))
    end
    show(io, MIME"text/html"(), app)
    println(io, """
        </body>
    </html>
    """)
end