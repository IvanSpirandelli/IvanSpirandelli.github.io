using WGLMakie, Bonito, FileIO, Observables, JLD2, MorphoMol
WGLMakie.activate!()

open("interactive_plot.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session
    function get_two_subunit_phase_diagram()

        conf_grad = cgrad(:Dark2_3, 3, categorical=true)
        confc1 = conf_grad[1]
        confc2 = conf_grad[2]
        confc3 = conf_grad[3]
        confmap = [confc1, confc3, confc2]

        @load "jld2s/weighted_recalculation_aggregated.jld2" all_persistence_weights min_states min_measures_normalized min_thetas pinput
        all_persistence_weights_o = deepcopy(all_persistence_weights)
        min_states_o = deepcopy(min_states)
        min_measures_normalized_o = deepcopy(min_measures_normalized)
        min_thetas_o = deepcopy(min_thetas);

        template_radii = pinput["template_radii"]

        @load "jld2s/weighted_aggregated.jld2" all_persistence_weights min_states min_measures_normalized min_thetas pinput
        all_persistence_weights = [all_persistence_weights; all_persistence_weights_o]
        min_states = [min_states; min_states_o]
        min_measures_normalized = [min_measures_normalized; min_measures_normalized_o]
        min_thetas = [min_thetas; min_thetas_o];

        exp_template_centers = MorphoMol.EXPERIMENTAL_ASSEMBLY["6r7m"]["template_centers"]
        exp_template_radii = MorphoMol.EXPERIMENTAL_ASSEMBLY["6r7m"]["template_radii"]
        exp_state = MorphoMol.EXPERIMENTAL_ASSEMBLY["6r7m"]["state"]

        points = MorphoMol.get_point_vector_realization(exp_state, exp_template_centers)
        radii = vcat([exp_template_radii for i in 1:2]...)
        dgms = MorphoMol.Energies.get_weighted_alpha_shape_persistence_diagram(points, radii)
        p0 = MorphoMol.Energies.get_total_persistence(dgms[1])
        p1 = MorphoMol.Energies.get_total_persistence(dgms[2])
        p2 = MorphoMol.Energies.get_total_persistence(dgms[3])

        exp_measure = [p0, p1, p2]

        p1s = -2.0:0.1:1.0
        p2s = -2.0:0.1:1.0

        dist_to_exp = [Inf for i in 1:length(p1s), j in 1:length(p2s)]

        for i in 1:length(p1s)
            for j in 1:length(p2s)
                p1 = p1s[i]
                p2 = p2s[j]
                exp_eval = sum([1.0, p1, p2] .* exp_measure)
                min_eval = minimum([sum([1.0, p1, p2] .* min_measure) for min_measure in min_measures_normalized])
                dist_to_exp[i, j] = log(exp_eval - min_eval)
            end
        end
        cb_range = (minimum(dist_to_exp), maximum(dist_to_exp))

        n_mol = pinput["n_mol"]
        n_atoms_per_mol = length(pinput["template_radii"])

        coma = :viridis

        pad = 1

        fx = Figure(fontsize = 16)
        hm_gl = GridLayout(fx[1:2,1:2])

        ax_2_phase = Axis(hm_gl[1:2,1:2], aspect = 1, xlabel = L"\lambda_1", ylabel = L"\lambda_2", ylabelrotation = 0.0)
        hm = heatmap!(ax_2_phase, p1s, p2s, dist_to_exp, colormap = coma)

        p1s = -2.0:0.1:1.0
        p2s = -2.0:0.1:1.0

        minimizers = [0 for i in 1:length(p1s), j in 1:length(p2s)]

        for i in 1:length(p1s)
            for j in 1:length(p2s)
                p1 = p1s[i]
                p2 = p2s[j]
                min_dex = argmin([sum([1.0, p1, p2] .* min_measure) for min_measure in min_measures_normalized])
                minimizers[i, j] = min_dex
            end
        end
        uniques = unique(minimizers)

        aggregation_map = Dict(u => u for u in uniques)
        if pinput["n_mol"] <= 6
            template_centers = pinput["template_centers"]
            aggregation_threshold = 2.0

            for i in 1:length(uniques)
                x1 = min_states[uniques[i]]
                for j in i+1:length(uniques)
                    x_comp = min_states[uniques[j]]
                    theta = MorphoMol.get_configuration_distance(
                        MorphoMol.convert_flat_state_to_tuples(x1), 
                        MorphoMol.convert_flat_state_to_tuples(x_comp), 
                        template_centers, template_centers, 
                        MorphoMol.get_theta_of_pair
                        )
                    if theta < aggregation_threshold
                        aggregation_map[uniques[j]] = uniques[i]
                    end
                end
            end
            for _ in 1:10
                for k in keys(aggregation_map)
                    aggregation_map[k] = aggregation_map[aggregation_map[k]]
                end
            end
            aggregation_map
        end

        aggregation_map[2799] = 82
        aggregation_map[4619] = 82
        aggregation_map[3828] = 82
        aggregation_map[2008] = 82
        aggregation_map[41] = 82

        aggregation_map[1913] = 549
        aggregation_map[1790] = 549

        aggregation_map[548] = 3518
        aggregation_map[718] = 3518

        for k in keys(aggregation_map)
            aggregation_map[k] = aggregation_map[aggregation_map[k]]
        end

        aggregated_minimizers = [aggregation_map[minimizers[i, j]] for i in 1:length(p1s), j in 1:length(p2s)]
        aggregated_uniques = unique(aggregated_minimizers)

        function flip_state(state)
            [state[7:12]; state[1:6]]
        end

        n_mol = 2
        n_atoms_per_mol = length(pinput["template_radii"])

        aggregated_uniques = [549, 998, 3518, 82]
        states_to_flip = []
        configurations = []
        radii = vcat([pinput["template_radii"] for i in 1:pinput["n_mol"]]...)

        configurations = []
        radii = vcat([pinput["template_radii"] for i in 1:pinput["n_mol"]]...)

        for state in min_states[aggregated_uniques]
            conf = MorphoMol.get_point3f_realization(state, pinput["template_centers"])
            com = sum(conf) / length(conf)
            conf_centered = [p - com for p in conf]
            push!(configurations, conf_centered)
        end

        remaps = Dict(aggregated_uniques[i] => findfirst(x -> x == aggregated_uniques[i], aggregated_uniques) for i in 1:length(aggregated_uniques))
        remapped_minimizer = [remaps[aggregated_minimizers[i, j]] for i in 1:length(p1s), j in 1:length(p2s)];

        function draw_borders(i,j, p1s, p2s, remapped_minimizer, ax)
            offset = 0.05
            if (i > 1 && remapped_minimizer[i-1, j] != remapped_minimizer[i, j])
                lines!(ax, [p1s[i-1] + offset, p1s[i] - offset], [p2s[j] - offset, p2s[j] + offset], color = :white, linewidth = 1.5)
            end
            if (j > 1 && remapped_minimizer[i, j-1] != remapped_minimizer[i, j])
                lines!(ax, [p1s[i] - offset, p1s[i] + offset], [p2s[j-1] + offset, p2s[j] - offset], color = :white, linewidth = 1.5)
            end
        end

        for i in 1:size(remapped_minimizer)[1]
            for j in 1:size(remapped_minimizer)[2]
                draw_borders(i,j, p1s, p2s, remapped_minimizer, ax_2_phase)
            end
        end

        COLORS = [(0.8, 0.12, 0.5), (0.98, 0.54, 0.14), (0.06, 0.3, 0.36), (0.34, 0.54, 0.98)]
        colors = [i for i in 1:2 for _ in 1:length(exp_template_radii)]

        conf_one_gl = GridLayout(fx[1, 3])
        conf_one_ax = LScene(conf_one_gl[1, 1], show_axis = false)
        meshscatter!(conf_one_ax, configurations[4], markersize = vcat([template_radii for _ in 1:n_mol]...), color = colors, colormap = confmap)

        conf_two_gl = GridLayout(fx[1, 4])
        conf_two_ax = LScene(conf_two_gl[1, 1], show_axis = false)
        meshscatter!(conf_two_ax, configurations[2], markersize = vcat([template_radii for _ in 1:n_mol]...), color = colors, colormap = confmap)

        conf_three_gl = GridLayout(fx[2, 3])
        conf_three_ax = LScene(conf_three_gl[1, 1], show_axis = false)
        meshscatter!(conf_three_ax, configurations[3], markersize = vcat([template_radii for _ in 1:n_mol]...), color = colors, colormap = confmap)

        conf_four_gl = GridLayout(fx[2, 4])
        conf_four_ax = LScene(conf_four_gl[1, 1], show_axis = false)
        meshscatter!(conf_four_ax, configurations[1], markersize = vcat([template_radii for _ in 1:n_mol]...), color = colors, colormap = confmap)

        for (label, layout) in zip(["A", "B", "C", "D"], [conf_one_gl, conf_two_gl, conf_three_gl, conf_four_gl])
            Label(layout[1, 1, TopLeft()], label,
                fontsize = 12,
                font = :bold,
                padding = (0, pad, pad, 0),
                halign = :right)
        end

        text!(ax_2_phase, -1.3, 0.4, text = "A", align = (:center, :center), 
            fontsize = 18, font = :bold, color = :white)
        text!(ax_2_phase, 0.5, -0.5, text = "B", align = (:center, :center), 
            fontsize = 18, font = :bold, color = :white)
        text!(ax_2_phase, -1.0, -1.0, text = "C", align = (:center, :center), 
            fontsize = 18, font = :bold, color = :white)
        text!(ax_2_phase, -1.8, -1.8, text = "D", align = (:center, :center), 
            fontsize = 18, font = :bold, color = :white)
        text!(ax_2_phase, -0.75, -1.95, text = "D", align = (:center, :center), 
            fontsize = 18, font = :bold, color = :white)

        fx
    end
    # n = 10
    # index_slider = Slider(1:n)
    # volume = rand(n, n, n)
    # slice = map(index_slider) do idx
    #     return volume[:, :, idx]
    # end
    # fig = Figure()
    # ax, cplot = contour(fig[1, 1], volume)
    # rectplot = linesegments!(ax, Rect(-1, -1, 12, 12), linewidth=2, color=:red)
    # on(index_slider) do idx
    #     translate!(rectplot, 0,0,idx)
    # end
    # heatmap(fig[1, 2], slice)
    # slider = DOM.div("z-index: ", index_slider, index_slider.value)
    #return Bonito.record_states(session, DOM.div(slider, fig))
    fig = get_two_subunit_phase_diagram()
    DOM.div(WGLMakie.WithConfig(fig; use_html_widget=true))
    end
    show(io, MIME"text/html"(), app)
    # or anything else from Bonito, or that can be displayed as html:
    println(io, """
        </body>
    </html>
    """)
end

