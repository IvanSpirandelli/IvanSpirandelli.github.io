# Interface Visualization for Webpage
# Shows protein dimer interface (top) and random point cloud interface (bottom)
# Each with point cloud on left, interface with filtration slider on right

using WGLMakie, Bonito, JLD2, Rotations
using WGLMakie.Makie.GeometryBasics
WGLMakie.activate!(resize_to=:parent)

"""
    rotate_points_deg(points, roll, pitch, yaw)

Rotate a vector of Point3f by Euler angles (in degrees).
Uses ZYX convention: yaw (Z) -> pitch (Y) -> roll (X).
"""
function rotate_points_deg(points::Vector{Point3f}, roll::Real, pitch::Real, yaw::Real)
    R = RotZYX(deg2rad(yaw), deg2rad(pitch), deg2rad(roll))
    return [Point3f(R * Vec3f(p...)) for p in points]
end

"""
Generate colored mesh from interface surface at a given filtration level
"""
function generate_colored_mesh(vertices, filtration, max_value=Inf)
    faces = [TriangleFace(Int.(e[1])) for e in filtration if length(e[1]) == 3 && e[2] <= max_value]
    if isempty(faces)
        return nothing, Float64[]
    end
    mesh = GeometryBasics.Mesh([Point3f(v) for v in vertices], faces)
    mesh_colors = [e[2] for e in filtration if length(e[1]) == 1]
    return mesh, mesh_colors
end

"""
Prepare rotated meshes for all sampled filtration levels
"""
function prepare_meshes(interface_vertices, interface_filtration, levels, rotation, n_sample=50)
    # Sample filtration levels
    n_sample = min(n_sample, length(levels))
    sampled_indices = round.(Int, range(1, length(levels), length=n_sample))
    sampled_levels = levels[sampled_indices]

    rotated_meshes = []
    rotated_mesh_colors = []

    for lvl in sampled_levels
        mesh, mesh_colors = generate_colored_mesh(interface_vertices, interface_filtration, lvl)
        if mesh !== nothing
            verts_3f = collect(GeometryBasics.coordinates(mesh))
            verts_rotated = rotate_points_deg(verts_3f, rotation...)
            rotated_mesh = GeometryBasics.Mesh(verts_rotated, GeometryBasics.faces(mesh))
            push!(rotated_meshes, rotated_mesh)
            push!(rotated_mesh_colors, mesh_colors)
        end
    end

    return rotated_meshes, rotated_mesh_colors
end

open("4bmg_interface.html", "w") do io
    println(io, """
    <html>
        <head>
            <style>
                body { margin: 0; padding: 0; }
            </style>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session
        # ============= PROTEIN INTERFACE =============
        @load "jld2s/4bmg_dimer_interface.jld2" points color_labels radii interface_vertices interface_filtration levels

        protein_rotation = (-65, 0, 35)  # roll away, yaw
        protein_points_3f = [Point3f(p) for p in points]
        protein_points_rotated = rotate_points_deg(protein_points_3f, protein_rotation...)
        protein_radii = radii
        protein_color_labels = color_labels

        protein_meshes, protein_mesh_colors = prepare_meshes(
            interface_vertices, interface_filtration, levels, protein_rotation, 50
        )
        protein_n_levels = length(protein_meshes)

        # ============= RANDOM INTERFACE =============
        @load "jld2s/random_interface.jld2" points color_labels radii interface_vertices interface_filtration levels

        random_rotation = (0, 0, 0)
        random_points_3f = [Point3f(p) for p in points]
        random_points_rotated = rotate_points_deg(random_points_3f, random_rotation...)
        random_radii = radii
        random_color_labels = color_labels

        random_meshes, random_mesh_colors = prepare_meshes(
            interface_vertices, interface_filtration, levels, random_rotation, 50
        )
        random_n_levels = length(random_meshes)

        # ============= FIGURE LAYOUT =============
        # Two columns: protein (left), random (right)
        # Each column: interface (top), point cloud (bottom)
        fig = Figure(size=(800, 800))

        # --- Left column: Protein ---
        protein_interface_scene = LScene(fig[1, 1], show_axis=false,
            scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

        protein_slider = Slider(1:protein_n_levels, value=protein_n_levels)

        protein_current_mesh = map(protein_slider) do idx
            protein_meshes[idx]
        end
        protein_current_colors = map(protein_slider) do idx
            protein_mesh_colors[idx]
        end
        protein_colorrange = (minimum(protein_mesh_colors[end]), maximum(protein_mesh_colors[end]))

        mesh!(protein_interface_scene, protein_current_mesh,
            color=protein_current_colors,
            colorrange=protein_colorrange,
            colormap=:viridis)

        protein_pc_scene = LScene(fig[2, 1], show_axis=false)
        meshscatter!(protein_pc_scene, protein_points_rotated,
            markersize=protein_radii,
            color=protein_color_labels,
            colormap=:Dark2_4)

        # --- Right column: Random ---
        random_interface_scene = LScene(fig[1, 2], show_axis=false,
            scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

        random_slider = Slider(1:random_n_levels, value=random_n_levels)

        random_current_mesh = map(random_slider) do idx
            random_meshes[idx]
        end
        random_current_colors = map(random_slider) do idx
            random_mesh_colors[idx]
        end
        random_colorrange = (minimum(random_mesh_colors[end]), maximum(random_mesh_colors[end]))

        mesh!(random_interface_scene, random_current_mesh,
            color=random_current_colors,
            colorrange=random_colorrange,
            colormap=:viridis)

        random_pc_scene = LScene(fig[2, 2], show_axis=false)
        meshscatter!(random_pc_scene, random_points_rotated,
            markersize=0.25,
            color=random_color_labels,
            colormap=:Dark2_4)

        # ============= SLIDERS =============
        slider_style = DOM.style("""
            .slider-container {
                display: flex;
                align-items: center;
                padding: 10px 20px;
            }
            .slider-label {
                font-size: 14px;
                margin-right: 10px;
                min-width: 80px;
            }
            input[type="range"] {
                -webkit-appearance: none;
                appearance: none;
                width: 100%;
                height: 8px;
                background: #e0e0e0;
                border-radius: 4px;
                outline: none;
                cursor: pointer;
            }
            input[type="range"]::-webkit-slider-thumb {
                -webkit-appearance: none;
                appearance: none;
                width: 12px;
                height: 12px;
                background: #4063d8;
                border: none;
                border-radius: 50%;
                cursor: pointer;
            }
            input[type="range"]::-moz-range-thumb {
                width: 12px;
                height: 12px;
                background: #4063d8;
                border: none;
                border-radius: 50%;
                cursor: pointer;
            }
        """)

        slider_row = DOM.div(
            slider_style,
            DOM.div(
                DOM.span("Protein:", class="slider-label"),
                protein_slider,
                class="slider-container",
                style="flex: 1;"
            ),
            DOM.div(
                DOM.span("Random:", class="slider-label"),
                random_slider,
                class="slider-container",
                style="flex: 1;"
            ),
            style="display: flex; gap: 20px;"
        )

        Bonito.record_states(session, DOM.div(slider_row, fig))
    end
    show(io, MIME"text/html"(), app)
    println(io, """
        </body>
    </html>
    """)
end

println("Generated 4bmg_interface.html")
