# Exotic Clusters - Split into two visualizations
function ply_to_points_manual(filepath::String)::Vector{Point3f}
    points = Point3f[]

    # Open the file and read it line by line
    open(filepath, "r") do f
        # Skip the header until "end_header" is found
        for line in eachline(f)
            if strip(line) == "end_header"
                break
            end
        end

        # Read the data lines
        for line in eachline(f)
            parts = split(line)
            if length(parts) >= 3 # Ensure there are at least 3 values
                # Parse x, y, z and create a Point3f
                x = parse(Float32, parts[1])
                y = parse(Float32, parts[2])
                z = parse(Float32, parts[3])
                push!(points, Point3f(x, y, z))
            end
        end
    end

    return points
end

using WGLMakie, Bonito, MorphoMol, Distances, JLD2, Rotations
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

# Function to create figure with two clusters
function create_clusters_figure(cluster_files, cluster_rs_values; rotations=nothing)
    conf_grad = cgrad(:Dark2_3, 3, categorical=true)
    confc1 = conf_grad[1]
    confc2 = conf_grad[2]
    confc3 = conf_grad[3]
    template_radii = [1.0]

    fig = Figure(size=(800, 800))

    for (idx, (filepath, rs)) in enumerate(zip(cluster_files, cluster_rs_values))
        cluster_gl = GridLayout(fig[idx, 1])
        points = ply_to_points_manual(filepath)

        # Apply rotation if specified (rotations is a vector of (roll, pitch, yaw) tuples in degrees)
        if rotations !== nothing
            roll, pitch, yaw = rotations[idx]
            points = rotate_points_deg(points, roll, pitch, yaw)
        end

        n = length(points)

        # Full cluster view
        cluster_ax = LScene(cluster_gl[1, 1], show_axis = false)
        connections = Vector{Point3f}()
        for i in 1:n
            for j in (i+1):n
                if euclidean(points[i], points[j]) < 2.0 + 0.5 * rs
                    push!(connections, points[i])
                    push!(connections, points[j])
                end
            end
        end
        meshscatter!(cluster_ax, points, markersize = vcat([template_radii for _ in 1:n]...), color = confc1)
        meshscatter!(cluster_ax, points, markersize = vcat([template_radii .+ rs for _ in 1:n]...), color = RGBAf(0.0, 0.5, 0.3, 0.15))

        # Contact graph view
        graph_ax = LScene(cluster_gl[1, 2], show_axis = false)
        meshscatter!(graph_ax, points, markersize = vcat([template_radii .* 0.15 for _ in 1:n]...), color = confc2)
        linesegments!(graph_ax, connections, color = :black, linewidth = 2)
    end

    fig
end

# Generate first HTML (clusters 1 and 2)
open("exotic_clusters_1.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session
        cluster_files = [
            "plys/10_helix_0.455_0.464.ply",
            "plys/12_classic_ufo_0.4625_0.45.ply"
        ]
        cluster_rs = [0.3, 0.4625]

        fig = create_clusters_figure(cluster_files, cluster_rs; rotations=[(-25,35,0), (0,150,0)])
        DOM.div(WGLMakie.WithConfig(fig; use_html_widget=true))
    end
    show(io, MIME"text/html"(), app)
    println(io, """
        </body>
    </html>
    """)
end

# Generate second HTML (clusters 3 and 4)
open("exotic_clusters_2.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    app = App() do session::Session
        cluster_files = [
            "plys/15_bilinski_dodecahedron.ply",
            "plys/17_starbase_0.4_0.45.ply"
        ]
        cluster_rs = [0.225, 0.4]

        fig = create_clusters_figure(cluster_files, cluster_rs; rotations=[(0,0,0), (0,150,0)])
        DOM.div(WGLMakie.WithConfig(fig; use_html_widget=true))
    end
    show(io, MIME"text/html"(), app)
    println(io, """
        </body>
    </html>
    """)
end
