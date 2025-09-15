
# 1. Import the necessary libraries
using WGLMakie, Makie, FileIO

# Set WGLMakie as the active backend
WGLMakie.activate!()


# 2. Create the figure and plot
fig = Figure()

# 3. Create an interactive element (a slider)
slider = Slider(
    fig[2, 1]
)
slider_value = slider.value # This is an "Observable"

ax = Axis3(fig[1, 1],
    title = lift(val -> "Surface Plot with param = $(round(val, digits=2))", slider_value)
)
# Define the surface function using the slider's value
xs = -10:0.5:10
ys = -10:0.5:10
z_data = lift(slider_value) do val
    [sin(sqrt(x^2 + y^2) - val) for x in xs, y in ys]
end

surface!(ax, xs, ys, z_data)

# 4. Add the slider to the layout
fig[2, 1] = slider

# Display the plot locally to test it (optional)
# display(fig)

# 5. Save the figure as a self-contained HTML file
# This is the crucial step!
save("interactive_plot.html", fig)

println("Successfully saved interactive plot to interactive_plot.html")