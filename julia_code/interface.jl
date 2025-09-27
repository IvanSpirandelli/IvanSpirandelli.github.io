
# 1. Import the necessary libraries
using WGLMakie, Makie, FileIO, JSServe

# Set WGLMakie as the active backend
WGLMakie.activate!()


# 2. Create the figure and plot
fig = Figure()

# 3. Create an interactive element (a slider)
slider = Makie.Slider(fig[2, 1])
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


# 3. Create a JSServe Page object configured for offline export
page = JSServe.Page(
    exportable=true, # an exportable page can be written to a file
    offline=true     # this inlines all JS/CSS dependencies
)

# 4. Insert your Makie figure into the page's content
# This is an alternative and more direct way than the previous DOM.div
page(fig)

open("interactive_plot.html", "w") do io
    println(io, JSServe.page_html(page))
end

println("Successfully saved interactive plot to interactive_plot.html")