using WGLMakie
using FileIO
using Observables
using Bonito

WGLMakie.activate!()

open("index.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    # Then, you can just inline plots or whatever you want :)
    # Of course it would make more sense to put this into a single app
    app = App() do session::Session
        n = 10
        index_slider = Slider(1:n)
        volume = rand(n, n, n)
        slice = map(index_slider) do idx
            return volume[:, :, idx]
        end
        fig = Figure()
        ax, cplot = contour(fig[1, 1], volume)
        rectplot = linesegments!(ax, Rect(-1, -1, 12, 12), linewidth=2, color=:red)
        on(index_slider) do idx
            translate!(rectplot, 0,0, idx)
        end
        heatmap(fig[1, 2], slice)
        slider = DOM.div("z-index: ", index_slider, index_slider.value)
        return Bonito.record_states(session, DOM.div(slider, fig))
    end
    show(io, MIME"text/html"(), app)
    # or anything else from Bonito, or that can be displayed as html:
    println(io, """
        </body>
    </html>
    """)
end