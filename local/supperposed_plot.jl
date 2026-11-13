using CairoMakie
using JLD2
using LaTeXStrings

include("FieldLines.jl")

t_f = 0.01
iota = 4

deltas = [0.05, 0.1, 0.2, 0.3, 0.4]
# delta = 0.5
# iotas = [2.2, 2.4, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8]
iotas = [4.0]

for delta in deltas
    for iota in iotas
        δ1 = delta
        δ2 = delta
        xₛ = 1.0

        function B(X, x, p, t, iota=iota)
            m1 = 1
            n1 = -1
            m2 = 1
            n2 = 1

            X[1] = δ1*x[1]*(2-x[1])*0.25*m1*sin(m1*x[2]-n1*t) +
                δ2*x[1]*(2-x[1])*0.25*m2*sin(m2*x[2]-n2*t)

            X[2] = iota*(x[1]-xₛ) +
                δ1*(2-2*x[1])*0.25*cos(m1*x[2]-n1*t) +
                δ2*(2-2*x[1])*0.25*cos(m2*x[2]-n2*t)
        end

        poindata = FieldLines.construct_poincare(B, [0.0, 2.0], [0, 2π])

        k_values = [4.0, 6.0, 8.0, 10.0, 12.0]
        colors = [:red, :orange, :green, :purple, :blue]

        hs = Figure(
            size=(2200, 1400),
            fontsize=45,
            figure_padding=(5, 5, 35, 45)
        )

        # Label(
        #     hs[1, 1],
        #     "δ₁ = δ₂ = $(δ1), m₁ = 1, m₂ = 1, n₁ = -n₂ = 1",
        #     fontsize=50,
        #     font=:bold,
        #     padding=(0, 0, 5, 5)
        # )


        ax_p = Axis(
            hs[2, 1],
            xlabel="x",
            ylabel="y",
            xlabelsize=60,
            ylabelsize=60,
            xticklabelsize=55,
            yticklabelsize=55,
            xaxisposition=:bottom,
            yaxisposition=:left
        )

        xlims!(ax_p, 0, 2)
        ylims!(ax_p, 0, 2π)

        ax_u = Axis(
            hs[2, 1],
            xlabel="",
            ylabel="T(x)",
            yaxisposition=:right,
            xaxisposition=:top,
            xticklabelsvisible=false,
            xticksvisible=false,
            xgridvisible=false,
            ygridvisible=false,
            xlabelsize=60,
            ylabelsize=60,
            yticklabelsize=55,
            backgroundcolor=:transparent
        )

        xlims!(ax_u, 0, 2)

        # Make the plot cell occupy almost everything
        rowgap!(hs.layout, 2)
        colgap!(hs.layout, 2)

        rowsize!(hs.layout, 1, Fixed(55))       # title height
        rowsize!(hs.layout, 2, Relative(1.0))   # plot gets the rest
        colsize!(hs.layout, 1, Relative(1.0))   # full width

        scatter!(
            ax_p,
            poindata.ψ[:],
            poindata.θ[:],
            markersize=3.0,
            color=:black,
            alpha=0.35
        )

        legend_elements = []
        legend_labels = []

        for (i, exp_k) in enumerate(k_values)

            infile = "C:/Users/mirei/OneDrive - epfl.ch/MA3/SdS/local/num_results_loc/FaADE/2I_k$(exp_k)_tf$(t_f)_d$(δ1)_$(δ2)_res201_m1_dT_iota$(iota).jld2"

            x_mid = JLD2.load(infile, "x_mid")
            u_mid = JLD2.load(infile, "u_mid")

            lineplot = lines!(
                ax_u,
                x_mid,
                u_mid,
                color=colors[i],
                linewidth=5
            )

            push!(legend_elements, lineplot)
            push!(legend_labels, L"\kappa_\perp/\kappa_\parallel = 10^{%$(-Int(exp_k))}")
        end

        hidexdecorations!(ax_u, grid=false)

        Legend(
            hs[2, 1], 
            legend_elements,
            legend_labels,
            title="Parallel transport",
            titlesize=60,
            labelsize=60,
            framevisible=true,
            halign=:right,
            valign=:bottom,
            margin=(10, 20, 20, 10)
        )

        display(hs)

        save(
            "C:/Users/mirei/OneDrive - epfl.ch/MA3/SdS/local/Figures/FaADE/2I_tf$(t_f)_d$(δ1)_$(δ2)_m1_iota$(iota)_superposed.png",
            hs
        )
    end
end