# Green-Ampt Problem Solver in Julia


# Parameters
L = .5          # Length of soil column (meters)
Nx = 100         # Number of spatial grid points
dx = L / Nx      # Grid spacing
dt = 0.01       # Smaller time step
T_final = 1.0    # Final simulation time (seconds)

# Soil properties (you can adjust these)
Ks = 2.06829e-5 #Tu2       # Saturated hydraulic conductivity (m/s)
θs = 0.4971    #0.432         # Saturated water content
θr = 0.05      # Residual water content
α = 7.242 / 10     # van Genuchten parameter
n = 1.06062    # van Genuchten parameter

#Umrechnung in 1/cm in 1/m

# Initial conditions
h0 = 0.0         # Initial water pressure head

θ = zeros(Nx)
fill!(θ, θs)  # Initialize theta values to saturated water content

# Set boundary conditions
h_left = 1.0     # Left boundary: constant water pressure head
h_right = 1.0   # Right boundary: constant water pressure head

# Time-stepping loop
for t in 0:dt:T_final
    for i in 2:Nx-1
        # Discretized Green-Ampt equation
        i_next = Int64(min(Nx, i + 1))
        if i_next <= Nx
            h_next = h0 + (θ[i] - θr) * (h_left - h0) / (θ[i_next] - θr)
            if h_next > 0.0
                θ[i_next] = θr + (θ[i] - θr) * (1 + (α * h_next)^n) ^ (-1 / n)
            end
        end
    end
    # Update the boundary conditions at every time step
    h_left = 0.5 * sin(π * (t / T_final))
    h_right = 0.5 * sin(π * (t / T_final))
end

# Print results
println("Final water content profile:")
for i in 1:Nx
    println("x = ", i*dx, "\tθ = ", θ[i])
end

plot(θ,yaxis=:log10,label="Θ",title="Tu2 Results")

#plot(θ,yaxis=:log2,label="Θ")

# #heatmap(θ, c=:matter, title="Green-Ampt Results", xaxis="", yaxis="")
# plot(xaxis=1:length(θ),x=θ, 
#     yaxis=:log,
#     c=:matter, title="Green-Ampt Results")
# boxplot(θ)


##############
function plot_retention_curves(df)
    # Extract the relevant columns from the DataFrame
    ths = df[:, :ths]
    thr = df[:, :thr]
    alpha = df[:, :alpha]
    npar = df[:, :npar]
    nm = df[:, :Name]
    mpar = df[:, :mpar]
    ksat = df[:, :ksat]
    # tort = df[:, :tort]

    # Create an empty plot
    p1 = plot(title="Retention Curves")

    # Plot each retention curve
    for i in 1:size(df, 1)
        θ = zeros(Nx)
        fill!(θ, ths[i])
        h_left = 1.0
        h_right = 1.0

        for t in 0:dt:T_final
            for j in 2:Nx-1
                i_next = Int64(min(Nx, j + 1))
                if i_next <= Nx
                    h_next = h0 + (θ[j] - thr[i]) * (h_left - h0) / (θ[i_next] - thr[i])
                    if h_next > 0.0
                        θ[i_next] = thr[i] + (θ[j] - thr[i]) * (1 + (alpha[i] * h_next)^npar[i]) ^ (-1 / npar[i])
                    end
                end
            end
            h_left = 0.5 * sin(π * (t / T_final))
            h_right = 0.5 * sin(π * (t / T_final))
        end

        lab = nm[i]

        plot!(θ, yaxis=:identity, 
        legend = :outertopright,
        legendcolumns = 2,
            label="$lab")
    end

    return p1
    
end

nom = hsub|>dropmissing
sort!(nom,:Name)
plot_retention_curves(nom)


function theta(Psi, kf, alpha, n, m, l)
    numerator = (1 - (alpha * abs(Psi))^(n-1) * (1 + (alpha * abs(Psi))^n)^(-m))^2
    denominator = (1 + (alpha * abs(Psi))^n)^(m*l)
    return kf * (numerator / denominator)
end

function pc2(df)
    # Extract the relevant columns from the DataFrame
    ths = df[:, :ths]
    thr = df[:, :thr]
    alpha = df[:, :alpha]
    npar = df[:, :npar]
    nm = df[:, :Name]
    mpar = df[:, :mpar]
    ksat = df[:, :ksat]
    tort = df[:, :tort]

    # Parameters
    #L = .5          # Length of soil column (meters)
    Nx = 100        # Number of spatial grid points
    #dx = L / Nx      # Grid spacing
    dt = 0.1       # time step
    T_final = 1.0  # Final simulation time

    # Create an empty plot
    p1 = plot(title="Retention Curves")
    tdf = []
    # Plot each retention curve
    for i in 1:size(df, 1)
        θ = zeros(Nx)
        fill!(θ, ths[i])
        h_left = 0.0
        h_right = 0.0

        for t in 0:dt:T_final
            for j in 2:Nx-1
                i_next = Int64(min(Nx, j + 1))
                if i_next <= Nx
                    h_next = h0 + (θ[j] - thr[i]) * (h_left - h0) / (θ[i_next] - thr[i])
                    if h_next > 0.0
                        θ[i_next] = theta(h_next, ksat[i], alpha[i], npar[i], mpar[i], tort[i])
                    end
                end
            end
            h_left = 0.5 * sin(π * (t / T_final))
            h_right = 0.5 * sin(π * (t / T_final))
        end

        lab = string.(nm[i])

        plot!(θ, yaxis=:log10, 
        legend = :outertopright,
        legendcolumns = 2,
            label="$lab")
        push!(tdf,
        DataFrame(lab=θ))        
    end



    display(p1)
    return tdf
end

plotly()
ndf = pc2(nom)
plot_retention_curves(nom)
