using PyPlot

function RK4( f::Function, vars::Vector{Float64}, params::Vector{Float64}
    , h::Float64)

    k1 = f(vars, params)
    vars_k2 = vars .+ 0.5 * h .* k1
    k2 = f(vars_k2, params)
    vars_k3 = vars .+ 0.5 * h .* k2
    k3 = f(vars_k3, params)
    vars_k4 = vars .+ h .* k3
    k4 = f(vars_k4, params)

    slope = (k1 .+ 2.0.*k2 .+ 2.0.*k3 .+ k4) ./ 6.0
    @. vars = vars + (h * slope)
    return vars
end

function erkmodel(vars::Vector{Float64}, params::Vector{Float64})
    fgf = params[1]
    R   = vars[1]
    M   = vars[2]
    E   = vars[3]
    dR  = βR1 * (1/(1 +(E/kR1)^n))*(fgf/(kR2 + fgf)) - (βR2*R)/(kR3 + R)
    dM  = βM1 * (R^m)/(kM1^m + R^m) - (βM2*M)/(kM2 + M)
    dE  = βE1 * M^2/(kE1^2 + M^2) - (βE2*E)/(kE2 + E) 
    return dR, dM, dE
end

const kR1 = 9.0
const kR2 = 0.1
const kR3 = 8.0
const βR1 = 2.0
const βR2 = 0.5
const n   = 2

const m   = 6
const kM1 = 10.0
const kM2 = 10.0
const βM1 = 1.0
const βM2 = 0.75

const kE1 = 15.0
const kE2 = 4.0
const βE1 = 1.0
const βE2 = 0.75

params = [2.0]
state  = [0.0, 0.0 ,0.0]
h = 0.5

### FIRST SIMULATION
# Run a simulation to find the ERK minimum within a cycle of oscillation
times = collect(range(0,step=h,stop=5000))
R = zeros(length(times))
M = zeros(length(times))
E = zeros(length(times))
for t in eachindex(times)   
    state .= RK4(erkmodel, state, params, h)
    R[t] = state[1]
    M[t] = state[2]
    E[t] = state[3]
end

### SECOND SIMULATION 
# Select the values for the start of the plotting simulation
# This step also helps getting rid of transitory dynamics at the start of the simulation
idx   = argmin(E[round(Int64, length(E)/2):end])
state = [R[round(Int64, length(E)/2):end][idx], M[round(Int64, length(E)/2):end][idx] ,E[round(Int64, length(E)/2):end][idx]]
times = collect(range(0,step=h,stop=2500))
R = zeros(length(times))
M = zeros(length(times))
E = zeros(length(times))

# Integration loop
for t in eachindex(times)   
    R[t] = state[1]
    M[t] = state[2]
    E[t] = state[3]
    state .= RK4(erkmodel, state, params, h)
end

### PLOTTING
PyPlot.close("all")
fig, ax = plt.subplots(figsize=(9,4))

# Convert time to minutes
times = times./60 
ax.plot(times, R, label="Raf", linestyle="--", linewidth=2)
ax.plot(times, M, label="MEK", linestyle="--", linewidth=2)
ax.plot(times, E, label="ERK", linewidth=3, color="green")
ax.tick_params(axis="x",labelsize=20)
ax.set_xticks([0,5,10,15,20,25,30,35,40])
ax.set_xticklabels([0,"",10,"",20,"",30, "",40])
ax.set_yticks([])
ax.axes.get_yaxis().set_visible(false)
plt.legend(fontsize=20, bbox_to_anchor=[1.27,1.045])
plt.tight_layout()
plt.savefig("figures/figure3/time_trace.svg", format="svg")
gcf()
