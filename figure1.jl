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

function NegFeedDelOsc(vars::Vector{Float64}, params::Vector{Float64})
    ydelay=params[1]
    dx = beta1*k1^n / (k1^n + ydelay^n) - delta1*vars[1]
    dy = beta2*vars[1]^m / (k2^m+ vars[1]^m) - delta2*vars[2]
    return dx, dy
end

const beta1  = 10.0
const delta1 = 1.0
const k1     = 0.5
const n      = 2
const beta2  = 10.0
const delta2 = 1.0
const k2     = 5.0
const m      = 2.0

h = 0.1
params=[0.0]
times = collect(range(0,step=h,stop=50))

### PARAMETERS FOR THE FIRST PANEL
tau1 = 0.7
stepsd1   = round(Int64,tau1/h)
delayedy1 = zeros(stepsd1)
state1 = [1.0, 1.0]
X1 = zeros(length(times))
Y1 = zeros(length(times))

### PARAMETERS FOR THE SECOND PANEL
tau2 = 2.0
stepsd2   = round(Int64,tau2/h)
delayedy2 = zeros(stepsd2)
state2 = [1.0, 1.0]
X2 = zeros(length(times))
Y2 = zeros(length(times))

### PARAMETERS FOR THE THIRD PANEL
tau3 = 4.0
stepsd3   = round(Int64,tau3/h)
delayedy3 = zeros(stepsd3)
state3 = [1.0, 1.0]
X3 = zeros(length(times))
Y3 = zeros(length(times))

### INTEGRATION
for t in eachindex(times)   
    ## Integrate first condition
    params[1] = delayedy1[end]
    state1 .= RK4(NegFeedDelOsc, state1, params, h)
    X1[t] = state1[1]
    Y1[t] = state1[2]
    delayedy1[2:end] .= delayedy1[1:end-1]
    delayedy1[1]=Y1[t]

    ## Integrate second condition
    params[1] = delayedy2[end]
    state2 .= RK4(NegFeedDelOsc, state2, params, h)
    X2[t] = state2[1]
    Y2[t] = state2[2]
    delayedy2[2:end] .= delayedy2[1:end-1]
    delayedy2[1]=Y2[t]

    ## Integrate third condition
    params[1] = delayedy3[end]
    state3 .= RK4(NegFeedDelOsc, state3, params, h)
    X3[t] = state3[1]
    Y3[t] = state3[2]
    delayedy3[2:end] .= delayedy3[1:end-1]
    delayedy3[1]=Y3[t]
end

### PLOTTING
PyPlot.close("all")
fig, ax = PyPlot.subplots(3,1,figsize=(7,10),sharex=true)

## Plot first panel
ax[1].plot(times, X1, label="x", linewidth=3,color=(0.0,0.4,1.0))
ax[1].plot(times, Y1, label="y", linewidth=3,color=(1.0,0.4,0.))
ax[1].set_yticks([])
ax[1].set_xticks(range(0,stop=times[end],step=5))
ax[1].set_xticklabels([0,"",10,"",20,"",30,"",40,"",50])
ax[1].tick_params(axis="both", which="major", labelsize=30)
ax[1].spines["right"].set_visible(false)
ax[1].spines["top"].set_visible(false)
rng = round(Int64, 5.0/h)
xmax = maximum(X1[rng:end])
xmin = minimum(X1[rng:end])
ymax = maximum(Y1[rng:end])
ymin = minimum(Y1[rng:end])
ax[1].set_xlim(times[rng],times[end])
ax[1].set_ylim(ymin-0.1*xmax, xmax+0.1*xmax)
ax[1].legend(fontsize=20, bbox_to_anchor=[1.2,1.25])

## Plot second panel
ax[2].plot(times, X2, label="x", linewidth=3,color=(0.0,0.4,1.0))
ax[2].plot(times, Y2, label="y", linewidth=3,color=(1.0,0.4,0.))
ax[2].set_yticks([])
ax[2].set_xticks(range(0,stop=times[end],step=5))
ax[2].set_xticklabels([0,"",10,"",20,"",30,"",40,"",50])
ax[2].tick_params(axis="both", which="major", labelsize=30)
ax[2].spines["right"].set_visible(false)
ax[2].spines["top"].set_visible(false)
rng = round(Int64, 5.0/h)
xmax = maximum(X2[rng:end])
xmin = minimum(X2[rng:end])
ymax = maximum(Y2[rng:end])
ymin = minimum(Y2[rng:end])
ax[2].set_xlim(times[rng],times[end])
ax[2].set_ylim(ymin-0.1*xmax, xmax+0.1*xmax)

## Plot third panel
ax[3].plot(times, X3, label="x", linewidth=3,color=(0.0,0.4,1.0))
ax[3].plot(times, Y3, label="y", linewidth=3,color=(1.0,0.4,0.))
ax[3].set_yticks([])
ax[3].set_xticks(range(0,stop=times[end],step=5))
ax[3].set_xticklabels([0,"",10,"",20,"",30,"",40,"",50])
ax[3].tick_params(axis="both", which="major", labelsize=30)
ax[3].spines["right"].set_visible(false)
ax[3].spines["top"].set_visible(false)
rng = round(Int64, 5.0/h)
xmax = maximum(X3[rng:end])
xmin = minimum(X3[rng:end])
ymax = maximum(Y3[rng:end])
ymin = minimum(Y3[rng:end])
ax[3].set_xlim(times[rng],times[end])
ax[3].set_ylim(ymin-0.1*xmax, xmax+0.1*xmax)

plt.tight_layout()
plt.savefig("figures/figure1/time_trace.svg", format="svg")

gcf()