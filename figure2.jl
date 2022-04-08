using PyCall

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

function dll1model(vars::Vector{Float64}, params::Vector{Float64})
    HP = vars[1]
    HR = vars[2]
    HF = vars[3]
    MP = vars[4]
    MR = vars[5]
    MF = vars[6]
    DP = vars[7]
    DR = vars[8]

    if params[5]==1.0
        dHP = k1*HR - k2*HP*HF - k3*HP
        dHR = k4/(1+HP^2) - k5*HR
        dHF = k6/(1+HP^2) - k2*HP*HF - k7*HF
        dMP = k8*MR - k9*MP*MF - k10*MP
        dMR = k11/(1+HP^2) - k12*MR
        dMF = k13*(1+HP^2) - k9*MP*MF - k14*MF
        dDP = k15*DR - k16*DP
        dDR = k17/(1+HP^2)*(k19 + (MP^2)/k20) - k18*DR
    else
        dHP = - (HP - params[1])/tau1
        dHR = 0.0
        dHF = 0.0
        dMP = - (MP - params[2])/tau2
        dMR = 0.0
        dMF = 0.0
        max = params[3]*1.2 - params[4]
        dp  = DP - params[4]
        dDP = dp*(1 - dp/max)/tau3
        dDR = 0.0
    end
    return dHP, dHR, dHF, dMP, dMR, dMF, dDP, dDR
end

# Model constants
const k1 = 0.15
const k2 = 0.011
const k3 = 0.0155
const k4 = 0.25
const k5 = 0.014
const k6 = 10.0
const k7 = 0.15
const k8 = 0.3
const k9 = 0.01
const k10 = 0.007
const k11 = 0.025
const k12 = 0.00385
const k13 = 2.0
const k14 = 0.015
const k15 = 0.05
const k16 = 0.0085
const k17 = 0.01
const k18 = 0.01155
const k19 = 10.0
const k20 = 5.0

# Differentiation cartoon constants
const tau1 = 20.0
const tau2 = 50.0
const tau3 = 20.0
const timedif = 780.0

### FIRST SIMULATION
# Run a simulation to find the ERK minimum within a cycle of oscillation
params=[0.0, 0.0, 0.0, 0.0, true]
state =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
h = 0.001
times = collect(range(0,step=h,stop=1000))
HP = zeros(length(times))
HR = zeros(length(times))
HF = zeros(length(times))
MP = zeros(length(times))
MR = zeros(length(times))
MF = zeros(length(times))
DP = zeros(length(times))
DR = zeros(length(times))

for t in eachindex(times)   
    state .= RK4(dll1model, state, params, h)
    HP[t] = state[1]
    HR[t] = state[2]
    HF[t] = state[3]
    MP[t] = state[4]
    MR[t] = state[5]
    MF[t] = state[6]
    DP[t] = state[7]
    DR[t] = state[8]
end

### SECOND SIMULATION 
# Select the values for the start of the plotting simulation
# This step also helps getting rid of transitory dynamics at the start of the simulation
params=[0.0, 0.0, 0.0, 0.0, true]
idx1 = round(Int64, length(times)/2)
idx2 = length(times)
params[1] = minimum(HP[idx1:idx2])
params[2] = maximum(MP[idx1:idx2])
params[3] = maximum(DP[idx1:idx2])
params[4] = minimum(DP[idx1:idx2])-0.1
times = collect(range(0,step=h,stop=timedif+400))
HP = zeros(length(times))
HR = zeros(length(times))
HF = zeros(length(times))
MP = zeros(length(times))
MR = zeros(length(times))
MF = zeros(length(times))
DP = zeros(length(times))
DR = zeros(length(times))

for t in eachindex(times)
    #check if we are in the differentiation period
    if times[t] > timedif
        params[5] = false
    end
    state .= RK4(dll1model, state, params, h)
    HP[t] = state[1]
    HR[t] = state[2]
    HF[t] = state[3]
    MP[t] = state[4]
    MR[t] = state[5]
    MF[t] = state[6]
    DP[t] = state[7]
    DR[t] = state[8]
end

### PLOTTING
PyPlot.close("all")
fig, ax = PyPlot.subplots(figsize=(10,4))

ax1 = ax.twinx()
ax1.plot(times./60, HP, color=(0.0,100.0/255,1.0,1.0), linewidth=3, label="Hes1")
ax.fill_between([timedif./60, times[end]./60], [maximum(MP)*1.05, maximum(MP)*1.05], [minimum(MP)*0.95, minimum(MP)*0.95], color=(0.8, 0.8, 0.8), alpha=1.0)
ax.plot(times./60, MP, color=(0.3,0.6,0.3,1.0), linewidth=3, label="MyoD")
ax1.plot([], [], color=(0.3,0.6,0.3,1.0), linewidth=3, label="MyoD")
ax1.plot(times./60, DP, color=(1.0,0.3,0.3,1.0), linewidth=3, label="Dll1")

ax1.axes.get_yaxis().set_visible(false)
ax.set_xticks([0,2,4,6,8,10,12])
ax.set_xticklabels([0,"",4,"",8,"",12])
ax.tick_params(axis="both", which="major", labelsize=20)
ax.set_yticks([])
ax.plot(timedif.*ones(1000)./60, range(minimum(MP)*0.95,length=1000, stop=maximum(MP)*1.05), color="black", linestyle="--", linewidth=3)
ax1.set_ylim(minimum(HP)*0.9, maximum(HP)*1.1)
ax.set_ylim(minimum(MP)*0.9, maximum(MP)*1.11)
plt.legend(fontsize=20, bbox_to_anchor=[1.22,1.0])
ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)
ax1.spines["right"].set_visible(false)
ax1.spines["top"].set_visible(false)
plt.tight_layout()
plt.savefig("figures/figure2/time_trace.svg", format="svg")
gcf()
