### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 6f73203a-88b7-43cd-97bf-845384540559
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ 0bf4c192-cee6-11ef-1f91-1fd13e6da49e
using PlutoUI

# ╔═╡ 5c029511-b19e-4542-821b-90e55c892d1e
using DifferentialEquations, ModelingToolkit, Plots, CSV, Tidier, ComponentArrays, Optimization, OptimizationOptimJL, OptimizationNLopt

# ╔═╡ 85db720a-2909-4073-957a-bb118babed7d
TableOfContents()

# ╔═╡ 9a18f8ae-c74c-48b2-b1e4-a3adf377113f
md"# mAb mPBPK example"

# ╔═╡ 8a712e88-d36b-4587-90ab-04fcdd2933ab
md"## Libraries"

# ╔═╡ e99f2f6c-ea60-41be-9021-64f7417e29e6
md"## Model development"

# ╔═╡ 434bca6d-cd0f-4a07-91bf-668a20df2699
# Reproduction of mouse model from Cao et al., 2014 “Incorporating Target-Mediated Drug Disposition in a Minimal Physiologically-Based Pharmacokinetic Model for Monoclonal Antibodies.” Journal of Pharmacokinetics and Pharmacodynamics 41 (4): 375–87.

#The model implemented here represents Model B in the publication, with the modification that the steady-state assumption was not used, and everything is expressed in concentrations instead.

# This version includes infusion dosing into the plasma compartment

function PBPK(; name)
	@independent_variables t 
	Dt = Differential(t)

	pars = @parameters begin
		s_per_h = 3600.0
		s_per_d = 3600.0*24.0
		KD_FAPa = 3.73e-11*(1e9) # DB:FAPalpha dissociation constant in cyno (M to nM)
	    k_on = 7.7e5*(1e-9)      # DB:FAPalpha binding rate constant in cyno
	    k_int = log(2)/(24.0 * s_per_h)
	    Thalf_FAPa = 24.0
	    k_syn = 0.3/s_per_h
		ISF = 0.579        # cyno total interstitial volume (L)
	    V_plasma = 0.0966  # cyno plasma volume (L)
	    V_lymph = 0.193    # cyno lymph volume (L)
	    L = 0.275*(1/s_per_d) #0.2*(1/s_per_d),  # cyno total lymph flow (L/d to L/s)
	    sig_tight = 0.945    # Sigma_1, vascular reflection coefficient for tight tissue (Table 1) 
	    sig_leaky = 0.687   # Sigma_2, vascular reflection coefficient for leaky tissue (Table 1) 
	    sig_L = 0.2        # lymphatic capillary reflection coefficient (unitless), predefined as in 'several previous PBPK models'
	    Kp = 0.8           # available fraction of ISF for mAb binding (unitless)
	    CL_p = 0.35*2.6*(10^-3)*(1/s_per_h)/2 # clearance from plasma (mL/h/kg to L/s)
	    dose_mgkg = 1.0   # dose in mg/kg
		cyno_WT = 2.6   # cyno weight = 3 kg
		nmol_per_mol = 1.0e9
		DB_MW = 150000.0    # molecular weight
		dose = dose_mgkg*(cyno_WT)*(1/DB_MW)*(1/1e3)*(nmol_per_mol)
		infusion_dur = 30*60 # duration of infusion (30 min to s)
	    infusion_rate(t) = dose/infusion_dur   # rate of infusion (nmol/s)
	    k_off = KD_FAPa*k_on
	    k_deg = log(2)/(Thalf_FAPa * s_per_h)
	    V_tight = 0.30108 #0.65 * ISF * Kp,
	    V_leaky = 0.16212 #0.35 * ISF * Kp,
		L_tight = 0.33 * L
    	L_leaky = 0.67 * L
	end

	vars = @variables begin
		C_plasma(t) = 0.0
		R_tight(t) = k_syn/k_deg
		C_tight(t) = 0.0
		AR_tight(t) = 0.0
		R_leaky(t) = k_syn/k_deg
		C_leaky(t) = 0.0
		AR_leaky(t) = 0.0
		C_lymph(t) = 0.0
	end
	
	eqs = [
		# Antibody concentration in plasma - Cp in 2014 paper, Eq (1)
	    Dt(C_plasma) ~ 1/V_plasma * (-L_leaky * (1 - sig_leaky) * C_plasma - L_tight * (1 - sig_tight) * C_plasma + L * C_lymph - CL_p * C_plasma + infusion_rate)
	    
	    # Receptor concentration in tight tissue compartment 
	    Dt(R_tight) ~ k_syn - R_tight*k_deg - k_on*C_tight*R_tight + k_off*AR_tight
	  
	    # Antibody concentration in tight tissue compartment
	    Dt(C_tight) ~ 1/V_tight * (L_tight * (1 - sig_tight) * C_plasma - L_tight * (1 - sig_L) * C_tight) - k_on*C_tight*R_tight + k_off*AR_tight
	
	    # Antibody:receptor complex concentration in leaky tissue compartment
	    Dt(AR_tight) ~ k_on*C_tight*R_tight - k_off*AR_tight - AR_tight*k_int
	
	    # Receptor concentration in leaky tissue compartment
	    Dt(R_leaky) ~ k_syn - R_leaky*k_deg - k_on*C_leaky*R_leaky + k_off*AR_leaky
	
	    # Antibody concentration in leaky tissue compartment
	    Dt(C_leaky) ~ 1/V_leaky * (L_leaky * (1 - sig_leaky) * C_plasma - L_leaky * (1 - sig_L) * C_leaky) - k_on*C_leaky*R_leaky + k_off*AR_leaky
	    
	    # Antibody:receptor complex concentration in leaky tissue compartment
	    Dt(AR_leaky) ~ k_on*C_leaky*R_leaky - k_off*AR_leaky - AR_leaky*k_int
	
	    # Antibody concentration in lymph
	    Dt(C_lymph) ~ 1/V_lymph * (L_tight * (1 - sig_L) * C_tight + L_leaky * (1 - sig_L) * C_leaky - L * C_lymph)
	]

	ODESystem(eqs, t, vars, pars; name=name, discrete_events=[[30*60] => [infusion_rate ~ 0.0]])
end

# ╔═╡ 6228f0fb-7385-4de4-b681-1b03cf621644
@mtkbuild pbpk = PBPK()

# ╔═╡ 4ce983fe-c11d-448c-ad72-61b7c8ec84d0
md"## Simulation"

# ╔═╡ fdec998d-7ab7-4b3a-86b2-52672a1d3da0
md"### Create ODE problem"

# ╔═╡ 33cb45b4-4c89-4f78-9284-0a4373a52673
prob = ODEProblem(pbpk, [], (0.0, 504.5*3600.0), [])

# ╔═╡ de34585d-6675-4901-82bd-a04c7e3b91d3
md"### Solve"

# ╔═╡ 4bbf1703-c5ea-4f48-8f57-99735862e2cc
sol = solve(prob, saveat = 60.0)

# ╔═╡ 65189443-8cd8-49cb-b458-226a03fa238b
md"### Plot"

# ╔═╡ 591357ca-8044-4656-8705-96d7f075f41b
plot(sol, idxs = pbpk.C_plasma)

# ╔═╡ d45886de-7193-4aaa-aa06-db5e088c3bdc
md"## Optimization"

# ╔═╡ aa333e4c-566d-42b4-a268-d4febe2c9768
md"### Read observed data"

# ╔═╡ dc053f7f-ada9-4c1c-838c-f3b12d736ac1
dat = CSV.read("7E3-plasma-conc.csv", DataFrame)

# ╔═╡ f43f6694-7fc0-4f55-80f1-3f1cb6d42055
md"### Initial simulation"

# ╔═╡ f8eb1cea-7d63-4bbd-9cf9-59df6b6aa6c0
begin
	timepoints = dat.time .* 3600.0
	sol2 = solve(prob, Tsit5(), saveat=timepoints)
	C_plasma_ugml = sol2[pbpk.C_plasma] * prob.ps[pbpk.DB_MW] * 1e-6
end

# ╔═╡ d61ed4a8-704b-4956-affb-fff5f9d62803
prob.ps[pbpk.DB_MW]

# ╔═╡ a5468278-f78f-439a-ab49-ae3a0957a0e3
begin
	scatter(dat.time, dat.conc)
	plot!(sol2.t/3600.0, sol2[pbpk.C_plasma])
end

# ╔═╡ f5408321-7ddb-4acd-8dcd-086f0382af03
sol2.t[4:end] ./ 3600

# ╔═╡ 577574f3-99ab-4ecc-bc6f-eddec689cb36
dat.time

# ╔═╡ ccb9efea-a4aa-46a0-aa40-f76c2f63f3db
prob.ps[:sig_tight]

# ╔═╡ 353bb679-70ea-41ad-b0ac-f30ab23a4111
md"### Create objective function"

# ╔═╡ c13c3283-6f2d-4efe-ad1a-2fdc3f74618c
function loss(x, p, pred=false)
    # extract info
    odeprob, dat = p  # ODEProblem stored as parameters to avoid using global variables
    # remake the problem, passing in our new parameter object
	newprob = remake(odeprob; u0 = Dict(), p = [:CL_p => x[1], :sig_tight => x[2], :sig_leaky => x[3], :k_syn => x[4]])
    #newprob = remake(odeprob; u0 = Dict(), p = Dict([:CL_p => x.CL_p, :sig_tight => x.sig_tight, :sig_leaky => x.sig_leaky, :k_syn => x.k_syn]))
    sol = solve(newprob, Tsit5(), saveat = dat.time .* 3600.0);
    if pred; return sol; end
    return sum(abs2, dat.conc .- sol[pbpk.C_plasma][4:end])
end

# ╔═╡ e2c11847-6374-4e5b-837b-6d17cd58a858
pp = remake(prob; p=[pbpk.CL_p => 2.0])

# ╔═╡ bbf75564-36f2-4694-9cb5-e40011dcf070
pp.ps[pbpk.CL_p]

# ╔═╡ 668ea1d4-ac4f-4982-9b0d-ed69288707b7
md"### Initial parameters"

# ╔═╡ 3188b60f-07cb-467e-af80-2c30cd2e9858
# define initial parameters
θ = ComponentArray(CL_p = prob.ps[:CL_p], sig_tight = prob.ps[:sig_tight], sig_leaky = prob.ps[:sig_leaky], k_syn = prob.ps[:k_syn])

# ╔═╡ f98f1780-c890-4034-ab36-66b9e7095fcc
# define constant parameters
p_const = (prob, dat);

# ╔═╡ dd8c0dc2-7a2a-45b7-8c3b-7dd7212ac555
begin
	# define callback to track losses
	losses = Float64[]
	
	callback = function (p, l)
	  push!(losses, l)
	  if length(losses)%50==0
	      println("Current loss after $(length(losses)) iterations: $(losses[end])")
	  end
	  return false
	end
end

# ╔═╡ 02b14f77-ddad-4519-aa2d-933e243896ca
# define optimization function and problem then optimize
begin
	## derivate-free 
	lb = [0.0, 0.0, 0.0, 0.0]
    ub = [1e-3, 0.999999, 0.999999, 0.3/3600.0]
	optfn = OptimizationFunction(loss)
	optprob = OptimizationProblem(optfn, θ, p_const, lb=lb, ub=ub)
	p_optim_df = solve(optprob, NLopt.LN_SBPLX())
end

# ╔═╡ ba1540b3-be85-4288-ae0a-d37775062f0c
0.3/3600

# ╔═╡ Cell order:
# ╟─0bf4c192-cee6-11ef-1f91-1fd13e6da49e
# ╟─85db720a-2909-4073-957a-bb118babed7d
# ╟─9a18f8ae-c74c-48b2-b1e4-a3adf377113f
# ╟─8a712e88-d36b-4587-90ab-04fcdd2933ab
# ╠═6f73203a-88b7-43cd-97bf-845384540559
# ╠═5c029511-b19e-4542-821b-90e55c892d1e
# ╟─e99f2f6c-ea60-41be-9021-64f7417e29e6
# ╠═434bca6d-cd0f-4a07-91bf-668a20df2699
# ╠═6228f0fb-7385-4de4-b681-1b03cf621644
# ╟─4ce983fe-c11d-448c-ad72-61b7c8ec84d0
# ╟─fdec998d-7ab7-4b3a-86b2-52672a1d3da0
# ╠═33cb45b4-4c89-4f78-9284-0a4373a52673
# ╟─de34585d-6675-4901-82bd-a04c7e3b91d3
# ╠═4bbf1703-c5ea-4f48-8f57-99735862e2cc
# ╟─65189443-8cd8-49cb-b458-226a03fa238b
# ╠═591357ca-8044-4656-8705-96d7f075f41b
# ╟─d45886de-7193-4aaa-aa06-db5e088c3bdc
# ╟─aa333e4c-566d-42b4-a268-d4febe2c9768
# ╠═dc053f7f-ada9-4c1c-838c-f3b12d736ac1
# ╟─f43f6694-7fc0-4f55-80f1-3f1cb6d42055
# ╠═f8eb1cea-7d63-4bbd-9cf9-59df6b6aa6c0
# ╠═d61ed4a8-704b-4956-affb-fff5f9d62803
# ╠═a5468278-f78f-439a-ab49-ae3a0957a0e3
# ╠═f5408321-7ddb-4acd-8dcd-086f0382af03
# ╠═577574f3-99ab-4ecc-bc6f-eddec689cb36
# ╠═ccb9efea-a4aa-46a0-aa40-f76c2f63f3db
# ╟─353bb679-70ea-41ad-b0ac-f30ab23a4111
# ╠═c13c3283-6f2d-4efe-ad1a-2fdc3f74618c
# ╠═e2c11847-6374-4e5b-837b-6d17cd58a858
# ╠═bbf75564-36f2-4694-9cb5-e40011dcf070
# ╟─668ea1d4-ac4f-4982-9b0d-ed69288707b7
# ╠═3188b60f-07cb-467e-af80-2c30cd2e9858
# ╠═f98f1780-c890-4034-ab36-66b9e7095fcc
# ╠═dd8c0dc2-7a2a-45b7-8c3b-7dd7212ac555
# ╠═02b14f77-ddad-4519-aa2d-933e243896ca
# ╠═ba1540b3-be85-4288-ae0a-d37775062f0c
