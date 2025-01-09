### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 724f1c36-0183-4686-8025-f8fc1d610911
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ 63f71b7c-cebb-11ef-0b10-f5b7374ea6a5
using PlutoUI

# ╔═╡ 519fa022-4be2-49aa-b0c6-ddbbd5cacac9
using DifferentialEquations, ModelingToolkit, Unitful, Plots, Random, Distributions;

# ╔═╡ baf3b164-37a7-4af2-850f-113432e91637
TableOfContents()

# ╔═╡ f79f36eb-4068-4548-8b92-15e570df2bdf
md"# Libraries"

# ╔═╡ cddad652-d401-4c49-8902-a33230f6a0b8
md"# Simple PBPK Example"

# ╔═╡ 6f3b37eb-e97e-47d8-b626-e7c2954df3bd
# ╠═╡ show_logs = false
LocalResource(
	"figs/simple-pbpk.png",
	:width => 500
)

# ╔═╡ 5810984e-508d-4a85-8be0-0809e7f19458
function PBPK(; name)
    @independent_variables t [description = "time", unit = u"hr"]
    Dt = Differential(t)
    
    pars = @parameters begin
		# volumes (L)
		Vmu = 29, [description = "muscle volume", unit = u"L"]
		Var = 1.65, [description = "arterial blood volume", unit = u"L"]
		Vve = 3.9, [description = "venous blood volume", unit = u"L"]
		Vlu = 0.5, [description = "lung volume", unit = u"L"]
		Vli = 1.8, [description = "liver volume", unit = u"L"]

		# blood flows (L/hr)
		Qmu = 0.17 * 6.5 * 60, [description = "muscle flow", unit = u"L*hr^-1"]
		Qli = 0.245 * 60 * 6.5, [description = "liver flow", unit = u"L*hr^-1"]
		Qlu = 6.5 * 60, [description = "lung flow", unit = u"L*hr^-1"]

		# partition coefficients
		Kpmu = 1.0, [description = "muscle:plasma"]
		Kpli = 1.0, [description = "liver:plasma"]
		Kpre = 1.0, [description = "rest:plasma"]
		Kplu = 1.0, [description = "lung:plasma"]

		# other
		BP = 1, [description = "blood:plasma"]
		WEIGHT = 73, [description = "weight", unit = u"L"]
		Cl_hepatic = 10.0, [description = "hepatic clearance", unit = u"L*hr^-1"]
		fup = 0.5, [description = "unbound fraction"]
    end

    vars = @variables begin
		(MUSCLE(t) = 0.0), [description = "muscle amount", unit = u"mg"]
        (LUNG(t) = 0.0), [description = "lung amount", unit = u"mg"]
		(LIVER(t) = 0.0), [description = "liver amount", unit = u"mg"]
		(REST(t) = 0.0), [description = "rest of body amount", unit = u"mg"]
		(VEN(t) = 10.0), [description = "liver amount", unit = u"mg"]
		(ART(t) = 0.0), [description = "liver amount", unit = u"mg"]
		CP(t), [description = "plasma conc", unit = u"mg*L^-1"]
	end

	# transformed parameters
	Vre = WEIGHT - (Vli + Vlu + Vve + Var + Vmu)
	Qre = Qlu - (Qli + Qmu)
	
	# concentrations
	Carterial = ART / Var
	Cmuscle = MUSCLE / Vmu
	Cvenous = VEN / Vve
	Cliver = LIVER / Vli
	Crest = REST / Vre
	Clung = LUNG / Vlu
	
    observed = [
        CP ~ Cvenous / BP
    ]

	eqs = [
		Dt(MUSCLE) ~ Qmu * (Carterial - Cmuscle / (Kpmu / BP)),
		Dt(REST) ~ Qre * (Carterial - Crest / (Kpre / BP)),
		Dt(LUNG) ~ Qlu * (Cvenous - Clung / (Kplu / BP)),
		Dt(LIVER) ~ Qli * (Carterial - Cliver / (Kpli / BP)) - Cl_hepatic * fup * (Cliver / (Kpli/BP)),
		Dt(ART) ~ Qlu * (Clung / (Kplu / BP) - Carterial),
		Dt(VEN) ~ Qmu * (Cmuscle / (Kpmu/BP)) + Qli * (Cliver / (Kpli/BP)) + Qre * (Crest / (Kpre/BP)) - Qlu * Cvenous
	]

	ODESystem(eqs, t, vars, pars; name=name, observed=observed) 
end

# ╔═╡ e1fafda4-0dfa-45ef-bb66-00084217ab48
@mtkbuild pbpk = PBPK()

# ╔═╡ 1b070a9f-f584-43ef-a55b-08b821b675c3
# ╠═╡ show_logs = false
prob = ODEProblem(pbpk, [], (0.0, 24.0), [])

# ╔═╡ 6be467dc-d381-44fd-a547-8cd0d47e1ebd
sol = solve(prob, Tsit5())

# ╔═╡ 3d8413b9-680c-48fa-be96-38cb5822adef
plot(sol)

# ╔═╡ 01c88ade-449d-41f0-9532-0ca34a3e7b1b
plot(sol, idxs = :CP)

# ╔═╡ eb74aabb-e3fb-4e88-ac2f-df20c39dc932
plot(sol, idxs = :CP, yscale = :log10)

# ╔═╡ 6bfa8874-2cc1-43d4-802a-1889eceb43ad
idx_ven = ModelingToolkit.variable_index(pbpk, :VEN)

# ╔═╡ bb52d7de-235f-4d91-adcc-4647c5a773cf
begin
	affect!(integrator) = integrator.u[idx_ven] += 10.0
	cb = PresetTimeCallback([24.0, 48.0], affect!)
end

# ╔═╡ 1d0f13e0-1e43-4591-a0fb-3cb61558b1fc
begin
	prob_cb = remake(prob, tspan=(0.0,72.0))
	sol_cb = solve(prob_cb, callback=cb)
end

# ╔═╡ 1da695d6-7efd-4450-a4c1-fd6af5ee864b
plot(sol_cb)

# ╔═╡ 29043929-0261-466e-ac2e-0e949e48ce7f
begin
	## create simulation function
	Random.seed!(123)  # set seed for reproducibility
	function prob_func(prob, i, repeat)
	    remake(prob; p = [:Cl_hepatic => rand(LogNormal(log(10.0), 0.3))])
	end

	## create an ensemble problem and solve
	ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
	ensemble_sol = solve(ensemble_prob, Tsit5(), trajectories = 10)
end

# ╔═╡ c0f7b4dc-975a-44bc-b762-a4b1dfb9fa3c
plot(ensemble_sol, idxs=:MUSCLE)

# ╔═╡ Cell order:
# ╟─63f71b7c-cebb-11ef-0b10-f5b7374ea6a5
# ╟─baf3b164-37a7-4af2-850f-113432e91637
# ╟─f79f36eb-4068-4548-8b92-15e570df2bdf
# ╠═724f1c36-0183-4686-8025-f8fc1d610911
# ╠═519fa022-4be2-49aa-b0c6-ddbbd5cacac9
# ╟─cddad652-d401-4c49-8902-a33230f6a0b8
# ╟─6f3b37eb-e97e-47d8-b626-e7c2954df3bd
# ╠═5810984e-508d-4a85-8be0-0809e7f19458
# ╠═e1fafda4-0dfa-45ef-bb66-00084217ab48
# ╠═1b070a9f-f584-43ef-a55b-08b821b675c3
# ╠═6be467dc-d381-44fd-a547-8cd0d47e1ebd
# ╠═3d8413b9-680c-48fa-be96-38cb5822adef
# ╠═01c88ade-449d-41f0-9532-0ca34a3e7b1b
# ╠═eb74aabb-e3fb-4e88-ac2f-df20c39dc932
# ╠═6bfa8874-2cc1-43d4-802a-1889eceb43ad
# ╠═bb52d7de-235f-4d91-adcc-4647c5a773cf
# ╠═1d0f13e0-1e43-4591-a0fb-3cb61558b1fc
# ╠═1da695d6-7efd-4450-a4c1-fd6af5ee864b
# ╠═29043929-0261-466e-ac2e-0e949e48ce7f
# ╠═c0f7b4dc-975a-44bc-b762-a4b1dfb9fa3c
