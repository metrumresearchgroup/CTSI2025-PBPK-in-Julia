### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8b19cc51-8ddf-4f4d-8e85-06953641d9f7
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ d242c44a-ce2e-11ef-3c3a-b366eaada965
using PlutoUI

# ╔═╡ 25751995-8656-4346-8dca-2df4c342b095
using DifferentialEquations, ModelingToolkit, Unitful, Plots, CSV, Tidier, ComponentArrays, Optimization, OptimizationOptimJL

# ╔═╡ 1f85fbdc-7a64-4def-8e8f-2aafec12bdbb
TableOfContents()

# ╔═╡ 98a5a8f5-79ef-4cff-ab2e-592b5358a90f
md"# PBPK example"

# ╔═╡ 3179d612-eb4e-4007-b3c5-18bbc2d8e4ed
md"
Let's build a PBPK model for voriconazole
"

# ╔═╡ 9d8005dc-89bd-42e4-a0b5-abcb75d3e80b
# ╠═╡ show_logs = false
LocalResource(
	"figs/pbpk.png",
	:width => 400
)

# ╔═╡ e308fc4a-157d-486b-8d06-317cbad4283d
md"## Libraries"

# ╔═╡ e8c04c18-9ed6-4e06-a7a9-362b6bb55de2
md"## Model development"

# ╔═╡ 57966d02-b00e-4e65-b46d-43d9d343123d
PBPK = function(; name)
	@independent_variables t, [description = "time"]
	Dt = Differential(t)

	pars = @parameters begin
		# volumes (L); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
		Vad = 18.2
  		Vbo = 10.5 
  		Vbr = 1.45
  		Vgu = 0.65
  		Vhe = 0.33
  		Vki = 0.31
  		Vli = 1.8
  		Vlu = 0.5
  		Vmu = 29
  		Vsp = 0.15
  		Vbl = 5.6

		# blood flows (L/h); Cardiac output = 6.5 (L/min); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981 
		Qad = 0.05*6.5*60
  		Qbo = 0.05*6.5*60
  		Qbr = 0.12*6.5*60
  		Qgu = 0.15*6.5*60 
  		Qhe = 0.04*6.5*60
  		Qki = 0.19*6.5*60
  		Qmu = 0.17*6.5*60
  		Qsp = 0.03*6.5*60
  		Qha = 0.065*6.5*60  
  		Qlu = 6.5*60 

		# partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
  		Kpad = 9.89   # adipose:plasma
  		Kpbo = 7.91   # bone:plasma
  		Kpbr = 7.35   # brain:plasma
  		Kpgu = 5.82   # gut:plasma
  		Kphe = 1.95   # heart:plasma
  		Kpki = 2.9    # kidney:plasma
  		Kpli = 4.66   # liver:plasma
  		Kplu = 0.83   # lungs:plasma
  		Kpmu = 2.94   # muscle:plasma; optimized
  		Kpsp = 2.96   # spleen:plasma
  		Kpre = 4      # calculated as average of non adipose Kps
  		BP = 1.0      # blood:plasma ratio

		# other parameters
  		WEIGHT = 73
  		ka = 0.849  #2.0  #0.849   # absorption rate constant (/hr) 
  		fup = 0.42   # fraction of unbound drug in plasma
  
  		# in vitro hepatic clearance parameters http://dmd.aspetjournals.org/content/38/1/25.long
  		fumic = 0.711 # fraction of unbound drug in microsomes
  		MPPGL = 30.3  # adult mg microsomal protein per g liver (mg/g)
  		VmaxH = 40    # adult hepatic Vmax (pmol/min/mg)
  		KmH = 9.3     # adult hepatic Km (uM)
  
  		# renal clearance  https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
  		CL_Ki = 0.096 # (L/hr) renal clearance
	end

	vars = @variables begin
		GUTLUMEN(t) = 200.0 
		GUT(t) = 0.0 
		ADIPOSE(t) = 0.0
		BRAIN(t) = 0.0
		HEART(t) = 0.0
		BONE(t) = 0.0
  		KIDNEY(t) = 0.0
		LIVER(t) = 0.0
		LUNG(t) = 0.0
		MUSCLE(t) = 0.0
		SPLEEN(t) = 0.0
		REST(t) = 0.0
  		ART(t) = 0.0
		VEN(t) = 0.0
		CP(t)
	end

	# additional volume derivations
  	Vve = 0.705*Vbl  # venous blood
    Var = 0.295*Vbl  # arterial blood
    Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+Vgu+Vbl)  # volume of rest of the body compartment
  
    # additional blood flow derivation
    Qli = Qgu + Qsp + Qha
    Qre = Qlu - (Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr)
  
    # intrinsic hepatic clearance calculation
  	CL_Li = ((VmaxH/KmH)*MPPGL*Vli*1000*60*1e-6) / fumic  # (L/hr) hepatic clearance

	# Calculation of tissue drug concentrations (mg/L)
  	Cadipose = ADIPOSE/Vad
    Cbone = BONE/Vbo
    Cbrain = BRAIN/Vbr 
    Cheart = HEART/Vhe 
    Ckidney = KIDNEY/Vki
    Cliver = LIVER/Vli 
    Clung = LUNG/Vlu 
    Cmuscle = MUSCLE/Vmu
    Cspleen = SPLEEN/Vsp
    Crest = REST/Vre
    Carterial = ART/Var
    Cvenous = VEN/Vve
    Cgut = GUT/Vgu

	observed = [
		CP ~ Cvenous/BP
	]

	eqs = [
  		Dt(GUTLUMEN) ~ -ka*GUTLUMEN,
  		Dt(GUT) ~ ka*GUTLUMEN + Qgu*(Carterial - Cgut/(Kpgu/BP)),
  		Dt(ADIPOSE) ~ Qad*(Carterial - Cadipose/(Kpad/BP)),
  		Dt(BRAIN) ~ Qbr*(Carterial - Cbrain/(Kpbr/BP)),
  		Dt(HEART) ~ Qhe*(Carterial - Cheart/(Kphe/BP)),
  		Dt(KIDNEY) ~ Qki*(Carterial - Ckidney/(Kpki/BP)) - CL_Ki*(fup*Ckidney/(Kpki/BP)),
  		Dt(LIVER) ~ Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - CL_Li*(fup*Cliver/(Kpli/BP)),
  		Dt(LUNG) ~ Qlu*(Cvenous - Clung/(Kplu/BP)),
  		Dt(MUSCLE) ~ Qmu*(Carterial - Cmuscle/(Kpmu/BP)),
  		Dt(SPLEEN) ~ Qsp*(Carterial - Cspleen/(Kpsp/BP)),
  		Dt(BONE) ~ Qbo*(Carterial - Cbone/(Kpbo/BP)),
  		Dt(REST) ~ Qre*(Carterial - Crest/(Kpre/BP)),
  		Dt(VEN) ~ Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
    Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
    Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qre*(Crest/(Kpre/BP)) - Qlu*Cvenous,
  		Dt(ART) ~ Qlu*(Clung/(Kplu/BP) - Carterial)
	]

	ODESystem(eqs, t, vars, pars; name=name, observed=observed)
end

# ╔═╡ 106b1b40-4458-453d-b5bd-931969c9471b
@mtkbuild pbpk = PBPK()

# ╔═╡ d63b718e-1847-4c92-90c7-6b507293582f
md"## Simulation"

# ╔═╡ 7c4f551f-2667-4e6b-955a-568a993261a8
md"### Create ODE problem"

# ╔═╡ 93e5fd52-20aa-48dd-8cbc-d8c5f934f70b
prob = ODEProblem(pbpk, [], (0.0, 24.0), [])

# ╔═╡ e9d019eb-1244-47ca-b1d7-74dae4db3975
md"### Solve"

# ╔═╡ 858cda2e-7fde-4120-a60f-d104e429b5e8
sol = solve(prob, Tsit5())

# ╔═╡ 8c5cbd36-ba16-44c4-a9bd-27927301765b
plot(sol, idxs = :CP)

# ╔═╡ 8039ca65-bb01-4428-b3c5-6de0c866912c
md"## Optimization"

# ╔═╡ aa86aafd-e409-4180-a51b-15b74aae8ffd
md"### Read observed data"

# ╔═╡ 8aa63f09-d953-4459-a50a-3a2665722b7c
# https://link.springer.com/article/10.1007/s00228-010-0869-3
dat = CSV.read("shi2010-fig1.csv", DataFrame)

# ╔═╡ b4891fca-ccf1-40f8-b437-f9b8fab4c1f5
md"### Initial simulation"

# ╔═╡ 943810c3-7f3a-4502-9970-199f1e40578c
sol_init = solve(prob, Tsit5(), saveat = dat.time)

# ╔═╡ a675ead7-4d59-406d-a37f-66dda9d0ee13
begin
	scatter(dat.time, dat.conc)
	plot!(sol_init, idxs = :CP)
end

# ╔═╡ 6b4a2a16-6856-417a-9a91-e00f71c62b91
md"### Create loss function"

# ╔═╡ 651760f9-9faa-47a4-a5b9-8b3860390c36
function loss(x, p, pred=false)
    # extract info
    odeprob, dat = p  # ODEProblem stored as parameters to avoid using global variables
    # remake the problem, passing in our new parameter object
    newprob = remake(odeprob; p = Dict([:ka => x.ka, :Kpmu => x.Kpmu, :Kpli => x.Kpli, :Kpad => x.Kpad, :BP => x.BP]))
    sol = solve(newprob, Tsit5(), saveat = dat.time);
    if pred; return sol; end
	return sum(abs2, dat.conc .- sol[:CP])
end

# ╔═╡ 8f473f67-ff48-44cb-a948-84f1c2315588
# define initial param values
θ = ComponentArray(ka = prob.ps[:ka], Kpmu = prob.ps[:Kpmu], Kpli = prob.ps[:Kpli], Kpad = prob.ps[:Kpad], BP = prob.ps[:BP])

# ╔═╡ 8831e5b7-5359-4f29-9095-36808fa19c35
# define constant parameters
p_const = (prob, dat);

# ╔═╡ b23f7e84-eb0c-40a6-893a-dc040e13d357
md"### Optimize"

# ╔═╡ e28f9182-a414-44ac-9c76-6c855bae1f6c
begin
	# define optimization function and problem then optimize 
	optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
	optprob = OptimizationProblem(optfn, θ, p_const, lb = [0.1,0.01,0.1,0.1,0.5], ub = [5.0,20.0,20.0,30.0,2.0])
	p_optim = solve(optprob, Optim.NelderMead())
end

# ╔═╡ ebd2719a-f5fc-4ba0-849d-cee4421e8b37
p_optim.retcode

# ╔═╡ d0ffcaa0-0569-4457-8eb0-aa3e829b42b9
println(θ)

# ╔═╡ c26dc0d1-e27d-4b21-8881-f87e99b79e08
println(p_optim.u)

# ╔═╡ 7e417e58-455f-4a22-b196-a9b82c1484ee
md"### Validate"

# ╔═╡ 7764c9e9-b66b-4679-8512-40d7de180c8a
pred_before = loss(θ, p_const, true)

# ╔═╡ 584f8a03-231f-47a3-bab6-1391e412b448
pred_after = loss(p_optim.u, p_const, true)

# ╔═╡ ec600ce9-68c3-4c84-98bc-447a943a7ac2
begin
	scatter(dat.time, dat.conc, label = "data")
	plot!(pred_before, idxs = :CP, label = "initial", linestyle = :dash)
	plot!(pred_after, idxs = :CP, label = "optimized")
end

# ╔═╡ Cell order:
# ╟─d242c44a-ce2e-11ef-3c3a-b366eaada965
# ╟─1f85fbdc-7a64-4def-8e8f-2aafec12bdbb
# ╟─98a5a8f5-79ef-4cff-ab2e-592b5358a90f
# ╟─3179d612-eb4e-4007-b3c5-18bbc2d8e4ed
# ╟─9d8005dc-89bd-42e4-a0b5-abcb75d3e80b
# ╟─e308fc4a-157d-486b-8d06-317cbad4283d
# ╠═8b19cc51-8ddf-4f4d-8e85-06953641d9f7
# ╠═25751995-8656-4346-8dca-2df4c342b095
# ╟─e8c04c18-9ed6-4e06-a7a9-362b6bb55de2
# ╠═57966d02-b00e-4e65-b46d-43d9d343123d
# ╠═106b1b40-4458-453d-b5bd-931969c9471b
# ╟─d63b718e-1847-4c92-90c7-6b507293582f
# ╟─7c4f551f-2667-4e6b-955a-568a993261a8
# ╠═93e5fd52-20aa-48dd-8cbc-d8c5f934f70b
# ╟─e9d019eb-1244-47ca-b1d7-74dae4db3975
# ╠═858cda2e-7fde-4120-a60f-d104e429b5e8
# ╠═8c5cbd36-ba16-44c4-a9bd-27927301765b
# ╟─8039ca65-bb01-4428-b3c5-6de0c866912c
# ╟─aa86aafd-e409-4180-a51b-15b74aae8ffd
# ╠═8aa63f09-d953-4459-a50a-3a2665722b7c
# ╟─b4891fca-ccf1-40f8-b437-f9b8fab4c1f5
# ╠═943810c3-7f3a-4502-9970-199f1e40578c
# ╠═a675ead7-4d59-406d-a37f-66dda9d0ee13
# ╟─6b4a2a16-6856-417a-9a91-e00f71c62b91
# ╠═651760f9-9faa-47a4-a5b9-8b3860390c36
# ╠═8f473f67-ff48-44cb-a948-84f1c2315588
# ╠═8831e5b7-5359-4f29-9095-36808fa19c35
# ╟─b23f7e84-eb0c-40a6-893a-dc040e13d357
# ╠═e28f9182-a414-44ac-9c76-6c855bae1f6c
# ╠═ebd2719a-f5fc-4ba0-849d-cee4421e8b37
# ╠═d0ffcaa0-0569-4457-8eb0-aa3e829b42b9
# ╠═c26dc0d1-e27d-4b21-8881-f87e99b79e08
# ╟─7e417e58-455f-4a22-b196-a9b82c1484ee
# ╠═7764c9e9-b66b-4679-8512-40d7de180c8a
# ╠═584f8a03-231f-47a3-bab6-1391e412b448
# ╠═ec600ce9-68c3-4c84-98bc-447a943a7ac2
