### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ c7030848-59db-42d9-b827-733f2afccd0f
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ 76902963-43fc-405d-9132-2ed8409fa56d
using PlutoUI

# ╔═╡ 01d1f036-99f2-4ef6-8ab3-f1e5180fd5c5
using DifferentialEquations, ModelingToolkit, Unitful, Plots, Random, Distributions;

# ╔═╡ 179bbe7a-16f1-466f-baa0-a132a4f916ce
TableOfContents()

# ╔═╡ c3f2054a-8df4-41ea-9062-e08fb4b1792b
md"
# PK example
"

# ╔═╡ a0e14355-0fbd-4b98-8790-e7ec50910a61
md"
Let's build a one-compartment PK model with oral absorption
"

# ╔═╡ 83701d28-60a8-4c11-8923-ef7be27d384d
# ╠═╡ show_logs = false
LocalResource(
	"figs/pk.png",
	:width => 400
)

# ╔═╡ 3a1c3888-1671-4bcc-bd0f-884d6e0e9365
md"
## Libraries
"

# ╔═╡ c4a4b27e-d053-4f70-aa2a-5ab7538780ee
md"
## Model development
"

# ╔═╡ ac7b62aa-29ec-4b61-a01e-8a582b10e928
function PK(; name)
    
    @independent_variables t [description = "time", unit = u"hr"]
    Dt = Differential(t)
    
    pars = @parameters begin
        CL = 4.0, [description = "clearance", unit = u"L*hr^-1"]
        V = 35.0, [description = "volume", unit = u"L"]
        ka = 2.0, [description = "absorption rate constant", unit = u"hr^-1"]
    end

    vars = @variables begin
		(depot(t) = 320.0), [description = "Theophylline amount in depot compartment", unit = u"mg"]
        (cent(t) = 0.0), [description = "Theophylline amount in central compartment", unit = u"mg"]
		c_cent(t), [description = "Theophylline concentration in central compartment", unit = u"mg*L^-1"]
	end
	
    observed = [
        c_cent ~ cent/V
    ]

	eqs = [
		Dt(depot) ~ -ka * depot,
		Dt(cent) ~ ka * depot - (CL / V) * cent
	]

	ODESystem(eqs, t, vars, pars; name=name, observed=observed, checks = ~ModelingToolkit.CheckUnits)  # 100124 unit checking is broken at this point so will suppress for now
end

# ╔═╡ 3604e6c4-a29d-45ee-af90-ae6d480a185b
@mtkbuild pk = PK()

# ╔═╡ dd493723-8a28-41ca-9693-7ad328a4a618
ModelingToolkit.UnitfulUnitCheck.get_unit(pk.CL)

# ╔═╡ c1505fa8-b7fa-4f24-9a37-575cb87dcd89
ModelingToolkit.getdescription(pk.cent) 

# ╔═╡ b4c202a6-5561-46ab-98a8-ed157c6d9742
ModelingToolkit.unknowns(pk)

# ╔═╡ 5cf34c06-cc4a-4baf-8c9b-1017d006c971
ModelingToolkit.equations(pk) 

# ╔═╡ bbbe3a0c-7553-46c4-9543-ce8afba41cb5
ModelingToolkit.observed(pk)

# ╔═╡ 7ffd9484-d960-4d7a-b9ae-ca09d0ba5e62
md"
## Simulation
"

# ╔═╡ df3ceb01-8b50-4771-89cb-7c54b77f5a03
md"
### Create ODE problem
"

# ╔═╡ 1a062802-ee3b-493c-b260-9a8c784a1f21
prob = ODEProblem(pk, [], (0.0, 25.0), [])

# ╔═╡ c49719ee-7f27-4eaf-aeb4-ee551e91427c
md"
### Solve
"

# ╔═╡ 445b8741-6716-4684-8c07-2e171f070abb
sol = solve(prob, Tsit5())

# ╔═╡ d2582ca6-5ff3-49d4-8241-9703f976dd85
md"### Plot"

# ╔═╡ 73f7a96d-5f66-4ab4-82bb-0627a6139f8a
plot(sol)

# ╔═╡ 841eaf44-a093-4cc9-a792-184d11622dac
plot(sol, idxs = :c_cent)

# ╔═╡ 1ae0468a-0684-405d-addb-a918030df470
md"### Callbacks"

# ╔═╡ 12777c22-e6b2-47fb-aff6-3a533faa1680
idx_depot = ModelingToolkit.variable_index(pk, :depot)

# ╔═╡ 27799b5f-c1a4-4539-ae05-5fa3f43ce122
begin
	affect!(integrator) = integrator.u[idx_depot] += 320.0
	cb = PresetTimeCallback([24.0, 48.0], affect!)
end

# ╔═╡ 7f8f500f-3784-4b10-89a3-6212a18ead24
begin
	prob_cb = remake(prob, tspan=(0.0,72.0))
	sol_cb = solve(prob_cb, callback=cb)
end


# ╔═╡ fcc3498d-aae0-4ffd-8b09-b1ba19935be9
plot(sol_cb)

# ╔═╡ c26b08ef-bbcf-47c8-8a7f-3e2a23bd9c06
md"### Update parameter"

# ╔═╡ b2f61109-3de3-4b8e-9ab7-457d32afd598
prob.ps[:ka]

# ╔═╡ 62da16d4-f06c-4f89-88b3-5fc80e07f0fc
prob_up = remake(prob; p = [:ka => 0.5])

# ╔═╡ 96df73e4-19f7-4a4a-b70e-aac0554a6252
prob_up.ps[:ka]

# ╔═╡ 77cb307b-4286-4e5d-b4a4-c1f421922682
begin
	sol_up = solve(prob_up, Tsit5())
	plot(sol_up, idxs = :c_cent)
end

# ╔═╡ 3b614ec7-564b-4bbe-9e02-442a4c93cff2
md"### Population simulation"

# ╔═╡ 360ee3c6-2c2e-48ff-b91e-44d96ae65fbf
begin
	## create simulation function
	Random.seed!(123)  # set seed for reproducibility
	function prob_func(prob, i, repeat)
	    remake(prob; p = [:CL => rand(LogNormal(log(4.0), 0.2))])
	end

	## create an ensemble problem and solve
	ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
	ensemble_sol = solve(ensemble_prob, Tsit5(), trajectories = 10)
end

# ╔═╡ a3496c93-c819-4074-b3ce-553677ac2520
# plotting a specific output
plot(ensemble_sol, xlab="Time(h)", ylab="Concentration (mg/L)", idxs = [:c_cent])

# ╔═╡ 0766abe9-13db-4466-9f23-32a009995e3b
md"#### Parallelization"

# ╔═╡ 0dfe27d1-77d0-429f-8d5d-107a06a6b8b7
@time ensemble_sol_serial = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = 10);

# ╔═╡ 6f0c09a8-74b4-4e0d-9a52-7144c9ec4d5f

@time ensemble_sol_parallel = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 10);


# ╔═╡ Cell order:
# ╟─76902963-43fc-405d-9132-2ed8409fa56d
# ╟─179bbe7a-16f1-466f-baa0-a132a4f916ce
# ╟─c3f2054a-8df4-41ea-9062-e08fb4b1792b
# ╟─a0e14355-0fbd-4b98-8790-e7ec50910a61
# ╟─83701d28-60a8-4c11-8923-ef7be27d384d
# ╟─3a1c3888-1671-4bcc-bd0f-884d6e0e9365
# ╠═c7030848-59db-42d9-b827-733f2afccd0f
# ╠═01d1f036-99f2-4ef6-8ab3-f1e5180fd5c5
# ╟─c4a4b27e-d053-4f70-aa2a-5ab7538780ee
# ╠═ac7b62aa-29ec-4b61-a01e-8a582b10e928
# ╠═3604e6c4-a29d-45ee-af90-ae6d480a185b
# ╠═dd493723-8a28-41ca-9693-7ad328a4a618
# ╠═c1505fa8-b7fa-4f24-9a37-575cb87dcd89
# ╠═b4c202a6-5561-46ab-98a8-ed157c6d9742
# ╠═5cf34c06-cc4a-4baf-8c9b-1017d006c971
# ╠═bbbe3a0c-7553-46c4-9543-ce8afba41cb5
# ╟─7ffd9484-d960-4d7a-b9ae-ca09d0ba5e62
# ╟─df3ceb01-8b50-4771-89cb-7c54b77f5a03
# ╠═1a062802-ee3b-493c-b260-9a8c784a1f21
# ╟─c49719ee-7f27-4eaf-aeb4-ee551e91427c
# ╠═445b8741-6716-4684-8c07-2e171f070abb
# ╟─d2582ca6-5ff3-49d4-8241-9703f976dd85
# ╠═73f7a96d-5f66-4ab4-82bb-0627a6139f8a
# ╠═841eaf44-a093-4cc9-a792-184d11622dac
# ╟─1ae0468a-0684-405d-addb-a918030df470
# ╠═12777c22-e6b2-47fb-aff6-3a533faa1680
# ╠═27799b5f-c1a4-4539-ae05-5fa3f43ce122
# ╠═7f8f500f-3784-4b10-89a3-6212a18ead24
# ╠═fcc3498d-aae0-4ffd-8b09-b1ba19935be9
# ╟─c26b08ef-bbcf-47c8-8a7f-3e2a23bd9c06
# ╠═b2f61109-3de3-4b8e-9ab7-457d32afd598
# ╠═62da16d4-f06c-4f89-88b3-5fc80e07f0fc
# ╠═96df73e4-19f7-4a4a-b70e-aac0554a6252
# ╠═77cb307b-4286-4e5d-b4a4-c1f421922682
# ╟─3b614ec7-564b-4bbe-9e02-442a4c93cff2
# ╠═360ee3c6-2c2e-48ff-b91e-44d96ae65fbf
# ╠═a3496c93-c819-4074-b3ce-553677ac2520
# ╟─0766abe9-13db-4466-9f23-32a009995e3b
# ╠═0dfe27d1-77d0-429f-8d5d-107a06a6b8b7
# ╠═6f0c09a8-74b4-4e0d-9a52-7144c9ec4d5f
