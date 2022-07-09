using ACE1
using IPFitting

al = IPFitting.Data.read_xyz("./gain.xyz", energy_key="energy", force_key="forces", virial_key="virial")[1:5:end];

r0 = 3.8 

ACE_B = ACE1.Utils.rpi_basis(species = [:Ga, :In],
                              N = 3,
                              r0 = r0,
                              rin = 0.6 * r0,
                              rcut = 5.5,
                              maxdeg = 6)
                              
Bpair = pair_basis(species = [:Ga, :In],
      r0 = r0,
      maxdeg = 6,
      rcut = 7.0,
      pcut = 1,
      pin = 0)  

B = JuLIP.MLIPs.IPSuperBasis([Bpair, ACE_B]);

dB = LsqDB("", B, al)

Vref = OneBody(:Ga => -.0001, :In => -.0001)

weights = Dict(
        "default" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 ),
        "solid" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 ),
        "liquid" => Dict("E" => 30.0, "F" => 1.0 , "V" => 1.0 ))


solver_type = :lsqr
laplace_precon = false 

if solver_type == :lsqr
        solver = Dict(
                "solver" => :lsqr,
                "lsqr_damp" => 1e-2,
                "lsqr_atol" => 1e-6)
elseif solver_type == :rrqr
        solver = Dict(
                "solver" => :rrqr,
                "rrqr_tol" => 1e-5)
elseif solver_type == :brr
        solver = Dict(
                "solver" => :brr,
                "brr_tol" => 1e-3)
elseif solver_type == :ard 
        solver= Dict(
                "solver" => :ard,
                "ard_tol" => 1e-3,
                "ard_n_iter" => 10000,
                "ard_threshold_lambda" => 10000)
end

if laplace_precon
        using LinearAlgebra
        rlap_scal = 3.0
        P = Diagonal(vcat(ACE1.scaling.(dB.basis.BB, rlap_scal)...))
        solver["P"] = P
end

@show weights
IP, lsqinfo = lsqfit(dB, solver=solver, weights=weights, Vref=Vref, error_table = true);


rmse_table(lsqinfo["errors"])

save_dict("./fit_pot.json", Dict("IP" => write_dict(IP), "info" => lsqinfo))

ACE1.ExportMulti.export_ACE("./fit_pot.yace", IP)
