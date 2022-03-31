module quasiExact 
	import ..sorptivity, ..wrc, ..kunsat
	import BlackBoxOptim, Optim
 	export CONVERT_3D_2_1D, HYDRO_2_INFILTRATION3D, OF_QUASIEXACT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_2_TIMEη
	#		= COMPUTE NORMALISED TIME: Tη =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_2_TIMEη(Sorptivity, T, ΔK)
			return Time_η = T * 2.0 * (ΔK / Sorptivity) ^ 2.0
		end # function: TIME_2_TIMEη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_1D
	# 		= TRANSFORMS INFILTRATION_3D TO INFILTRATION_1D =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, option, T; θini= infiltParam.θini[iZ])

		Δθ = hydroInfilt.θs[iZ] - θini

		Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt, option, option.infilt)

		for iT = 1:N_Infilt[iZ]
			∑Infilt_1D[iZ,iT] = ∑Infilt_3D[iZ,iT] - (T[iZ,iT] * infiltParam.γ[iZ] * Sorptivity ^ 2.0) / (infiltParam.RingRadius[iZ] * Δθ)
		end
	return ∑Infilt_1D
	end # function : INFILTRATION3D_2_1D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  INFILTTRATIONη_2_3D
	# 		TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
	#		Function compute infiltration-1d from normalized infiltration
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTTRATIONη_2_3D(Infilt_η, infiltParam, iZ, K_θini, Sorptivity, T, Time_η, ΔK, Δθ)
			
			INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, T, ΔK) = K_θini * T + Infilt_η * (Sorptivity ^ 2.0) / (2.0 * ΔK)

			ΔI_η = Time_η * infiltParam.γ[iZ]

		return ∑Infilt_3D = INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, T, ΔK) + ΔI_η * (Sorptivity ^ 4.0) / (2.0 * infiltParam.RingRadius[iZ] * Δθ * (ΔK ^ 2.0))
		end # function :  INFILTTRATIONη_2_3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_η
 	# 		TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_η(∑Infilt_3D, infiltParam, iZ, K_θini, Sorptivity, T, ΔK, Δθ; ϵ=eps())
			return Infilt_η = max((2.0 * ΔK / Sorptivity ^ 2.0) * (∑Infilt_3D - K_θini * T - infiltParam.γ[iZ] * TIME_2_TIMEη(Sorptivity, T, ΔK) * (Sorptivity^4.0) / (infiltParam.RingRadius[iZ] * Δθ * 2.0* (ΔK^2.0)) ), ϵ)
		end # function INFILTRATION3D_2_η


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_2_INFILTRATION3D
	# 		COMPUTE INFILTRATION_3D FROM OPTIMIZED HYDRAULIC PARAMETERS
	# 		Solving quasiexact solution
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
		function HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, option, T; Infilt_η_Max_Start=0.5)

			Infilt_η = fill(0.0::Float64, N_Infilt[iZ])

			Sorptivity = sorptivity.SORPTIVITY(infiltParam.θini[iZ], iZ, hydroInfilt, option, option.infilt)

			Se_Ini = wrc.θ_2_Se(infiltParam.θini[iZ], iZ, hydroInfilt)

			K_θini = kunsat.Se_2_KUNSAT(option.infilt, Se_Ini, iZ, hydroInfilt)

			ΔK = hydroInfilt.Ks[iZ] - K_θini

			Δθ = hydroInfilt.θs[iZ] - infiltParam.θini[iZ]
		
			# At t=1
				∑Infilt_3D[1] = 0.0
				Infilt_η[1] = 0.0
				Infilt_η_Min = 10^-8
				Infilt_η_Max = Infilt_η_Max_Start #Since T[1] = 0

			# ~~~~~~~~~~~~~~~~~~~~
			function OF_QUASIEXACTη(Infilt_η, infiltParam, iZ, Time_η)
				Left_Term = Time_η

				Right_Term = (1.0 / (1.0 - infiltParam.β[iZ])) * (Infilt_η - log((exp(infiltParam.β[iZ] * Infilt_η) + infiltParam.β[iZ] - 1.0) / infiltParam.β[iZ]))
				
				if Right_Term < 0.0
					return OF = 10000.0 * exp(Infilt_η)
				else
					return OF = abs(Left_Term - Right_Term)
				end
			end # function OF_QUASIEXACTη ~~~~~~~~~~~~~~~~~~~~


			for iT = 2:N_Infilt[iZ] # Looping for every time step
				Time_η = TIME_2_TIMEη(Sorptivity, T[iZ,iT], ΔK)

				# Solving for Infilt_η
					Optimization = Optim.optimize(Infilt_η -> OF_QUASIEXACTη(Infilt_η, infiltParam, iZ, Time_η), Infilt_η_Min, Infilt_η_Max, Optim.GoldenSection())

					Infilt_η[iT] = Optim.minimizer(Optimization)[1]

				# Deriving the new bounds such that infiltration increases with time & the slope decreases with time
					Infilt_η_Min = Infilt_η[iT]

				# Maximum infiltration rate for T+1: (Infilt[T2] - Infilt[T1]) / (T2 - T1) which is 1 seconds
					if iT <= N_Infilt[iZ] - 1
						Infilt_η_Max = Infilt_η[iT] + (T[iZ,iT+1]- T[iZ,iT]) * (Infilt_η[iT] - Infilt_η[iT-1]) / (T[iZ,iT] - T[iZ,iT-1])
					else
						Infilt_η_Max = Infilt_η[iT] + (Infilt_η[iT] - Infilt_η[iT-1])
						
					end

				# Transforming INFILTTRATIONη to INFILTRATION3D 
					∑Infilt_3D[iZ,iT] =  INFILTTRATIONη_2_3D(Infilt_η[iT], infiltParam, iZ, K_θini, Sorptivity, T[iZ,iT], Time_η, ΔK, Δθ)
			end
		return ∑Infilt_3D
		end # function: HYDRO_2_INFILTRATION3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_QUASIEXACT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	function OF_QUASIEXACT(∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, T; W=0.8)

		Se_Ini = wrc.θ_2_Se(infiltParam.θini[iZ], iZ, hydroInfilt)
		
		K_θini = kunsat.Se_2_KUNSAT(option.infilt, Se_Ini, iZ, hydroInfilt)

		Kr_θini = K_θini / hydroInfilt.Ks[iZ]
		
		ΔK = hydroInfilt.Ks[iZ] - K_θini
		
		Δθ = hydroInfilt.θs[iZ] - infiltParam.θini[iZ]
		
		Sorptivity = sorptivity.SORPTIVITY(infiltParam.θini[iZ], iZ, hydroInfilt, option, option.infilt)

		Left_Term = zeros(Float64, N_Infilt[iZ])

		Right_Term = zeros(Float64, N_Infilt[iZ])

		iT_TransSteady = infiltOutput.iT_TransSteady_Data[iZ]
		
		Of_Penalty = 0.0 ; Of_Stead = 0.0 ; Of_Trans = 0.0
		
		for iT = 2:N_Infilt[iZ]
			Time_η = TIME_2_TIMEη(Sorptivity, T[iZ,iT], ΔK)

			Infilt_η = INFILTRATION3D_2_η(∑Infilt_Obs[iZ,iT], infiltParam, iZ, K_θini, Sorptivity, T[iZ,iT], ΔK, Δθ)

			Left_Term[iT] = Time_η

			Right_Term[iT] = (1.0 / (1.0 - infiltParam.β[iZ])) * (Infilt_η - log((exp(infiltParam.β[iZ] * Infilt_η) + infiltParam.β[iZ] - 1.0) / infiltParam.β[iZ]))
			
			if Right_Term[iT] > 0.0
				if iT <= infiltOutput.iT_TransSteady_Data[iZ]
					Of_Trans += ((Left_Term[iT]) - (Right_Term[iT])) ^ 2.0
				else
					Of_Stead += (log10(Left_Term[iT]) - log10(Right_Term[iT])) ^ 2.0
				end #  iT <= infiltOutput.iT_TransSteady_Data
			else
				Of_Penalty += 1000.0 * exp(Infilt_η)
				Right_Term[iT] = 0.0
			end #  Right_Term[iT] > 0.0
		end #  Right_Term[iT] > 0.0

	return Wof = (W * Of_Trans / Float64(iT_TransSteady-1)) + ((1.0 - W) * Of_Stead / Float64(N_Infilt[iZ] - iT_TransSteady + 1)) + Of_Penalty
	end # function: OF_INFILT_2_HYDRO
	


end # module quasiExact