# =============================================================
#		MODULE: psdOpt
# =============================================================
module psdOpt
	import ..psdFunc

	# =========================================
	#       PSD_RUN_ALLMODEL
	# 		THIS WILL RUN FOR ALL MODELS
	# =========================================
		function PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)
			
			θ_Rpart = fill(0.0::Float64, (NiZ, N_Psd_Max))
			Ψ_Rpart = fill(0.0::Float64, (NiZ, N_Psd_Max))

			for iZ = 1:NiZ
				θ_Rpart[iZ,1:N_Psd[iZ]], Ψ_Rpart[iZ,1:N_Psd[iZ]] = psdFunc.PSD_MODEL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], option, param, paramPsd)
			end # for iZ = 1:NiZ

		return θ_Rpart, Ψ_Rpart
		end # function PSD_RUN_ALLMODEL

	
	# =============================================================
	#		MODULE: imp
	# =============================================================
	module imp
		import ...stats, ...wrc, ...psdFunc, ..psdOpt
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		imp FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)
				θ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
				
				for iZ = 1:NiZ
					if !(option.psd.∑Psd_2_ξ1)
					
						SearchRange = [ (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], paramPsd, hydro, option, param; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
						; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

						paramPsd.∑Psd_2_ξ2_β1[iZ] = BlackBoxOptim.best_candidate(Optimization)[1]
						paramPsd.∑Psd_2_ξ2_β2[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
						paramPsd.Subclay[iZ]      = BlackBoxOptim.best_candidate(Optimization)[3]

					elseif option.psd.∑Psd_2_ξ1	
						SearchRange =  [(param.psd.imp.ξ1_Min, param.psd.imp.ξ1_Max), (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], paramPsd, hydro, option, param; ξ1 = P[1], ∑Psd_2_ξ2_β1 = P[2] ,∑Psd_2_ξ2_β2 = P[3], Subclay = P[4])
						; SearchRange = SearchRange, NumDimensions = 4, TraceMode = :silent)

						paramPsd.ξ1[iZ]           = BlackBoxOptim.best_candidate(Optimization)[1]
						paramPsd.∑Psd_2_ξ2_β1[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
						paramPsd.∑Psd_2_ξ2_β2[iZ] = BlackBoxOptim.best_candidate(Optimization)[3]
						paramPsd.Subclay[iZ]      = BlackBoxOptim.best_candidate(Optimization)[4]
					end # if option.psd.
				end	# for iZ = 1:NiZ
				
				# Compute the optimal values
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)

				# Optimal ξ2
				for iZ = 1:NiZ
					paramPsd.ξ2[iZ] = psdFunc.imp.∑PSD_2_ξ2(∑Psd[iZ,paramPsd.∑Psd_2_ξ2_Size[iZ]], param; ∑Psd_2_ξ2_β1=paramPsd.∑Psd_2_ξ2_β1[iZ], ∑Psd_2_ξ2_β2=paramPsd.∑Psd_2_ξ2_β2[iZ])
				end

				# Statistics
				paramPsd.Nse, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NSE_θΨ(hydro, N_Psd, NiZ, option.psd, θ_Rpart, Ψ_Rpart)

				return paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		imp FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, option, θs_Psd, θr_Psd, param, paramPsd, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		imp FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, NiZ, θs_Psd, θr_Psd, param, paramPsd, hydro; ξ1 = param.psd.imp.ξ1, ∑Psd_2_ξ2_β1 = param.psd.imp.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2 = param.psd.imp.∑Psd_2_ξ2_β2, Subclay = param.psd.imp.Subclay)

						θ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))

						for iZ = 1:NiZ
							paramPsd.ξ1[iZ]           = ξ1
							paramPsd.∑Psd_2_ξ2_β1[iZ] = ∑Psd_2_ξ2_β1
							paramPsd.∑Psd_2_ξ2_β2[iZ] = ∑Psd_2_ξ2_β2
							paramPsd.Subclay[iZ]      = Subclay
						end

						Of = 0.0
						for iZ = 1:NiZ
							θ_Rpart[iZ,1:N_Psd[iZ]], Ψ_Rpart[iZ,1:N_Psd[iZ]] = psdFunc.PSD_MODEL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], option, param, paramPsd)

							θΨ =  fill(0.0::Float64, (N_Psd[iZ]))
							for iRpart = 1:N_Psd[iZ]
								# Observed data
								θΨ[iRpart] = wrc. Ψ_2_θDual(option.psd,Ψ_Rpart[iZ, iRpart], iZ, hydro)
							end

							Of += stats.NSE_MINIMIZE(θΨ[1:N_Psd[iZ]], θ_Rpart[iZ,1:N_Psd[iZ]]; Power=2)
						end # for iZ = 1:NiZ
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


				if !(option.psd.∑Psd_2_ξ1)	
					SearchRange = [ (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, NiZ, θs_Psd, θr_Psd, param, paramPsd, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
					; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

					for iZ=NiZ
						paramPsd.∑Psd_2_ξ2_β1[iZ] = BlackBoxOptim.best_candidate(Optimization)[1]
						paramPsd.∑Psd_2_ξ2_β2[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
						paramPsd.Subclay[iZ]      = BlackBoxOptim.best_candidate(Optimization)[3]
					end

				elseif option.psd.∑Psd_2_ξ1
					SearchRange =  [(param.psd.imp.ξ1_Min, param.psd.imp.ξ1_Max), (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, NiZ, θs_Psd, θr_Psd, param, paramPsd, hydro; ξ1 = P[1], ∑Psd_2_ξ2_β1 = P[2] ,∑Psd_2_ξ2_β2 = P[3], Subclay = P[4])
					; SearchRange = SearchRange, NumDimensions = 4, TraceMode = :silent)

					for iZ=NiZ
						paramPsd.ξ1[iZ]           = BlackBoxOptim.best_candidate(Optimization)[1]
						paramPsd.∑Psd_2_ξ2_β1[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
						paramPsd.∑Psd_2_ξ2_β2[iZ] = BlackBoxOptim.best_candidate(Optimization)[3]
						paramPsd.Subclay[iZ]      = BlackBoxOptim.best_candidate(Optimization)[4]
					end
				end # if option.psd.

				# Compute the optimal values
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)

				# Optimal ξ2
				for iZ = 1:NiZ
					paramPsd.ξ2[iZ] = psdFunc.imp.∑PSD_2_ξ2(∑Psd[iZ,paramPsd.∑Psd_2_ξ2_Size[iZ]], param; ∑Psd_2_ξ2_β1=paramPsd.∑Psd_2_ξ2_β1[iZ], ∑Psd_2_ξ2_β2=paramPsd.∑Psd_2_ξ2_β2[iZ])
				end

				# Statistics
				paramPsd.Nse, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NSE_θΨ(hydro, N_Psd, NiZ, option.psd, θ_Rpart, Ψ_Rpart)
				

				return paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		imp FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, paramPsd, hydro, option, param; ξ1 = paramPsd.ξ1[iZ], ∑Psd_2_ξ2_β1 = paramPsd.∑Psd_2_ξ2_β1[iZ], ∑Psd_2_ξ2_β2 = paramPsd.∑Psd_2_ξ2_β2[iZ], Subclay = paramPsd.Subclay[iZ])

					paramPsd.ξ1[iZ]           = ξ1
					paramPsd.∑Psd_2_ξ2_β1[iZ] = ∑Psd_2_ξ2_β1
					paramPsd.∑Psd_2_ξ2_β2[iZ] = ∑Psd_2_ξ2_β2
					paramPsd.Subclay[iZ]      = Subclay

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(iZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd)

					# For every Ψ_Rpart
					Of = 0.0
					θΨ = zeros(Float64, N_Psd)
					for iRpart = 1:N_Psd
						# Observed data
						θΨ[iRpart] = wrc. Ψ_2_θDual(option.psd, Ψ_Rpart[iRpart], iZ, hydro)
					end

					Of = stats.NSE_MINIMIZE(θΨ[1:N_Psd], θ_Rpart[1:N_Psd]; Power=2)
				end  # function OF_SINGLE_SOIL
		
	end  # module imp
	# ............................................................


	# =============================================================
	#		MODULE: chang
	# =============================================================
	module chang
	 	import ...stats, ...wrc, ....psdFunc, ..psdOpt
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		chang FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, paramPsd, hydro)
				θ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
				
				for iZ = 1:NiZ
					SearchRange = [ (param.psd.chan.ξ1_Min, param.psd.chan.ξ1_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], paramPsd, hydro, option, param; ξ1=P[1])
					; SearchRange = SearchRange, NumDimensions=1, TraceMode = :silent)

					paramPsd.ξ1[iZ] = BlackBoxOptim.best_candidate(Optimization)[1]
				end	# for iZ = 1:NiZ
				
				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)

				# Statistics
				paramPsd.Nse, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NSE_θΨ(hydro, N_Psd, NiZ, option.psd, θ_Rpart, Ψ_Rpart)

				return paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		chang FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, option, θs_Psd, θr_Psd, param, paramPsd, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		chang FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, NiZ, θs_Psd, θr_Psd, param, paramPsd, hydro; ξ1 = paramPsd.ξ1)
						θ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (NiZ, N_Psd_Max))

						for iZ = 1:NiZ
							paramPsd.ξ1[iZ] = ξ1
						end

						Of = 0.0
						for iZ = 1:NiZ
							θ_Rpart[iZ,1:N_Psd[iZ]], Ψ_Rpart[iZ,1:N_Psd[iZ]] = psdFunc.PSD_MODEL(iZ, Psd[iZ,1:N_Psd[iZ]], ∑Psd[iZ,1:N_Psd[iZ]], Rpart[iZ,1:N_Psd[iZ]], N_Psd[iZ], θs_Psd[iZ], θr_Psd[iZ], option, param, paramPsd)

							θΨ = zeros(Float64, N_Psd[iZ])
							for iRpart = 1:N_Psd[iZ]
								# Observed data
								θΨ[iRpart] = wrc.Ψ_2_θDual(option.psd,Ψ_Rpart[iZ, iRpart], iZ, hydro)
							end

							Of += stats.NSE_MINIMIZE(θΨ[1:N_Psd[iZ]], θ_Rpart[iZ,1:N_Psd[iZ]]; Power=2)
						end # for iZ = 1:NiZ
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				SearchRange = [ (param.psd.chang.ξ1_Min, param.psd.chang.ξ1_Max) ]

				Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, NiZ, θs_Psd, θr_Psd, param, paramPsd, hydro; ξ1 = P[1])
				; SearchRange = SearchRange, NumDimensions=1, TraceMode = :silent)

				for iZ=1:NiZ
					paramPsd.ξ1[iZ] = BlackBoxOptim.best_candidate(Optimization)[1]
				end

				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd, hydro)

				# Statistics
				paramPsd.Nse, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NSE_θΨ(hydro, N_Psd, NiZ, option.psd, θ_Rpart, Ψ_Rpart)

				return paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		chang FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, paramPsd, hydro, option, param; ξ1 = paramPsd.ξ1[iZ])

					paramPsd.ξ1[iZ] = ξ1

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(iZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, option, param, paramPsd)

					# For every Ψ_Rpart
					θΨ = zeros(Float64, N_Psd)
					for iRpart = 1:N_Psd
						# Observed data
						θΨ[iRpart] = wrc.Ψ_2_θDual(option.psd,Ψ_Rpart[iRpart], iZ, hydro)
					end

					Of = stats.NSE_MINIMIZE(θΨ[1:N_Psd], θ_Rpart[1:N_Psd]; Power=2)
				end  # function OF_SINGLE_SOIL
		
		
	end  # module chang
	# ............................................................


end  # module psdOpt
# ............................................................