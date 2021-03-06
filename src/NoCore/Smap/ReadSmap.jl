	# =============================================================
	#		module non Core : smap
	# =============================================================
	module readSmap
	   import ..tool
   	import Polynomials, DelimitedFiles
		   export SMAP, ROCKFRAGMENT_WETTABLE_STRUCT

			RockFragment_Max = 0.9

			function SMAP(IdSelect, NiZ, Path)
				println("    ~  $(Path) ~")

				# Read data
					Data = DelimitedFiles.readdlm(Path, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]
				# Sort data
					Data = sortslices(Data, dims=1)
				
				IsTopsoil, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "IsTopsoil", NiZ; N_Point_Max=1)
				IsTopsoil = 	Int64.(IsTopsoil[1:NiZ])

				Soilname, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Soilname", NiZ; N_Point_Max=1)
			
				RockClass, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "RockClass", NiZ; N_Point_Max=1)
				
				Smap_Depth, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "depth_mm", NiZ; N_Point_Max=1)

				RockFragment, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Stone_Prop", NiZ; N_Point_Max=1)

				RockFragment = min.(RockFragment_Max, RockFragment)

				Smap_RockDepth, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "RockDepth_mm", NiZ; N_Point_Max=1)

				Smap_MaxRootingDepth, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "MaxRootingDepth_mm", NiZ; N_Point_Max=1)

				Smap_SmapFH, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "SmapFH", NiZ; N_Point_Max=1)

				Smap_PermeabilityClass, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "PermeabilityClass", NiZ; N_Point_Max=1)
		
			return IsTopsoil, RockClass, RockFragment, Smap_Depth, Smap_MaxRootingDepth, Smap_PermeabilityClass, Smap_RockDepth, Smap_SmapFH, Soilname
			end  # function: SMAP


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION :SMAP
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	struct SMAP_STRUCT
		# 		Smap_Depth        ::Vector{Float64}
		# 		IsTopsoil    ::Vector{Int64}
		# 		Soilname     ::Vector{String}
		# 		RockFragment ::Vector{Float64}
		# 		RockClass    ::Vector{String}
		# 		Smap_RockDepth    ::Vector{Float64}
		# 		Smap_MaxRootingDepth ::Vector{Float64}
		# 	end

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ROCKFRAGMENT_WETTABLE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			struct ROCKFRAGMENT_WETTABLE_STRUCT
				# RockClass::Array{String}
				RockClass_Dict::Dict{String, Int64} 
				??_Rf::Array{Float64} 
				??_Rf::Array{Float64}
				N_??::Array{Int64}
				N_RockClass::Int64
				RockClass_Polynomial_Array::Array{} 
			end
			function ROCKFRAGMENT_WETTABLE(Path)
				println("    ~  $(Path) ~")
				
				# Read data
					Data = DelimitedFiles.readdlm(Path, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]
				# Sort data
					RockClass, N_RockClass = tool.readWrite.READ_HEADER_FAST(Data, Header, "RockClass")

					RockClass_Unique = unique(RockClass)
					
					N_RockClass = length(RockClass_Unique)

				# Dictionary
					RockClass_Dict = Dict("a"=>9999)
					for i=1:N_RockClass
						RockClass_Dict[RockClass_Unique[i]] = i
					end

				# Read data
					?????, N??? = tool.readWrite.READ_HEADER_FAST(Data, Header, "H[mm]")
					?????, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Theta[0-1]") 

					??_Rf = zeros(Int64,(N_RockClass, 100))
					??_Rf = zeros(Float64,(N_RockClass, 100))
					N_?? = zeros(Int64,(N_RockClass))

					iRockClass=1 ; i??=1
					for i=1:N???
						if RockClass[i] == RockClass_Unique[iRockClass]
							??_Rf[iRockClass,i??] = ?????[i]
							??_Rf[iRockClass,i??] = ?????[i]
						else
							N_??[iRockClass]  = i?? -1
							iRockClass += 1
							i?? = 1
							??_Rf[iRockClass,i??] = ?????[i]
							??_Rf[iRockClass,i??] = ?????[i]
						end
						i?? += 1
					end # for i=1:N

					N_??[iRockClass]  = i?? - 1

				RockClass_Polynomial_Array = []
				for iRockClass=1:N_RockClass
					RockClass_Polynomial = Polynomials.fit(log1p.(??_Rf[iRockClass,1:N_??[iRockClass]]), ??_Rf[iRockClass,1:N_??[iRockClass]])
					X = log1p.(??_Rf[iRockClass,1:N_??[iRockClass]])

					Coeffs = Polynomials.coeffs(RockClass_Polynomial)
				
					RockClass_Polynomial_Array = push!(RockClass_Polynomial_Array, [Coeffs])
				end

			return rfWetable = ROCKFRAGMENT_WETTABLE_STRUCT(RockClass_Dict, ??_Rf, ??_Rf, N_??, N_RockClass, RockClass_Polynomial_Array)	
			end  # function: ROCKFRAGMENT_WETTABLE
		
	end  # module: smap
	# ...........................................................