using AntennaFieldRepresentations
using Test

@testset verbose = true "AntennaFieldRepresentations.jl" begin

    @testset "Regressions" verbose = true begin
        @testset "Freeze Tests" begin
            include(joinpath("test_Regressions", "freeze_tests.jl"))
        end
    end

    # @testset "DipoleInteractions" verbose=true begin
    #     @testset "Dipole interactions" begin
    #         include(joinpath("test_DipoleInteractions", "test_interactions.jl"))
    #     end    
    #     @testset "Dipole farfields" begin
    #         include(joinpath("test_DipoleInteractions", "test_farfields.jl"))
    #     end
    # end
    @testset "SphericalVectorModeFields" verbose = true begin
        @testset "Legendre Polynomials" begin
            # include("test_Pl.jl")     
            include(joinpath("test_SphericalVectorModeFields", "test_Plm.jl"))
            include(joinpath("test_SphericalVectorModeFields", "test_sphPlm.jl"))
        end
        @testset "Spherical Wave Functions" begin
            include(joinpath("test_SphericalVectorModeFields", "test_R_dependencies.jl"))
            include(joinpath("test_SphericalVectorModeFields", "test_Fslm.jl"))
            include(joinpath("test_SphericalVectorModeFields", "test_Kslm.jl"))
        end
        @testset "Electric and Magnetic Field" begin
            include(joinpath("test_SphericalVectorModeFields", "test_planewave.jl"))
        end
    end
    @testset "FastSpherical" verbose = true begin
        @testset "fastsphericalforward" begin
            include(joinpath("test_FastSpherical", "test_fastsphericalforward.jl"))
        end
        @testset "fastsphericalinverse" begin
            include(joinpath("test_FastSpherical", "test_fastsphericalinverse.jl"))
        end
        @testset "fastsphericaladjoint" begin
            include(joinpath("test_FastSpherical", "test_fastsphericaladjoint.jl"))
        end
    end

    # @testset "PlaneWaveRepresentations" verbose=true begin
    #     @testset "Conversions" begin
    #         include(joinpath("test_PlaneWaveRepresentations", "test_convert.jl"))
    #     end

    #     @testset "Legendre functions" begin
    #         include(joinpath("test_PlaneWaveRepresentations", "test_Pl.jl"))
    #     end

    #     @testset "Interpolation" begin
    #         include(joinpath("test_PlaneWaveRepresentations", "test_interpolation.jl"))
    #     end

    #     @testset "Operations" begin
    #         include(joinpath("test_PlaneWaveRepresentations", "test_operations.jl"))
    #     end

    #     @testset "Transmission" begin
    #         include(joinpath("test_PlaneWaveRepresentations", "test_transmission.jl"))
    #     end
    # end

    @testset verbose = true "Conversions of Field Representations" begin

        @testset "ConvertRepresentationMaps" begin
            include(joinpath("test_Conversions", "test_convertrepresentationmaps.jl"))
        end

        # @testset "Conversions Dipole -> Spherical" begin 
        #     include(joinpath("test_Conversions", "test_spherical_dipole_conversion.jl"))
        # end

        # @testset "Conversions Dipole -> PlaneWave" begin
        #     include(joinpath("test_Conversions", "test_planewave_dipole_conversion.jl"))
        # end

        # @testset "Conversions Spherical <-> PlaneWave" begin
        #     include(joinpath("test_Conversions", "test_planewave_spherical_conversion.jl"))
        # end
    end


    # @testset verbose=true "Interaction Between Field Representations" begin

    #     @testset "Interactions Dipole <-> Spherical" begin
    #         include(joinpath("test_Conversions", "test_spherical_dipole_interaction.jl"))
    #     end

    #     @testset "Interactions Dipole <-> PlaneWave" begin
    #         include(joinpath("test_Conversions", "test_planewave_dipole_interaction.jl"))
    #     end

    #     @testset "Interactions Spherical <-> PlaneWave" begin
    #         include(joinpath("test_Conversions", "test_planewave_spherical_interaction.jl"))
    #     end
    # end

    # @testset "Wacker algorithm" begin
    #     include(joinpath("test_Wacker", "test_wacker.jl"))
    # end

    # @testset verbose = true "MLFMM" begin
    #     include(joinpath("test_MLFMM", "test_dipoles_mlfmm_forward.jl"))
    #     include(joinpath("test_MLFMM", "test_dipoles_mlfmm_matrix.jl"))
    # end
    @testset "Formatting" begin
        include("test_formatting.jl")
    end

    @testset "Doctests" begin
        include("test_doctests.jl")
    end
end
