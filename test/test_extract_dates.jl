using Test
using Dates

@testset "Date Extraction (`SBAS.extract_dates`)" begin

    @testset "Single Date Paths" begin
        @test SBAS.extract_dates("path/to/20231015.geo") == [Date(2023, 10, 15)]
        @test SBAS.extract_dates("/data/insar/20210101.tif") == [Date(2021, 1, 1)]
    end

    @testset "Date Pair Paths" begin
        pair_result = SBAS.extract_dates("path/to/20231015_20231027.unwrap")
        @test pair_result == [Date(2023, 10, 15), Date(2023, 10, 27)]
        @test length(pair_result) == 2
    end

    @testset "Invalid Inputs & Edge Cases" begin
        @test isempty(SBAS.extract_dates("path/to/no_date_here.txt"))
        @test isempty(SBAS.extract_dates("path/to/2023101.geo")) # 7 digits
        @test isempty(SBAS.extract_dates(""))                   # Empty string
    end

    @testset "Type Safety" begin
        res = SBAS.extract_dates("path/to/20231015.geo")
        @test res isa Vector{Date}
        @test SBAS.extract_dates("invalid") isa Vector{Date}
    end

end
