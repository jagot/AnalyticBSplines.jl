using AnalyticBSplines
using Test

import AnalyticBSplines: collapse!

@testset "Intervals" begin
    @testset "Subsets" begin
        @test Interval(1,2) ⊂ Interval(0,3)
        @test Interval(1,2) ⊆ Interval(0,3)
        @test !(Interval(1,2) ⊂ Interval(1,2))
        @test Interval(1,2) ⊆ Interval(1,2)
    end

    @testset "Set differences" begin
        @test Interval(0,1) \ Interval(2,3) == ([Interval(0,1)],[Interval(2,3)])
        @test Interval(0,3) \ Interval(1,2) == ([Interval(0,1),Interval(2,3)],[])
        @test Interval(0,3) \ Interval(-1,2) == ([Interval(2,3)],[Interval(-1,0)])
        @test Interval(0,3) \ Interval(1,4) == ([Interval(0,1)],[Interval(3,4)])
        @test Interval(0,3) \ Interval(0,3) == ([],[])
        @test Interval(0,3) \ Interval(0,4) == ([],[Interval(3,4)])
        @test Interval(0,3) \ Interval(-1,3) == ([],[Interval(-1,0)])
        @test Interval(0,3) \ Interval(-1,4) == ([],[Interval(-1,0),Interval(3,4)])
        @test Interval(0,3) \ Interval(1//2,2) == ([Interval(0,1//2),Interval(2//1,3)],[])
    end

    @testset "Merge intervals" begin
        @test collapse!([Interval(0,1), Interval(1,2)]) == [Interval(0,2)]
        @test collapse!([Interval(1,2), Interval(0,1)]) == [Interval(0,2)]
        @test collapse!([Interval(0,2), Interval(1,3)]) == [Interval(0,3)]
        @test collapse!([Interval(0,1), Interval(2,3)]) == [Interval(0,1), Interval(2,3)]
    end
end

@testset "Piecewise polynomials" begin
    @testset "Addition" begin
        pp1 = PPoly([
            (0,1) => [1,2],
            (1//1,2) => [3],
        ])
        pp2 = PPoly([
            (0,1) => [1,2],
            (3//2,3) => [3,0,2],
        ])
        pp3 = PPoly([
            (0,1) => [2,4],
            (1,3//2) => [3],
            (3//2,2) => [6,0,2],
            (2,3) => [3,0,2]
        ])
        @test pp1 + pp2 == pp3

        @test PPoly([
            (0,1//2) => [0,2]
        ]) + PPoly([
            (1//2,1) => [2,-2]
        ]) == PPoly([
            (0,1//2) => [0,2],
            (1//2,1) => [2,-2]
        ])

        @test PPoly([
            (0,5) => [0,2]
        ]) + PPoly([
            (1,2) => [2,-2],
            (3,4) => [8,0,7]
        ]) == PPoly([
            (0,1) => [0,2],
            (1,2) => [2],
            (2,3) => [0,2],
            (3,4) => [8,2,7],
            (4,5) => [0,2]
        ])
    end
end
