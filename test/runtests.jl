using Test
using AxisArrays
using SampledSignals
using Statistics
using Unitful
using DSP
using PlotAxes
using CorticalSpectralTemporalResponses

x = SampleBuf(sin.(2π .* 1000 .* range(0,stop=1,length=8000)),8000)
X = filt(audiospect,x)

as = 

err = 0.05
x_hat = filt(inv(audiospect,target_error=err,max_iterations=1000),X)
as_x_hat = filt(audiospect,x_hat)

@testset "CorticalSpectralTemporalResponses" begin

@testset "Spectrogram" begin
  @test_throws ErrorException filt(audiospect,SampleBuf(collect(1:10),4000))
  @test mean(X[:,0.9kHz ..1.1kHz]) > mean(X[:,1.9kHz .. 2.1kHz])
  @test sum((as_x_hat .* (mean(X) / mean(as_x_hat)) .- X).^2) ./
    sum(X.^2) <= err
  X_small = filt(Audiospect(freq_step=2),x)
  @test nfrequencies(X_small) == 64
  @test mean(X_small[:,0.9kHz ..1.1kHz]) > mean(X_small[:,1.9kHz .. 2.1kHz])
  msg = r"It's recommended that you have a frame length of at least 8.0 samples"
  @test_logs((:warn,msg), Audiospect(Δt = 0.5ms))

  data = AxisArray(rand(2000),Axis{:time}(range(0,step=1/8000,length=2000)))
  @test filt(audiospect,data) isa CorticalSpectralTemporalResponses.AuditorySpectrogram
  data = AxisArray(rand(2000),Axis{:time}(range(0,step=1/3000,length=2000)))
  @test_throws ErrorException("Expected samplerate of 8000 Hz.") filt(audiospect,data) 

  @test eltype(filt(audiospect,collect(1:10))) == float(Int)

  @test length(frequencies(X)) > 0
  @test minimum(frequencies(X)) > 0Hz
  @test maximum(frequencies(X)) > 3kHz
  @test nfrequencies(X) == 128
  @test nfrequencies(X[:,1:20]) == 20
  @test abs(maximum(times(X)) - duration(x)*s) < 2Δt(X)
  @test Δt(X) ≈ times(X)[2] - times(X)[1]
  @test Δf(X) ≈ frequencies(X)[2]/frequencies(X)[1]
  @test 3 < length(freq_ticks(X)[1]) < 7
  @test floor.(freq_ticks(X)[1]./10) == freq_ticks(X)[1]./10
end

# TODO: add tests for plotting cortical data as well
scalef = scalefilter()
S_cr = filt(scalef,X)
S_X̂ = filt(inv(scalef),S_cr)

ratef = ratefilter()
R_cr = filt(ratef,X)
R_X̂ = filt(inv(ratef),R_cr)

cort = cortical()
cr = filt(cort,X)

@testset "Single dimension cortical model" begin
  @test mean(abs,S_cr[:,:,0.9kHz ..1.1kHz]) >
    mean(abs,S_cr[:,:,1.9kHz .. 2.1kHz])
  @test mean(abs,R_cr[:,:,0.9kHz ..1.1kHz]) >
    mean(abs,R_cr[:,:,1.9kHz .. 2.1kHz])
  @test quantile(vec((S_X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test quantile(vec((R_X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test isempty(setdiff(scales(S_cr),default_scales))
  @test isempty(setdiff(rates(R_cr),default_rates))
  @test freq_ticks(S_cr) == freq_ticks(X)
  @test freq_ticks(R_cr) == freq_ticks(X)
end

@testset "Multi dimensional cortical model" begin
  X̂ = filt(inv(cort),cr)

  # TODO: fix tests
  @test mean(abs,cr[:,:,:,0.9kHz ..1.1kHz]) >
    mean(abs,cr[:,:,:,1.9kHz .. 2.1kHz])
  @test quantile(vec((X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test isempty(setdiff(scales(cr),default_scales))
  @test isempty(setdiff(rates(cr),default_rates))
  @test freq_ticks(cr) == freq_ticks(X)
end

@testset "Position and keyword arugment support" begin
  sc = [1cycoct, 2cycoct]
  rt = [1Hz, 2Hz]

  @test scalefilter(sc).data ==
    scalefilter(scales=sc).data
  @test ratefilter(rt).data == ratefilter(rates=rt).data
  cr1 = cortical(sc, rt)
  cr2 = cortical(scales=sc, rates=rt)
  @test cr1.scales.data == cr2.scales.data
  @test cr1.rates.data == cr2.rates.data

  msg = "Cannot specify both a position and keyword value for scales."
  @test_throws ErrorException(msg) scalefilter(sc.+1cycoct,scales=sc)
  msg = "Cannot specify both a position and keyword value for rates."
  @test_throws ErrorException(msg) ratefilter(rt.+1Hz,rates=rt)
  msg = "Cannot specify both a position and keyword value for scales."
  @test_throws ErrorException(msg) cortical(sc.+1cycoct,scales=sc)
  msg = "Cannot specify both a position and keyword value for rates."
  @test_throws ErrorException(msg) cortical(default_scales,rt.+1Hz,rates=rt)
end

@testset "Handles implicit and explicit unit." begin
  @test scalefilter([1cycoct]).data == scalefilter([1]).data
  @test ratefilter([1Hz]).data == ratefilter([1]).data
end

@testset "Data is plottable" begin
  @test size(PlotAxes.asplotable(X)[1],1) == length(X)
  @test size(PlotAxes.asplotable(S_cr,quantize=(1000,1000,20,20))[1],1) == length(S_cr)
  @test size(PlotAxes.asplotable(R_cr,quantize=(1000,1000,20,20))[1],1) == length(R_cr)
  @test size(PlotAxes.asplotable(cr,quantize=(1000,1000,20,20))[1],1) == length(cr)

  @test size(PlotAxes.asplotable(X,quantize=(5,5,5,5))[1],1) < length(X)
  @test size(PlotAxes.asplotable(S_cr,quantize=(5,5,5,5))[1],1) < length(S_cr)
  @test size(PlotAxes.asplotable(R_cr,quantize=(5,5,5,5))[1],1) < length(R_cr)
  @test size(PlotAxes.asplotable(cr,quantize=(5,5,5,5))[1],1) < length(cr)

  @test logrange(666) == log(666)
  @test PlotAxes.fn_prefix(logrange) == "log"
  @test logrange.(1:10) == log.(1:10)
end

@testset "Data is printed in console" begin
  iobuf = IOBuffer()
  display(TextDisplay(iobuf), X)
  @test length(String(take!(iobuf))) > 0

  iobuf = IOBuffer()
  display(TextDisplay(iobuf), S_cr)
  @test length(String(take!(iobuf))) > 0

  iobuf = IOBuffer()
  display(TextDisplay(iobuf), R_cr)
  @test length(String(take!(iobuf))) > 0

  iobuf = IOBuffer()
  display(TextDisplay(iobuf), cr)
  @test length(String(take!(iobuf))) > 0
end

@testset "Repeated axis cortical model" begin

  ratef = ratefilter(axis=:rateslow)
  crR = filt(ratef,cr)
  @test ndims(crR) == 5

  ratef = ratefilter()
  @test_throws(ArgumentError("axis name :rate is used more than once"),
    filt(ratef,cr))
end

end