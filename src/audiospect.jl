using AxisArrays
using MetaArrays
using DSP
using JLD2
using FileIO
using SampledSignals
using Statistics
using Unitful: s

export frequencies, times, nfrequencies, ntimes, delta_t, delta_f, Δt, Δf,
  frame_length, audiospect, freq_ticks, duration, hastimes, HasTimes, HasNoTimes,
  timedim, describe_axes, Audiospect, AudiospectInv, filt

########################################
# cochlear filters

(filter::Filter)(x) = DSP.filt(filter.B,filter.A,x)

########################################
# auditory spectrogram data type and its parameters
const fixed_fs = 8000

struct TimeAxis
  Δ::typeof(1.0s)
  decay::Float64
  fs::typeof(1.0Hz)
end

struct FreqAxis{T}
  nonlinear::Float64
  octave_shift::Float64
  step::Int
  cochlear::T
end

const WithAxes{Ax} = AxisArray{<:Any,<:Any,<:Any,Ax}
const AuditorySpectrogram =
  MetaAxisArray{<:WithAxes{<:Tuple{Axis{:time},Axis{:freq}}}}

resultname(x::AuditorySpectrogram) = "Auditory Spectrogram"

num_base_freqs(x::AxisMeta) = length(x.freq.cochlear.filters)-1
num_base_freqs(x::MetaAxisArray) = num_base_freqs(getmeta(x))
nfrequencies(x) = length(frequencies(x))

function frequencies(x::AxisMeta)
  result = (440.0Hz * 2.0.^(((1:num_base_freqs(x)).-31)./24 .+
    x.freq.octave_shift))
  result[1:x.freq.step:end]
end
frequencies(as::MetaUnion{AxisArray}) = axisvalues(AxisArrays.axes(as,Axis{:freq}))[1]

ntimes(x) = length(times(x))
times(as::MetaUnion{AxisArray}) = axisvalues(AxisArrays.axes(as,Axis{:time}))[1]
times(p::AxisMeta,x::AbstractArray) = (Base.axes(x,1) .- 1) .* Δt(p)

struct HasTimes end
struct HasNoTimes end
hastimes(x::AuditorySpectrogram) = HasTimes()
hastimes(x::MetaUnion{AxisArray}) = :time ∈ axisnames(x) ? HasTimes() : HasNoTimes()
timedim(x::MetaUnion{AxisArray}) = axisdim(x,Axis{:time})

duration(x::AbstractArray) = last(times(x)) - first(times(x))
duration(x::SampleBuf) = nframes(x) / samplerate(x)

delta_t(x) = Δt(x)
delta_f(x) = Δf(x)

frame_length(params::MetaAxisLike) = frame_length(params.time)
frame_length(ax::TimeAxis) = floor(Int,ax.Δ * ax.fs)
Δt(params::MetaAxisLike) = params.time.Δ
function Δf(params::MetaAxisLike)
  @assert !isempty(params.freq.cochlear.filters)
  2^(1/24)
end
SampledSignals.samplerate(params::MetaAxisLike) = ustrip(uconvert(Hz,params.time.fs))

Δt(x) = (ts = times(x); ts[2] - ts[1])

struct Audiospect{P <: AxisMeta}
  params::P
end

asseconds(x::Number) = x*s
asseconds(x::Quantity) = uconvert(s,x)
asHz(x::Number) = x*Hz
asHz(x::Quantity) = uconvert(Hz,x)

"""
    Audiospect(Δt=10ms,freq_step=1,decay_tc=8,nonlinear=-2,octave_shift=-1)

Keyword argument constructor for an auditory spectrogram frequency
filterbank. Can be applied to a time amplitude signal using `filt`, like so:

```julia
using CorticalSpectralTemporalResponses
using FileIO

as = Audiospect()
spect = filt(as,load("testsound.wav"))
```

To compute the inverse, you can call `inv` and run `filt` on the result, like so:

```julia
save("approxsound.wav",filter(inv(as),spect)))
```
"""
function Audiospect(;delta_t=10ms,Δt=delta_t,
                          freq_step=1,decay_tc=8,nonlinear=-2,
                          octave_shift=-1)

  time = TimeAxis(asseconds(Δt),Float64(decay_tc),asHz(fixed_fs))
  freq = FreqAxis(Float64(nonlinear),Float64(octave_shift),Int(freq_step),
    cochlear[])

  recommended_length = 2^(4+freq.octave_shift)
  if frame_length(time) < recommended_length
    @warn("It's recommended that you have a frame length of at least"*
          " $recommended_length samples, but you have $(frame_length(time)).")
  end

  Audiospect(AxisMeta(time = time,freq = freq))
end
struct DefaultAudiospect
end
const audiospect=DefaultAudiospect()

function DSP.filt(f::DefaultAudiospect,x::AbstractArray,progressbar=true)
  filt(Audiospect(),x,progressbar)
end

function DSP.filt(f::Audiospect,x::AbstractArray,progressbar=true)
  @warn "Assuming sample rate of input is $(fixed_fs)."
  filter_audiospect(x,f.params,progressbar)
end

inHz(x::Number) = x
inHz(x::Quantity) = ustrip(uconvert(Hz,x))

function DSP.filt(f::Audiospect,x::AxisArray,progressbar=true)
  @assert hastimes(x) isa HasTimes
  if !(inHz(1/step(times(x))) ≈ fixed_fs)
    error("Expected samplerate of $(fixed_fs) Hz.")
  end
  filter_audiospect(x,f.params,progressbar)
end

function DSP.filt(f::Audiospect,x::SampleBuf,progressbar=true)
  if samplerate(x) != fixed_fs
    error("Expected samplerate of $(fixed_fs) Hz.")
  end
  filter_audiospect(x,f.params,progressbar)
end

####################
# the actual computation of a spectrogram
function filter_audiospect(x::AbstractArray,params::AxisMeta,progressbar=true)
  frame_len  = frame_length(params)
  N = ceil(Int,length(x) / frame_len) # of frames
  x_ = if length(x) < N*frame_len
    [x; fill(zero(eltype(x)),N*frame_len - length(x))]
  else
    x
  end
  filter_audiospect__(x_,N,params,progressbar)
end

function filter_audiospect__(x::AbstractVector{T}, N, params::MetaAxisLike,
                           progressbar=true, internal_call=false) where {T}
  step = params.freq.step
  filters = params.freq.cochlear.filters[vcat(1:step:end-1,end)]
  M = length(filters)
  Y = fill(zero(float(T)),N, M-1)
  Y_haircell = !internal_call ? nothing : fill(zero(T),length(x),M-1)

  last_haircell = x |>
    params.freq.cochlear.filters[M] |>
    x -> ion_channels(x,params) |>
    x -> haircell_membrane(x,params)

  progress = progressbar ? Progress(desc="Auditory Spectrogram: ",M-1) : nothing
  for ch = (M-1):-1:1
    # initial haircell transduction
    y,last_haircell = x |>
      filters[ch] |>
      x -> ion_channels(x,params) |>
      x -> haircell_membrane(x,params) |>
      x -> lateral_inhibition(x,last_haircell)

    # rectification and temporal integration
    Y[:,ch] = y |> rectify |> temporal_integration(params,N)

    # save the intermediate result y if this is an internal call
    if internal_call; Y_haircell[:,ch] = y end
    next!(progress)
  end

  if internal_call
    Y,Y_haircell
  else
    f = Axis{:freq}(frequencies(params))
    t = Axis{:time}(times(params,Y))

    MetaAxisArray(params,AxisArray(Y,t,f))
  end
end

########################################
# inverse of auditory spectorgram

Base.@kwdef struct AudiospectInv
  max_iterations::Int = typemax(Int)
  target_error::Float64 = 0.05
end
Base.inv(x::DefaultAudiospect;kwds...) = AudiospectInv(;kwds...)
Base.inv(x::Audiospect;kwds...) = AudiospectInv(;kwds...)
const audiospect⁻ = AudiospectInv()

function DSP.filt(f::AudiospectInv,y_in::AuditorySpectrogram,
  progressbar=true)
  max_iterations = f.max_iterations; target_error = f.target_error

  @assert(max_iterations < typemax(Int) || target_error < Inf,
          "No stopping criterion specified (max_iterations or target_error).")
  M = length(y_in.freq.cochlear.filters)

  # expand y to include all frequencies
  y = zeros(eltype(y_in),size(y_in,1),M-1)
  f_ixs = minimum(frequencies(y_in)) .<= frequencies(y_in) .<= maximum(frequencies(y_in))
  # NOTE: this does not currently support frequency decimation
  # because I don't yet need it
  @assert sum(f_ixs) == size(y_in,2) "Unxpected frequency resolution."
  y[:,f_ixs] = y_in

  # generate initial guess
  x = inv_guess(y_in,y)
  N = ceil(Int,length(x) / frame_length(y_in)) # of frames

  # iteration setup
  min_err = Inf
  min_x = x

  prog = if !progressbar
    nothing
  elseif max_iterations < typemax(Int)
    Progress(max_iterations,"Inverting Spectrogram: ")
  else
    ProgressThresh(target_error,"Inverting Spectrogram: ")
  end

  for iteration in 1:max_iterations
    if min_err < target_error; break end
    y_in.freq.nonlinear == 0 && standardize!(x)

    ŷ,ŷ_haircell = filter_audiospect__(x,N,y_in,false,true)
    x = match_x(y_in,x,y,ŷ,ŷ_haircell)

    err = relative_error!(y,ŷ)

    if err < min_err
      min_x = x
      min_err = err
    elseif err-1 > min_err
      randomized_reset!(x)
    end

    x .*= 1.01

    if progressbar
      if max_iterations < typemax(Int)
        ProgressMeter.next!(prog;showvalues =
                            [(:error,string(100round(min_err,digits=4),"%"))])
      else
        ProgressMeter.update!(prog,min_err)
      end
    end
  end
  ProgressMeter.finish!(prog)

  SampleBuf(min_x,samplerate(y_in))
end

########################################
# visualization utility
function freq_ticks(as)
  a = minimum(frequencies(as))
  b = maximum(frequencies(as))
  step = 0.25

  helper(step) = round.(filter(f -> a <= f*Hz <= b,1000*2.0.^(-3:step:2)),
                        digits=-1)
  fbreaks = helper(step)
  while length(fbreaks) > 7
    fbreaks = helper(step *= 2)
  end

  fs = ustrip.(uconvert.(Hz,frequencies(as)))

  findices = mapslices(abs.(fbreaks .- fs'),dims=2) do row
    _, i = findmin(row)
    i
  end

  fbreaks,findices
end


################################################################################
# private helper functions

function ion_channels(x,params::MetaAxisLike)
  if params.freq.nonlinear > 0
    1.0 ./ (1.0 .+ exp.(.-x./params.freq.nonlinear))
  elseif params.freq.nonlinear == 0
    Float64.(x .> 0.0)
  elseif params.freq.nonlinear == -1
    max.(x,0.0)
  elseif params.freq.nonlinear == -2
    x
    # TODO: implement halfregu
  else
    error("Non linear factor of $fac not supported")
  end
end

function haircell_membrane(x,params)
  if params.freq.nonlinear  != -2
    β = exp(-1/((1//2)*2^(4+params.freq.octave_shift)))
    filt([1.0],[1.0; -β],x)
  else
    x
  end
end

rectify(x) = max.(x,0)

function temporal_integration(params::MetaAxisLike,N)
  frame_len = frame_length(params)
  if !iszero(params.time.decay)
    α = exp(-1/(params.time.decay*2^(4+params.freq.octave_shift)))
    function(x)
      y = filt([1.0],[1.0; -α],x)
      @views y[frame_len*(1:N)]
    end
  else    # short-term average
    frame_len == 1 ? identity : y -> mean(reshape(y, frame_len, N),1)'
  end
end

function lateral_inhibition(x,last_haircell)
  new_x = x .- last_haircell
  new_x,x
end

function inv_guess(params::MetaAxisLike,y::AbstractMatrix)
  # the initial guess only uses the first 48 channels
  f = ustrip.(uconvert.(Hz,frequencies(params)))[1:48]
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))
  t = steps ./ ustrip(uconvert(Hz,params.time.fs))

  x = dropdims(standardize!(sum(cos.(2π.*f'.*t) .* view(y,indices,1:48),dims=2)),
               dims=2)
end

function match_x(params::MetaAxisLike,x,y,ŷ,ŷ_haircell)
  M = length(params.freq.cochlear.filters)
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))

  ratios = map(ŷ,y) do ŷ,y
    !iszero(ŷ) ? y ./ ŷ : !iszero(y) ? 2 : 1
  end

  x .= 0
  ch_norm = params.freq.cochlear.norm
  for ch in 1:M-1
    if params.freq.nonlinear == -2
      y1 = ŷ_haircell[:,ch].*view(ratios,indices,ch)
    else
      y1 = ŷ_haircell[:,ch]
      posi = find(y1 .>= 0)
      y1[posi] .*= view(ratios,indices[posi],ch)

      negi = find(y1 .< 0)
      y1[negi] .*= maximum(y1[posi]) / -minimum(y1[negi])
    end

    x .+= reverse(params.freq.cochlear.filters[ch](reverse(y1))) ./ ch_norm
  end

  x
end

function randomized_reset!(x)
  x .= sign.(x) .+ rand(size(x))
  x .-= mean(x)
  x ./= std(x)
  x
end

function relative_error!(y,ŷ)
  y2sum = sum(y.^2)
  ymean = mean(y)

  ŷ .*= ymean/mean(ŷ)
  sum((ŷ .- y).^2) ./ y2sum
end

function standardize!(x)
  x .-= mean(x)
  x ./= std(x)
end