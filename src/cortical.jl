using AxisArrays
using FFTW

export rates, scales, nrates, nscales, default_rates, default_scales,
  cortical, cycoct, co, scalefilter, ratefilter

# re-express the spectral and cortical dimensions to have meta data specific to
# an axis and then make it possible to add an axis to an existing array (maybe
# have a flag to allow multiple axes of the same type)

@dimension Sc "Sc" Scale
@refunit cycoct "cyc/oct" CyclesPerOct Sc false

# NOTE: the `low` and `high` fields are used to determine which filters should
# be low- and high-pass, rather than band-pass
struct ScaleAxis
  low::typeof(1.0cycoct)
  high::typeof(1.0cycoct)
end

struct RateAxis
  low::typeof(1.0Hz)
  high::typeof(1.0Hz)
end

rates(x::MetaUnion{AxisArray}) =
  axisvalues(AxisArrays.axes(x,Axis{:rate}))[1]
nrates(x) = length(rates(x))

scales(x::MetaUnion{AxisArray}) =
  axisvalues(AxisArrays.axes(x,Axis{:scale}))[1]
nscales(x) = length(scales(x))

const default_rates = sort([-2 .^ (1:0.5:5); 2 .^ (1:0.5:5)]).*Hz
const default_scales = (2 .^ (-2:0.5:3)).*cycoct
const spect_rate = 24

struct DefaultValue
end
const default = DefaultValue()
withdefault(x,default) = x
withdefault(::DefaultValue,default) = default
function checkmatch(pos,kwd,name) 
  if pos !== kwd
    error("Cannot specify both a position and keyword value for $name.")
  end
end
checkmatch(::DefaultValue,x,str) = nothing

ascycoct(x::Number) = x*cycoct
ascycoct(x::Quantity) = uconvert(cycoct,x)
ascycoct(x::Quantity{<:Any,Sc}) = x

abstract type CorticalFilter
end
abstract type CorticalFilterInv
end

function DSP.filt(filter::CorticalFilter,y::MetaAxisArray)
  @assert all(∈(axisnames(y)),(:time,:freq))
  if axisname(filter) ∈ y
    ax = axisname(filter)
    error("Input already has an axis named `$ax`. If you intended ",
          "to add a new dimension, you will have to change the name of the ",
          "axis. When you define the filter you can specify the axis name ",
          "using the `axis` keyword argument.")
  end

  fir = FIRFiltering(y,Axis{fromaxis(filter)})
  cr = initfilter(y,filter)
  for (I,H) in list_filters(fir,cr,filter)
    cr[I] = view(apply(fir,H),Base.axes(y)...)
  end

  cr
end

function DSP.filt(cinv::CorticalFilterInv,cr::MetaAxisArray)
  z_cum = FFTCum(cr,axisnames(cinv))

  filters = list_filters(z_cum,cr,cinv)

  inner_dims = size(cr)[2:end-1]
  if length(inner_dims) != ndims(cinv)
    error("When computing the inverse you must invert all scale/rate "*
          "dimensions; partial inverses are not supported.")
  end

  for (I,filter) in filters
    addfft!(z_cum,view(cr,:,I,:),filter)
  end

  removed = removeaxes(getmeta(cr),axisnames(cinv)...)
  norm = normalize!(z_cum,cr,cinv.norm)
  MetaAxisArray(removed,AxisArray(norm,AxisArrays.axes(cr)[[1,end]]))
end

# rate filters
struct TimeRateFilter <: CorticalFilter
  data::Vector{typeof(1.0Hz)}
  bandonly::Bool
  axis::Symbol
end
axisname(x::TimeRateFilter) = x.axis
fromaxis(x::TimeRateFilter) = :time
AxisArrays.axisnames(x::TimeRateFilter) = (x.axis,)
Base.length(x::TimeRateFilter) = length(x.data)

"""
    ratefilter(rates,bandonly=false,axis=:rate)

Define a temporal-rate filter bank at the given rates (assumed to be in
Hz by default, but can be specified in a different unit via `Unitful`).
The rates can be either a positional or keyword argument. To compute
use inverse of the filter-bank call `inv`. Uses default_rates
if the rates are unspecified.

Apply the filter bank using `filt`, like so:

```julia
using FileIO
using CorticalSpectralTemporalResponses

x = load("testsound.wav")
X = filt(audiospect,x)
ratef = ratefilter()
R_cr = filt(ratef,X)
```

"""
function ratefilter(r=default;rates=withdefault(r,default_rates),
  bandonly=false,axis=:rate)
  checkmatch(r,rates,"rates")
  rates = asHz.(rates)

  if axis != :rate && !occursin("rate",string(axis))
    error("Rate axis name `$axis` must contain the word 'rate'.")
  end
  TimeRateFilter(rates,bandonly,axis)
end
list_filters(fir,cr,filter::TimeRateFilter) =
  ((Axis{axisname(filter)}(i), HR)
   for (i,HR) in enumerate(rate_filters(fir,cr,axisname(filter))))

# inverse of rate filters
struct TimeRateFilterInv <: CorticalFilterInv
  rates::TimeRateFilter
  norm::Float64
end
Base.ndims(x::TimeRateFilterInv) = 1
axisname(x::TimeRateFilterInv) = axisname(x.rates)
AxisArrays.axisnames(x::TimeRateFilterInv) = axisnames(x.rates)
Base.inv(rates::TimeRateFilter;norm=0.9) = TimeRateFilterInv(rates,norm)
list_filters(z_cum,cr,rateinv::TimeRateFilterInv) =
  ((i, HR)
   for (i,HR) in enumerate(rate_filters(z_cum,cr,axisname(rateinv),use_conj=true)))

# scale filters
struct FreqScaleFilter <: CorticalFilter
  data::Vector{typeof(1.0cycoct)}
  bandonly::Bool
  axis::Symbol
end
axisname(x::FreqScaleFilter) = x.axis
fromaxis(x::FreqScaleFilter) = :freq
AxisArrays.axisnames(x::FreqScaleFilter) = (x.axis,)
list_filters(fir,cs,scales::FreqScaleFilter) = 
  ((Axis{axisname(scales)}(i), [HS; zero(HS)]') 
   for (i,HS) in enumerate(scale_filters(fir,cs,axisname(scales))))

"""
    scalefilter(scales,bandonly=false,axis=:rate)

Define a spectral-scale filter bank at the given scales (assumed to be in
cycles/octave). The scales can be either a positional or keyword argument. To
compute use inverse of the filter-bank call `inv`. Uses default_scales if the 
scales are unspecified.

Apply the filter bank using `filt`, like so:

```julia
using FileIO
using CorticalSpectralTemporalResponses

x = load("testsound.wav")
X = filt(audiospect,x)
scalef = scalefilter()
S_cr = filt(scalef,X)
```

"""
function scalefilter(s=default;scales=withdefault(s,default_scales),
    bandonly=false,axis=:scale)
  checkmatch(s,scales,"scales")
  scales = ascycoct.(scales)
  
  if axis != :scale && !occursin("scale",string(axis))
    error("Scale axis name `$axis` must contain the word 'scale'.")
  end
  FreqScaleFilter(scales,bandonly,axis)
end

# inverse of scale filters
struct FreqScaleFilterInv <: CorticalFilterInv
  scales::FreqScaleFilter
  norm::Float64
end
Base.ndims(x::FreqScaleFilterInv) = 1
axisname(x::FreqScaleFilterInv) = axisname(x.scales)
AxisArrays.axisnames(x::FreqScaleFilterInv) = (axisname(x.scales),)
Base.inv(scales::FreqScaleFilter;norm=0.9) = FreqScaleFilterInv(scales,norm)
list_filters(z_cum,cr,scaleinv::FreqScaleFilterInv) =
  ((i, [HS; zero(HS)]') 
   for (i,HS) in enumerate(scale_filters(z_cum,cr,axisname(scaleinv))))

# combination of both scales and rates
struct ScaleRateFilter <: CorticalFilter
  scales::FreqScaleFilter
  rates::TimeRateFilter
end

"""
    cortical(scales,rates,bandonly=false,axis=:rate)

Define a joint scale-rate filter bank at the given scales and rates (these
are assumed to be in cycles per octave and Hz respectively, but can be
specified with different units via `Unitful`). The scales and rates can be
either both as positional arguments or both as keyword arguments. To compute
the inverse of the filter-bank call `inv`.

Apply the filter bank using `filt`, like so:

```julia
using FileIO
using CorticalSpectralTemporalResponses

x = load("testsound.wav")
X = filt(audiospect,x)
cortf = cortical()
cr = filt(cortf,X)
```
"""
function cortical(s=default,r=default;
    scales=withdefault(s,default_scales),
    rates=withdefault(r,default_rates),
    bandonly=false,axes=(:scale,:rate))

    checkmatch(s,scales,"scales")
    checkmatch(r,rates,"rates")
    scales = ascycoct.(scales)
    rates = asHz.(rates)

    ScaleRateFilter(scalefilter(scales,bandonly=bandonly,axis=axes[1]),
                    ratefilter(rates,bandonly=bandonly,axis=axes[2]))
end
cortical(scales::FreqScaleFilter,rates::TimeRateFilter) = 
  ScaleRateFilter(scales,rates)

function DSP.filt(cort::ScaleRateFilter,cr::MetaAxisArray,progresbar=true)
  filt(cort.rates,filt(cort.scales,cr))
end

# inverse of scale-rate filters
struct ScaleRateFilterInv <: CorticalFilterInv
  scales::FreqScaleFilterInv
  rates::TimeRateFilterInv
  norm::Float64
end
Base.ndims(x::ScaleRateFilterInv) = 2
Base.inv(cf::ScaleRateFilter;norm=0.9) =
  ScaleRateFilterInv(inv(cf.scales),inv(cf.rates),norm)
AxisArrays.axisnames(x::ScaleRateFilterInv) =
  (axisname(x.scales),axisname(x.rates))

function list_filters(z_cum,cr,cf::ScaleRateFilterInv)
  ((CartesianIndex(i,j), HR.*HS)
   for (i,HR) in list_filters(z_cum,cr,cf.rates)
   for (j,HS) in list_filters(z_cum,cr,cf.scales))
end

################################################################################
# private helper functions

function find_fft_dims(y)
  @assert axisdim(y,Axis{:freq}) == ndims(y)
  @assert axisdim(y,Axis{:time}) == 1
  find_fft_dims(size(y))
end
find_fft_dims(y::NTuple{N,Int}) where {N} =
  (nextprod([2,3,5],y[1]),y[2:end-1]...,nextprod([2,3,5],y[end]))

struct FIRFiltering{T,N,P}
  Y::Array{T,N}
  plan::P
end

function FIRFiltering(y,axis)
  dims = map(AxisArrays.axes(y)) do ax
    if AxisArrays.axes(y,axis) == ax
      2nextprod([2,3,5],length(ax))
    else
      length(ax)
    end
  end

  along = axisdim(y,axis)
  Y = fft(pad(y,dims),along)
  FIRFiltering(Y,plan_ifft(Y,along))
end
apply(fir::FIRFiltering,H) = fir.plan * (fir.Y .* H)
Base.size(x::FIRFiltering,i...) = size(x.Y,i...)
Base.ndims(x::FIRFiltering) = ndims(x.Y)

function initfilter(y,ratefilter::TimeRateFilter)
  bandonly = ratefilter.bandonly
  rates = sort(ratefilter.data)
  r = Axis{axisname(ratefilter)}(rates)
  ax = AxisArrays.axes(y)
  newax = ax[1],r,ax[2:end]...

  axar = AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
  arates = sort!(unique!(abs.(rates)))
  rate_axis = bandonly ? 
     RateAxis(-Inf*Hz,Inf*Hz) : RateAxis(first(arates),last(arates))
  axis_meta = addaxes(getmeta(y);Dict(axisname(ratefilter) => rate_axis)...)
  MetaAxisArray(axis_meta,axar)
end

function initfilter(y,scalefilter::FreqScaleFilter)
  bandonly = scalefilter.bandonly
  scales = sort(scalefilter.data)
  s = Axis{axisname(scalefilter)}(scales)
  ax = AxisArrays.axes(y)
  newax = ax[1:end-1]...,s,ax[end]

  axar = AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
  scale_axis = bandonly ? 
    ScaleAxis(-Inf*cycoct,Inf*cycoct) :
    ScaleAxis(first(scales),last(scales))
  axis_meta = addaxes(getmeta(y);Dict(axisname(scalefilter) => scale_axis)...)
  MetaAxisArray(axis_meta,axar)
end

# TODO: do this for rates as well
reshape_for(v::Array{T,3},cr::AxisArray{T,3}) where T = v
reshape_for(v::Array{T,4},cr::AxisArray{T,4}) where T = v
reshape_for(v::Array{T,3},cr::AxisArray{T,4}) where T =
    reshape(v,ntimes(cr),1,nfrequencies(cr))

# keeps track of cumulative sum of FIR filters
# in frequency-space so we can readily normalize the result.
struct FFTCum{T,P,N,M}
  z::Array{Complex{T},N}
  z_cum::Array{Complex{T},M}
  h_cum::Array{T,M}
  nfrequencies::Int
  ntimes::Int
  plan::P
end

# TODO: working on generalizing FFTCum to working
# for any number of scale and rate axes
function withoutdim(dims,without)
  dims[1:without-1]...,dims[without+1:end]...
end

function FFTCum(cr::MetaAxisArray,withoutaxes)
  dims = find_fft_dims((size(cr,1),size(cr,ndims(cr))))
  mult = fill(ndims(cr),1)
  mult[1] += any(ax -> startswith(string(ax),"scale"),withoutaxes)
  mult[end] += any(ax -> startswith(string(ax),"rate"),withoutaxes)
  z = zeros(eltype(cr),(dims .* mult)...)

  cumsize = (size(z,1),size(z,2))
  z_cum = zeros(eltype(z),cumsize)
  h_cum = zeros(real(eltype(z)),cumsize)
  plan = plan_fft(z)
  FFTCum(z,z_cum,h_cum,nfrequencies(cr),ntimes(cr),plan)
end

Base.size(x::FFTCum,i...) = size(x.z_cum,i...)
Base.ndims(x::FFTCum) = ndims(x.z_cum)

function addfft!(x::FFTCum,cr,h)
  @assert x.ntimes == size(cr,1)
  @assert x.nfrequencies == size(cr,2)

  x.z[1:x.ntimes,1:x.nfrequencies] = cr
  Z = x.plan * x.z
  x.h_cum .+= abs2.(h)
  x.z_cum .+= h .* Z

  x
end

function normalize!(x::FFTCum,cr,norm)
  x.h_cum[:,1] .*= 2
  old_sum = sum(x.h_cum[:,nfrequencies(cr)])
  x.h_cum .= norm.*x.h_cum .+ (1 .- norm).*maximum(x.h_cum)
  x.h_cum .*= old_sum ./ sum(view(x.h_cum,:,nfrequencies(cr)))
  x.z_cum ./= x.h_cum

  spectc = view((x.plan \ x.z_cum),1:ntimes(cr),1:nfrequencies(cr))
  max.(real.(2 .* spectc),0)
end

pad(x,lens) = pad(x,lens...)
function pad(x,lens::T...) where T <: Number
  @assert all(size(x) .<= lens)
  y = zeros(eltype(x),lens)
  y[Base.axes(x)...] = x
  y
end

# transforms a bandpass frequency response into either a high or low pass
# response (or leaves it untouched)
function askind(H,len,maxi,kind,nonorm)
  if kind == :band
    H
  else
    old_sum = sum(H)
    if kind == :low
      H[1:maxi-1] .= 1
    elseif kind == :high
      H[maxi+1:len] .= 1
    else
      error("Unexpected filter kind '$kind'.")
    end
    if !nonorm
      H .= H ./ sum(H) .* old_sum
    end

    H
  end
end

function scale_filters(Y,x,scaleax)
  N_f = size(Y,ndims(Y)) >> 1
  scaleparam = getproperty(x,scaleax)
  map(scales(x)) do scale
	  scale_filter(ustrip(uconvert(cycoct,scale)), N_f, spect_rate,
                 scale <= scaleparam.low ? :low : 
                 scale < scaleparam.high ? :band : :high)
  end
end

# create the frequency-scale filter (filter along spectral axis)
function scale_filter(scale,len,ts,kind)
  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)

  askind(H,len,argmax(H),kind,false)
end

function rate_filters(Y,x,rateax;use_conj=false)
  N_t = size(Y,1) >> 1
  rateparam = getproperty(x,rateax)

  map(rates(x)) do rate
    rate_filter(ustrip(uconvert(Hz,rate)), N_t, x.time.Δ,
                abs(rate) <= rateparam.low ? :low :
                abs(rate) < rateparam.high ? :band : :high,use_conj)
  end
end

# create the temporal-rate filter (filter along temporal axis)
function rate_filter(rate,len,Δt,kind,use_conj=false,return_partial=false)
  t = (0:len-1)*ustrip(uconvert(s,Δt))*abs(rate)
  h = @. sin(2π*t) * t^2 * exp(-3.5t)
  h .-= mean(h)

  H0 = view(fft(pad(h,2len)),1:len)
  A = angle.(H0)
  H = abs.(H0)

  maxH,maxi = findmax(H)
  H ./= maxH
  HR = askind(H,len,maxi,kind,true) .* exp.(A*im)

  if use_conj
    HR = conj.(HR)
  end

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  if return_partial
    HR,h
  else
    HR
  end
end
