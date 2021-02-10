using AxisArrays
using MetaArrays
using ImageInTerminal
using Colors

struct AxisMeta{Axs}
    axes::Axs # map from :axis -> metadata
end
AxisMeta(;kwds...) = AxisMeta(kwds.data)
Base.getproperty(x::AxisMeta,field::Symbol) =
  getproperty(getfield(x,:axes),field)

function validate(meta::AxisMeta,data::AxisArray)
    for field in fieldnames(typeof(Base.getfield(meta,:axes)))
        if field ∉ axisnames(data)
            error("Extra field `$field` in axis metadata.")
        end
    end

    for axis in axisnames(data)
        if axis ∉ fieldnames(typeof(Base.getfield(meta,:axes)))
            error("Missing field `$axis` in axis metadata.")
        end
    end
end

axismeta(data::AxisArray;keys...) = axismeta(keys.data,data)
axismeta(meta::NamedTuple,data::AxisArray) = axismeta(AxisMeta(meta),data)
function axismeta(meta::AxisMeta,data::AxisArray)
    validate(meta,data)
    MetaArray(meta,data)
end

const MetaAxisArray{A,M} = MetaArray{A,M} where {A<:AxisArray,M<:AxisMeta}
const MetaAxisLike = Union{AxisMeta,MetaAxisArray}
function MetaAxisArray(meta::AxisMeta,data::AxisArray)
    validate(meta,data)
    MetaArray(meta,data)
end

function addaxes(meta::AxisMeta;kwds...)
    AxisMeta(merge(Base.getfield(meta,:axes),kwds.data))
end

function removeaxes(meta::AxisMeta,rem...)
  axes = Base.getfield(meta,:axes)
  fields = setdiff(fieldnames(typeof(axes)),rem)
  vals = Tuple(map(f -> axes[f],fields))
  AxisMeta(NamedTuple{Tuple(fields),typeof(vals)}(vals))
end

function describe_axes(io::IO,x)
  for ax in AxisArrays.axes(x)
    println(io,string(AxisArrays.axisname(ax))," ",
            string(round(ustrip(ax.val[1]),digits=2)),
            " - ",string(round(ustrip(ax.val[end]),digits=2)),
            string(unit(ax.val[1])))
  end
end

resultname(x) = "Data with axes: "
function Base.show(io::IO,mime::MIME"text/plain",x::MetaAxisArray)
  if get(io, :compact, false)
    println(io,resultname(x))
  else
    println(io,resultname(x))
    describe_axes(io,x)
    if ndims(x) <= 2
      if eltype(x) <: Complex
        @warn "Ignoring phase in display of complex values."
        x = real.(abs.(x))
      end
      lo,hi = extrema(x)
      x .= (x .- lo) ./ (hi - lo)
      if ndims(x) == 2
        image = reverse(Gray.(getcontents(x)),dims=2)'
      else
        image = Gray.(getcontents(x))
      end
      ImageInTerminal.imshow(io,image,ImageInTerminal.colormode[1])
    end
  end
end

Base.show(io::IO, m::MIME"img/png", x::MetaAxisArray) = show(io,m,plotaxes(x))