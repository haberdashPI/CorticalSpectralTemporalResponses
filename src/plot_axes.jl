export logrange
using Compat

function PlotAxes.asplotable(x::AuditorySpectrogram,args...;quantize=(100,128),
    kwds...)
  args = replace(collect(args),:freq => (:freq => logrange))

  asplotable(getcontents(x),args...;quantize=quantize,kwds...)
end

default_quantize(x) = (100,)
default_quantize(x,y) = (100,128,)
default_quantize(x,y,args...) = (100,128,fill(10,length(args))...)

function PlotAxes.asplotable(x::MetaAxisArray;kwds...) 
  asplotable(x,axisnames(x)...;kwds...)
end
function PlotAxes.asplotable(x::MetaAxisArray,ax1,axes...;
  quantize=default_quantize(ax1,axes...),kwds...)
  # TODO: handle time or freq not being present
  args = (ax1,axes...)
  args = replace(collect(args),:freq => (:freq => logrange))
  infront = filter(!isnothing,indexin([:time,:freq => logrange],args))
  others = setdiff(1:length(args),infront)
  args = args[vcat(infront,others)]

  asplotable(getcontents(x),args...;quantize=quantize,kwds...)
end

            
            
