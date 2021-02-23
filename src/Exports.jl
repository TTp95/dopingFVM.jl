"""
Export file...
"""

macro publish(mod,name)
  quote
    using dopingFVM.$mod: $name; export $name
  end
end

@publish Structures M
