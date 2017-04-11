 #!/bin/env ruby

lines = IO.readlines(ARGV.shift)

header = lines.shift.split(' ')

species = "unknown_species"

lines.each do |line|
	elements = line.strip.split(' ')
        hit_rate = elements[4].to_f
        hits = elements[1].to_i
        this_organism = elements[0]
        if hits > 500000
        	species = this_organism
        end
end

puts species


