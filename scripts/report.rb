 #!/bin/env ruby

lines = IO.readlines(ARGV.shift)

header = lines.shift.split(' ')

species = "unknown_species"

bucket = {}

lines.each do |line|
	elements = line.strip.split(' ')
        hit_rate = elements[4].to_f
        hits = elements[1].to_i
        this_organism = elements[0]
      	bucket[this_organism] = hit_rate
end

current_max = 0.0
bucket.each do |organism,rate|
	if rate > current_max && rate > 0.3
		unless organism == "noMatch" && current_max > 0.3
			species = organism		
			current_max = rate
		end
	end
end

puts species


